# -*- coding: utf-8 -*-
"""
Groups the different stages of the processing described in the article :
"Bernhard Jenny (2021) Terrain generalization with line integral
convolution, Cartography and Geographic Information Science, 48:1, 78-92, DOI:
10.1080/15230406.2020.1833762
https://doi.org/10.1080/15230406.2020.1833762"

Including the preliminary processing of adaptive Gaussian blur applied to the DTM.
"""


# Extern imports
import numpy as np
import rasterio
from rasterio.windows import Window
from tqdm import tqdm
from scipy.ndimage import gaussian_filter
from functools import partial

# intern imports
from dem_lic.utils.morpho_dem import (
    calculate_maximal_curvature,
    fast_adaptive_gaussian_blur,
    initialize_flat_steep_grid,
    remove_small_flat_areas,
)




def extended_lic_weighted_altitude_lengthModulated(
    grid: np.ndarray,
    cellsize: float,
    f: np.ndarray,
    num_steps: int,
    sigma_modulated: bool = True
) -> np.ndarray:
    """Applies a Line Integral Convolution (LIC) with a dynamically modulated integration
    length (based on f) and a variable Gaussian kernel that reduces the influence of
    higher terrain along the line of integration.
    
    Parameters
    ----------
    grid : np.ndarray
        A 2D array representing the Digital Elevation Model (DEM).
    cellsize : float
        The size of each raster cell in meters.
    f : np.ndarray
        A 2D array (same shape as `grid`) in [0, 1], controlling the fraction of 
        the maximum integration steps allowed for each pixel (e.g., in flat or steep areas).
    num_steps : int
        The maximum number of integration steps for each pixel.
    sigma_modulated : bool, optional
        If True, the sigma for the weighting in the LIC algorithm is modulated by the altitude.
        If False, a fixed Gaussian sigma is used.
        
    Returns
    -------
    np.ndarray
        The DEM after the weighted LIC processing, where ridges are accentuated
        by decreasing the Gaussian kernel's sigma for portions that lie above
        the starting pixel's altitude.
    """

    # --- 0) Initial copies and shapes
    src = grid.copy()
    H, W = src.shape

    # The altitude of origin for each pixel (to compare with visited altitudes).
    start_alt = src.copy()

    # --- 1) Compute the local number of steps for each pixel from f
    f = np.clip(f, 0.0, 1.0)
    local_steps = np.round(f * num_steps).astype(int)
    max_steps = np.max(local_steps)

    # --- 2) Compute normalized gradients (vi, vj)
    DS = 0.5  # Integration speed (distance in pixel units per step)
    vi, vj = np.gradient(src, cellsize)  # Partial derivatives

    # Magnitude (avoid division by zero)
    mag = np.hypot(vi, vj)
    zero_mask = (mag == 0)
    mag[zero_mask] = 1.0
    mag *= (1.0 / DS)

    # Normalize gradients
    vi /= mag
    vj /= mag

    # --- 3) Pass 1: Compute range along each integration line
    range_altitude = np.zeros_like(src, dtype=float)

    for pass_id in range(2):
        if pass_id == 1:
            vi = -vi
            vj = -vj

        current_steps = local_steps.copy()
        active_cells = np.ones((H, W), dtype=bool)

        fi, fj = np.mgrid[:H, :W].astype(float)
        min_alt = start_alt.copy()
        max_alt = start_alt.copy()

        for step_idx in range(max_steps):
            if not np.any(active_cells):
                break

            fi[active_cells] += vi[active_cells]
            fj[active_cells] += vj[active_cells]

            pi = np.clip(fi.astype(int), 0, H - 1)
            pj = np.clip(fj.astype(int), 0, W - 1)

            # Update min and max altitude along the line
            min_alt = np.minimum(min_alt, src[pi, pj])
            max_alt = np.maximum(max_alt, src[pi, pj])

            # Update range altitude for each pixel
            range_altitude = np.maximum(range_altitude, max_alt - min_alt)

            current_steps[active_cells] -= 1
            active_cells = (current_steps > 0)

    range_altitude = np.maximum(range_altitude, 1e-6)  # Avoid division by zero

    # --- 4) Pass 2: Apply LIC with range-based modulation
    result = np.zeros_like(src, dtype=float)
    weights_sum = np.zeros_like(src, dtype=float)

    sigma_base = np.sqrt((num_steps**2 - 1) / 12.0)

    for pass_id in range(2):
        if pass_id == 1:
            vi = -vi
            vj = -vj

        current_steps = local_steps.copy()
        active_cells = np.ones((H, W), dtype=bool)

        fi, fj = np.mgrid[:H, :W].astype(float)

        for step_idx in range(max_steps):
            if not np.any(active_cells):
                break

            fi[active_cells] += vi[active_cells]
            fj[active_cells] += vj[active_cells]

            pi = np.clip(fi.astype(int), 0, H - 1)
            pj = np.clip(fj.astype(int), 0, W - 1)

            dh = src[pi, pj] - start_alt  # Compute altitude difference

            # Normalize dh by range along the integration line
            if sigma_modulated:
                f_alt = 1.0 - (dh / range_altitude)
                f_alt = np.clip(f_alt, 0.0, 1.0)
            else:
                f_alt = np.ones_like(dh)

            # Adjust sigma accordingly
            sigma_adjusted = np.maximum(sigma_base * f_alt, 1e-6)

            # Compute Gaussian weight
            distance = step_idx * DS
            weight = np.exp(-(distance**2) / (2 * sigma_adjusted**2))

            # Accumulate weighted values
            result[active_cells] += weight[active_cells] * src[pi[active_cells], pj[active_cells]]
            weights_sum[active_cells] += weight[active_cells]

            current_steps[active_cells] -= 1
            active_cells = (current_steps > 0)

    # --- 5) Normalize final result
    np.divide(result, np.maximum(weights_sum, 1e-12), out=result)

    return result




def LIC_iterations(
    grid: np.ndarray,
    # local_range_altitude: np.ndarray,
    cellsize: float,
    f: np.ndarray,
    num_steps: int,
    sigma_modulated: bool,
    n_iterations: int,
    sigma: float,
    k: float,
) -> np.ndarray:
    """Perform multiple iterations of the Line Integral Convolution (LIC) method to generalize a DEM
    by combining smoothing and feature enhancement.

    Parameters
    ----------
    grid : numpy.ndarray
        2D array representing the input Digital Elevation Model (DEM).
    # local_range_altitude : numpy ndarray
    #     2D array of the range of elevation in a windows centered on each point
    cellsize : float
        The spatial resolution of the DEM (in meters).
    f : numpy.ndarray
        Continuous grid used to modulate the integration length.
    num_steps : int
        Maximum length of the integration line.
    sigma_modulated : bool, optional
        If True, the sigma for the weighting in the LIC algorithm is modulated by the altitude.
        If False, a fixed Gaussian sigma is used.
    n_iterations : int
        Number of iterations of the LIC algorithm.
    sigma : float
        Standard deviation for the Gaussian blur applied to the curvature.
    k : float
        Weighting factor for the combination of grids.

    Returns
    -------
    numpy.ndarray
        2D array representing the DEM after LIC iterations.

    Notes
    -----
    This function performs the following steps for each iteration:
    1. Computes a modulated LIC 
    2. Calculates the maximal curvature of the resulting grid.
    3. Applies Gaussian smoothing to the curvature.
    4. Updates the DEM using a weighted combination of the original DEM and the LIC result.
    """
    for n in range(n_iterations):
        print(f"Iteration {n+1}/{n_iterations}")

        # Step 1: Compute a modulated LIC
        filtered_grid = extended_lic_weighted_altitude_lengthModulated(
            grid, cellsize, f, num_steps, sigma_modulated
        )

        # Step 2: Calculate maximal curvature of the filtered grid
        filtered_grid_max_c = calculate_maximal_curvature(
            filtered_grid, cellsize
        )

        # Step 3: Smooth the curvature using Gaussian filtering
        filtered_grid_max_c_blured = gaussian_filter(
            filtered_grid_max_c, sigma=sigma, mode="nearest"
        )

        # Step 4: Compute weighting factor from the smoothed curvature
        weight = k * np.abs(filtered_grid_max_c_blured)
        weight = np.nan_to_num(weight, nan=0.0, posinf=1.0, neginf=0.0)

        # Step 5: Update the DEM using the weighted combination
        grid = weight * grid + (1 - weight) * filtered_grid

    return grid


def correct_flat_area_values(
    b: np.ndarray,
    c: np.ndarray,
    sigma: float,
    epsilon: float = 1e-6,
) -> np.ndarray:
    """Adjusts the values in flat areas of a continuous grid after Gaussian blurring
    to improve the representation of transitions.

    Parameters
    ----------
    b : numpy.ndarray
        2D binary array representing flat areas (0) and steep areas (1).
    c : numpy.ndarray
        2D continuous array, the result of applying Gaussian blur to the binary grid.
    sigma : float
        Gaussian blur parameter for processing the differences.
    epsilon : float, optional
        A small value to avoid division by zero, by default 1e-6.

    Returns
    -------
    numpy.ndarray
        2D corrected array after adjusting values in flat areas.

    Notes
    -----
    The function modifies the continuous grid `c` by correcting values in flat areas (indicated by 0 in `b`).
    The adjustment is based on a blurred version of the difference between `b` and `c`, scaled by a
    computed factor `s`. This ensures smoother transitions while preserving flat and steep area characteristics.
    """
    # Step 1: Compute the difference between the binary and continuous grids
    diff = b - c

    # Step 2: Apply Gaussian blur to the difference
    blurred_diff = gaussian_filter(diff, sigma=sigma, mode="nearest")

    # Step 3: Calculate the scaling factor `s`
    max_diff = np.max(diff)
    max_blurred_diff = np.max(blurred_diff)
    s = max_diff / (max_blurred_diff + epsilon)  # Avoid division by zero

    # Step 4: Compute the corrected grid
    f = c + s * blurred_diff

    return f


def process_geotiff_in_block_with_overlap(
    MNT_input_path: str,
    output_path: str,
    processing_function,
    block_size: int = 2000,
    overlap: int = 20,
    **kwargs
) -> None:
    """Processes a GeoTIFF file in overlapping blocks and applies a user-defined processing function.

    Parameters
    ----------
    MNT_input_path : str
        Path to the input GeoTIFF file containing the Digital Elevation Model (DEM).
    output_path : str
        Path to the output GeoTIFF file.
    processing_function : function
        A function that processes a single block. It should accept a block (ndarray), cell size, 
        and any additional parameters.
    block_size : int, optional
        Size of non-overlapping blocks (in pixels) to process. Default is 2000.
    overlap : int, optional
        Size of the overlapping region (in pixels) between adjacent blocks. Default is 20.
    **kwargs
        Additional arguments passed to `processing_function`.

    Returns
    -------
    None
        The processed output is written directly to the specified GeoTIFF file.
    """

    # Open the input GeoTIFF and read metadata
    with rasterio.open(MNT_input_path) as MNT_src:
        profile = MNT_src.profile.copy()
        width, height = MNT_src.width, MNT_src.height
        cellsize = MNT_src.transform[0]

        # Update output profile
        profile.update(dtype=np.float32, compress="lzw", bigtiff="YES")

        print("Starting terrain generalization...")

        # Compute total number of blocks for progress tracking
        n_blocks_x = (width + block_size - 1) // block_size
        n_blocks_y = (height + block_size - 1) // block_size
        total_blocks = n_blocks_x * n_blocks_y

        # Create output GeoTIFF
        with rasterio.open(output_path, "w", **profile) as dst:
            with tqdm(total=total_blocks, desc="Processing blocks", unit="block") as pbar:
                for y in range(0, height, block_size):
                    for x in range(0, width, block_size):
                        # Define window with overlap
                        x_min, x_max = max(0, x - overlap), min(width, x + block_size + overlap)
                        y_min, y_max = max(0, y - overlap), min(height, y + block_size + overlap)

                        window = Window(x_min, y_min, x_max - x_min, y_max - y_min)

                        # Read the MNT block
                        MNT_block = MNT_src.read(1, window=window)


                        # Apply the processing function to the block
                        processed_block = processing_function(MNT_block, cellsize, **kwargs)

                        # Extract the central block (remove overlap)
                        if x == 0 and y != 0:
                            central_block = processed_block[overlap:block_size + overlap, 0:block_size]
                        elif y == 0 and x != 0:
                            central_block = processed_block[0:block_size, overlap:block_size + overlap]
                        elif y == 0 and x == 0:
                            central_block = processed_block[0:block_size, 0:block_size]
                        else:
                            central_block = processed_block[overlap:block_size + overlap, overlap:block_size + overlap]

                        # Define output window
                        output_window = Window(x, y, central_block.shape[1], central_block.shape[0])

                        # Write processed block to output file
                        dst.write(central_block.astype(np.float32), 1, window=output_window)

                        # Update progress bar
                        pbar.update(1)

    print(f"Processing completed. Output written to {output_path}")




def process_lic_extended(
    MNT_block: np.ndarray,
    cellsize: float,
    sigma_max: float = 5.0,
    slope_threshold: float = 6,
    num_bins: int = 10,
    min_area: int = 100,
    num_steps: int = 5,
    sigma_modulated: bool = True,
    n_iterations: int = 5,
    sigma_blur_maxcurv: float = 3.0,
    k: float = 2.5,    
) -> np.ndarray:
    """Processes a single block of a DEM with LIC and generalization techniques.

    Parameters
    ----------
    MNT_block : np.ndarray
        The block of the Digital Elevation Model.
    cellsize : float
        Size of the raster cells in meters.
    sigma_max : float
        Maximum kernel width for the adaptive Gaussian blur.
    slope_threshold : float
        Threshold for slope to distinguish flat and steep areas.
    num_bins : int
        Number of bins for sigma approximation in the adaptive blur.
    min_area : int
        Minimum size (in pixels) for flat areas to be preserved.
    num_steps : int
        Maximum integration length for the LIC algorithm.
    sigma_modulated : bool
        If True, sigma is modulated by altitude in LIC processing.
    n_iterations : int
        Number of iterations for the LIC algorithm.
    sigma_blur_maxcurv : float
        Gaussian blur parameter for maximum curvature.
    k : float
        Weighting factor for combining LIC grids.

    Returns
    -------
    np.ndarray
        Processed DEM block.
    """

    # Compute maximum curvature
    curvature_block = calculate_maximal_curvature(MNT_block, cellsize)

    # Apply adaptive Gaussian blur
    MNT_blurred_block = fast_adaptive_gaussian_blur(
        MNT_block, curvature_block, cellsize, sigma_max, num_bins
    )

    # Generate binary grid for flat/steep areas
    flat_steep_grid_block = initialize_flat_steep_grid(MNT_blurred_block, slope_threshold, cellsize)

    # Remove small flat areas
    cleaned_grid_block = remove_small_flat_areas(flat_steep_grid_block, min_area)

    # Smooth the binary grid
    continuous_grid_block = gaussian_filter(cleaned_grid_block.astype(float), sigma=5.0, mode="nearest")

    # Correct flat area values
    corrected_f = correct_flat_area_values(cleaned_grid_block, continuous_grid_block, sigma=2.0)

    # Perform multi-pass LIC processing
    processed_block = LIC_iterations(
        MNT_blurred_block,
        cellsize,
        corrected_f,
        num_steps,
        sigma_modulated,
        n_iterations,
        sigma_blur_maxcurv,
        k,
    )

    return processed_block


lic_extended_partial = partial(
    process_lic_extended,
    sigma_max=5.0,
    slope_threshold=6,
    num_bins=10,
    min_area=100,
    num_steps=5,
    sigma_modulated=True,
    n_iterations=5,
    sigma_blur_maxcurv=3.0,
    k=2.5
)


if __name__ == "__main__":
    """Main script for processing a GeoTIFF file with Line Integral Convolution (LIC) and
    adaptive Gaussian blur techniques for terrain generalization.

    """


    # File paths
    MNT_input_path = ""
    LIC_complete_output_path = ""
    
    
    # Call the processing function
    process_geotiff_in_block_with_overlap(
        MNT_input_path,
        LIC_complete_output_path,
        processing_function=lic_extended_partial,
        sigma_max = 2.0,
        num_steps = 5,
        sigma_modulated= True,
        n_iterations = 4,
        sigma_blur_maxcurv = 2.9,
    )

