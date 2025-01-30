# -*- coding: utf-8 -*-
"""
Groups the different stages of the processing described in the article :
"Bernhard Jenny (2021) Terrain generalization with line integral
convolution, Cartography and Geographic Information Science, 48:1, 78-92, DOI:
10.1080/15230406.2020.1833762
https://doi.org/10.1080/15230406.2020.1833762"

Including the preliminary processing of adaptive Gaussian blur applied to the DTM.
Use of a raster of relative altitudes
"""


# Extern imports
import numpy as np
import rasterio
from rasterio.windows import Window
from tqdm import tqdm
from scipy.ndimage import gaussian_filter
from typing import Dict

# intern imports
from dem_lic.utils.morpho_dem import (
    calculate_maximal_curvature,
    fast_adaptive_gaussian_blur,
    calculate_local_range,
    initialize_flat_steep_grid,
    remove_small_flat_areas,
    calculate_relative_altitude,
)


# def extended_lic_weighted_altitude(
#     grid: np.ndarray,
#     local_range_altitude: np.ndarray,
#     cellsize: float,
#     num_steps: int,
# ) -> np.ndarray:
#     """Applies a weighted Line Integral Convolution (LIC) on a DEM, using a Gaussian kernel 
#     whose sigma is adjusted based on local altitude range and the difference in elevation 
#     from each pixel's starting point. (chapter 3.4).

#     Parameters
#     ----------
#     grid : numpy.ndarray
#         A 2D array representing the Digital Elevation Model (DEM).
#     local_range_altitude : np.ndarray
#         A 2D array (same shape as grid) giving the local elevation range around each cell.
#         Typically computed by applying a max/min filter with a certain window size.
#     cellsize : float
#         The size of the raster cells in meters.
#     num_steps : int
#         The maximum length of the integration line.

#     Returns
#     -------
#     numpy.ndarray
#         A 2D array resulting from the LIC processing weighted by Gaussian adjustment.

#     Notes
#     -----
#     In this method, each pixel accumulates elevation values along the direction of maximum 
#     slope (gradient), in both forward and backward passes. The Gaussian weighting is adapted 
#     so that if the encountered altitude is significantly higher than the pixel's original 
#     altitude, the influence is reduced, thereby accentuating ridges.
#     """

#     src = grid.copy()
#     start_alt = src.copy() #elevation reference
#     sink = (slice(0, num_steps // 2), slice(0, num_steps // 2))
#     src0 = src[sink].copy()  # Store initial border values
#     H, W = src.shape

#     DS = 0.5  # Integration speed

#     # Step 1: Compute normalized gradients
#     vi, vj = np.gradient(src, cellsize)  # Partial derivatives in Y and X directions
#     res = np.hypot(vi, vj)  # Gradient magnitude
#     m = np.equal(res, 0)  # Mask to avoid division by zero
#     np.copyto(res, 1, where=m)  # Replace zeros with 1
#     np.multiply(res, 1.0 / DS, out=res)
#     for grad in (vi, vj):
#         np.divide(grad, res, out=grad)  # Normalize gradients

#     m2 = m.copy()
#     weights_sum = np.zeros_like(src, dtype=float)  # Sum of weights for normalization
#     res = np.zeros_like(src, dtype=float)  # Initialize result grid

#     # Step 2: Compute base sigma
#     l = num_steps  # Integration line length
#     sigma_base = np.sqrt((l**2 - 1) / 12)


#     for i in range(2):
#         if i:  # Reverse gradients for backward integration
#             for grad in (vi, vj):
#                 np.negative(grad, out=grad)
#         pi, pj = np.mgrid[slice(H), slice(W)]  # Initialize grid coordinates
#         fi, fj = (
#             _.astype(float) for _ in (pi, pj)
#         )  # Float coordinates for interpolation
#         for step in range(num_steps):
#             fi += vi[pi, pj]  # Update float coordinates
#             fj += vj[pi, pj]
#             np.trunc(fi, out=pi, casting="unsafe")  # Convert to integer indices
#             np.trunc(fj, out=pj, casting="unsafe")

#             # Handle out-of-bounds indices
#             np.less(pi, 0, m)
#             for a, v in [(pi, H), (pj, 0), (pj, W)]:
#                 if v != 0:
#                     np.greater_equal(a, v, m2)
#                 else:
#                     np.less(a, v, m2)
#                 np.logical_or(m, m2, m)
#             np.copyto(pi, 0, where=m)  # Reset out-of-bounds indices to 0
#             np.copyto(pj, 0, where=m)
            
#             range_val = local_range_altitude[i, j] 
#             dh = src[pi, pj] - start_alt
#             f_alt = 1.0 - dh/local_range_altitude
#             sigma_adjusted = np.maximum(
#                 f_alt * sigma_base, 1e-6
#             )  # Ensure sigma is non-zero


#             # Compute Gaussian weight with adjusted sigma
#             distance = step * DS
#             weight = np.exp(-(distance**2) / (2 * sigma_adjusted**2))  # Gaussian weight
#             res += weight * src[pi, pj]  # Add weighted values
#             weights_sum += weight  # Accumulate weights for normalization

#             np.logical_not(m, m2)  # Update mask for valid indices

#     # Step 4: Normalize result
#     weights_sum[sink] = 1
#     res[sink] = src0  # Reinstate initial border values
#     np.divide(res, weights_sum, out=res)  # Normalize by sum of weights

#     return res


def extended_lic_weighted_altitude(
    grid: np.ndarray,
    local_range_altitude: np.ndarray,
    cellsize: float,
    num_steps: int,
) -> np.ndarray:
    """Applies a weighted Line Integral Convolution (LIC) on a DEM, using a Gaussian kernel 
    whose sigma is adjusted based on local altitude range and the difference in elevation 
    from each pixel's starting point.
    
    Parameters
    ----------
    grid : np.ndarray
        A 2D array representing the Digital Elevation Model (DEM).
    local_range_altitude : np.ndarray
        A 2D array (same shape as grid) giving the local elevation range around each cell.
        Typically computed by applying a max/min filter with a certain window size.
    cellsize : float
        The size of the raster cells in meters.
    num_steps : int
        The maximum length of the integration line (number of steps in each direction).
    
    Returns
    -------
    np.ndarray
        The resulting 2D array after the LIC processing, weighted by a variable Gaussian kernel 
        that depends on local elevation range and the altitude difference along each flow line.
    
    Notes
    -----
    In this method, each pixel accumulates elevation values along the direction of maximum 
    slope (gradient), in both forward and backward passes. The Gaussian weighting is adapted 
    so that if the encountered altitude is significantly higher than the pixel's original 
    altitude, the influence is reduced, thereby accentuating ridges.
    """

    src = grid.copy()
    start_alt = src.copy()  # reference altitude for each pixel
    H, W = src.shape

    # Integration speed (step in pixels)
    DS = 0.5

    # 1. Normalized gradients
    vi, vj = np.gradient(src, cellsize)
    mag = np.hypot(vi, vj)
    mask_zero = (mag == 0)
    mag[mask_zero] = 1  # avoid division by zero
    mag *= (1.0 / DS)
    vi /= mag
    vj /= mag

    # 2. Initialize accumulators
    weights_sum = np.zeros_like(src, dtype=float)
    res = np.zeros_like(src, dtype=float)

    # Base sigma for the Gaussian, e.g. from a formula
    l = num_steps
    sigma_base = np.sqrt((l**2 - 1) / 12)

    # We'll create a 2D indexing grid for i, j
    i_coords, j_coords = np.mgrid[:H, :W]

    for direction in range(2):
        # Reverse gradient for the 2nd pass
        if direction == 1:
            vi = -vi
            vj = -vj

        # Current position arrays
        pi = i_coords.copy().astype(float)
        pj = j_coords.copy().astype(float)

        for step in range(num_steps):
            # Move in the direction of the slope
            pi += vi[i_coords, j_coords]
            pj += vj[i_coords, j_coords]

            # Truncate to integer indices
            np.trunc(pi, out=pi)
            np.trunc(pj, out=pj)

            # Clip out-of-bounds
            pi_clipped = np.clip(pi, 0, H-1).astype(int)
            pj_clipped = np.clip(pj, 0, W-1).astype(int)

            # Compute delta h relative to the starting altitude
            # shape (H, W)
            dh = src[pi_clipped, pj_clipped] - start_alt[i_coords, j_coords]

            # Retrieve the local altitude range for each pixel of origin
            range_val = local_range_altitude[i_coords, j_coords]
            range_val = np.maximum(range_val, 1e-6)  # avoid div zero

            # Factor that decreases if the visited altitude is higher
            f_alt = 1.0 - (dh / range_val)
            # Clip between 0 and 1 (for example)
            f_alt = np.clip(f_alt, 0.0, 1.0)

            sigma_adjusted = np.maximum(f_alt * sigma_base, 1e-6)

            # Distance from center in "pixel steps"
            distance = step * DS

            # Gaussian weight
            weight = np.exp(-(distance**2) / (2.0 * sigma_adjusted**2))

            res += weight * src[pi_clipped, pj_clipped]
            weights_sum += weight

    # Final normalization
    np.divide(res, weights_sum, out=res, where=(weights_sum != 0))

    return res


# def extended_lic_weighted_altitude_lengthModulated(
#     grid: np.ndarray,
#     relative_altitude: np.ndarray,
#     cellsize: float,
#     f: np.ndarray,
#     num_steps: int,
# ) -> np.ndarray:
#     """Applies a Line Integral Convolution (LIC) with integration line length dynamically modulated by a mask.
#     Chapter 3.6.

#     Parameters
#     ----------
#     grid : numpy.ndarray
#         A 2D array representing the Digital Elevation Model (DEM).
#     relative_altitude : numpy.ndarray
#         A 2D array of relative altitudes normalized between 0 and 1.
#     cellsize : float
#         The size of the raster cells in meters.
#     f : numpy.ndarray
#         Continuous modulation grid with values ranging between 0 and 1.
#     num_steps : int
#         The maximum length of the integration line.

#     Returns
#     -------
#     numpy.ndarray
#         A 2D array representing the smoothed DEM after LIC processing.

#     Notes
#     -----
#     This function dynamically adjusts the integration line length based on the modulation mask `f`,
#     allowing for finer control of the convolution process in flat or steep areas.
#     """

#     src = grid.copy()
#     H, W = src.shape

#     # Limit values of `f` between 0 and 1
#     f = np.clip(f, 0, 1)

#     # Precompute local integration lengths
#     local_steps = np.round(f * num_steps).astype(int)
#     max_steps = np.max(local_steps)  # Maximum number of steps
#     active_cells = np.ones(src.shape, dtype=bool)  # Initial mask for active cells

#     # Step 1: Compute normalized gradients
#     DS = 0.5  # Integration speed
#     vi, vj = np.gradient(src, cellsize)  # Partial derivatives
#     res = np.hypot(vi, vj)  # Gradient magnitude
#     m = np.equal(res, 0)  # Mask to avoid division by zero
#     np.copyto(res, 1, where=m)  # Replace zeros with 1
#     np.multiply(res, 1.0 / DS, out=res)
#     for grad in (vi, vj):
#         np.divide(grad, res, out=grad)  # Normalize gradients

#     # Initialize results
#     fi, fj = np.mgrid[slice(H), slice(W)].astype(float)  # Floating-point coordinates
#     weights_sum = np.zeros_like(src, dtype=float)  # Sum of weights for normalization
#     result = np.zeros_like(src, dtype=float)  # Initialize result array

#     # Invert the relative altitude grid (lower values become more significant)
#     relative_altitude = np.abs(relative_altitude - 1)

#     for i in range(2):
#         if i:  # Reverse gradients for backward integration
#             vi = -vi
#             vj = -vj

#         # Reset dynamic mask for this pass
#         active_cells = np.ones(src.shape, dtype=bool)

#         for step in range(max_steps):
#             if not np.any(active_cells):  # Stop if no active cells remain
#                 break

#             # Update coordinates for active cells
#             fi[active_cells] += vi[active_cells]
#             fj[active_cells] += vj[active_cells]

#             # Clip indices to image bounds
#             pi = np.clip(fi.astype(int), 0, H - 1)
#             pj = np.clip(fj.astype(int), 0, W - 1)

#             # # Use relative altitude to adjust weights
#             f_alt = relative_altitude[pi, pj]  # Local relative altitude values
            
#             # Calculate Gaussian weights
#             distance = step * DS
#             sigma_base = np.sqrt((num_steps**2 - 1) / 12)
#             sigma_adjusted = np.maximum(f_alt * sigma_base, 1e-6)
#             weight = np.exp(-(distance**2) / (2 * sigma_adjusted**2))

#             # Update results for active cells
#             result[active_cells] += (
#                 weight[active_cells] * src[pi[active_cells], pj[active_cells]]
#             )
#             weights_sum[active_cells] += weight[active_cells]

#             # Decrease remaining steps for active cells and update the active mask
#             local_steps[active_cells] -= 1
#             active_cells = local_steps > 0  # Update active cells

#     # Normalize results
#     np.divide(result, np.maximum(weights_sum, 1e-6), out=result)
#     return result



def extended_lic_weighted_altitude_lengthModulated(
    grid: np.ndarray,
    local_range_altitude: np.ndarray,
    cellsize: float,
    f: np.ndarray,
    num_steps: int,
) -> np.ndarray:
    """Applies a Line Integral Convolution (LIC) with a dynamically modulated integration
    length (based on f) and a variable Gaussian kernel that reduces the influence of
    higher terrain along the line of integration.
    
    Parameters
    ----------
    grid : np.ndarray
        A 2D array representing the Digital Elevation Model (DEM).
    local_range_altitude : np.ndarray
        A 2D array (same shape as `grid`) giving the local altitude range (e.g., 
        max-min in a neighborhood) for each pixel. Used to normalize the difference 
        in elevation along the integration line.
    cellsize : float
        The size of each raster cell in meters.
    f : np.ndarray
        A 2D array (same shape as `grid`) in [0, 1], controlling the fraction of 
        the maximum integration steps allowed for each pixel (e.g., in flat or steep areas).
    num_steps : int
        The maximum number of integration steps for each pixel.
    
    Returns
    -------
    np.ndarray
        The DEM after the weighted LIC processing, where ridges are accentuated
        by decreasing the Gaussian kernel's sigma for portions that lie above
        the starting pixel's altitude.
    
    Notes
    -----
    1. We perform two passes: forward and backward along the gradient vectors.
    2. The number of steps each pixel can take in each pass is `round(f[i,j] * num_steps)`.
    3. A variable Gaussian kernel is computed at each step. The difference 
       in altitude (dh) between the current position and the pixel's origin 
       is normalized by the local altitude range, and used to reduce sigma 
       for higher (amont) portions.
    4. The final values are normalized by the sum of the Gaussian weights.
    """


    # --- 0) Initial copies and shapes
    src = grid.copy()
    H, W = src.shape

    # The altitude of origin for each pixel (to compare with visited altitudes).
    # We'll call this "start_alt".
    start_alt = src.copy()

    # local_range_altitude should be > 0 to avoid division by zero
    # We'll ensure no zero-value (flat area) leads to a /0
    local_range_altitude = np.maximum(local_range_altitude, 1e-6)

    # --- 1) Compute the local number of steps for each pixel from f
    # Clip f between [0, 1]
    f = np.clip(f, 0.0, 1.0)
    # local_steps is how many steps each pixel can do, at most
    local_steps = np.round(f * num_steps).astype(int)
    max_steps = np.max(local_steps)

    # We'll keep track of "active" pixels that still have steps left
    # (we reset this mask for each forward/backward pass)
    active_cells = np.ones((H, W), dtype=bool)

    # --- 2) Compute normalized gradients (vi, vj)
    DS = 0.5  # Integration speed (distance in pixel units per step)
    vi, vj = np.gradient(src, cellsize)  # partial derivatives

    # Magnitude (avoid /0)
    mag = np.hypot(vi, vj)
    zero_mask = (mag == 0)
    mag[zero_mask] = 1.0
    mag *= (1.0 / DS)

    # Normalize
    vi /= mag
    vj /= mag

    # --- 3) Prepare accumulators
    result = np.zeros_like(src, dtype=float)    # sum of weighted alt
    weights_sum = np.zeros_like(src, dtype=float)  # sum of weights

    # We'll create coordinate grids for each pixel
    i_coords, j_coords = np.mgrid[:H, :W]  # shape (H, W)

    # We define a base sigma. You could choose e.g.:
    sigma_base = np.sqrt((num_steps**2 - 1) / 12.0)

    # --- 4) We'll do two passes: forward and backward
    for pass_id in range(2):
        # For the second pass, we invert the gradients
        if pass_id == 1:
            vi = -vi
            vj = -vj

        # We also reset the local steps array for this pass, 
        # because each pixel can do "f * num_steps" steps in forward 
        # and again in backward.
        current_steps = local_steps.copy()

        # Active mask resets
        active_cells[:] = True

        # We'll keep floating coords to track the position of each pixel
        fi = i_coords.astype(float)
        fj = j_coords.astype(float)

        # --- 4.1) For each step in [0..max_steps-1], 
        #          but we skip if no pixel can still move
        for step_idx in range(max_steps):
            if not np.any(active_cells):
                break

            # Update floating coords only for active pixels
            fi[active_cells] += vi[active_cells]
            fj[active_cells] += vj[active_cells]

            # Convert to integer indices (clip to bounds)
            pi = np.clip(fi.astype(int), 0, H - 1)
            pj = np.clip(fj.astype(int), 0, W - 1)

            # --- 4.2) Compute the difference in altitude 
            #    (current position altitude - origin altitude).
            #    shape( H, W )
            #    But for each pixel (i,j), the "origin altitude" is start_alt[i,j].
            #    We'll do that in a vectorized manner.
            dh = src[pi, pj] - start_alt[i_coords, j_coords]

            # --- 4.3) Normalize dh by local_range_altitude for the origin pixel
            range_val = local_range_altitude[i_coords, j_coords]
            # A linear factor:  f_alt = 1 - (dh / range_val)
            # => if dh>0 => factor < 1 => reduces sigma
            # => if dh<0 => factor>1 => might increase sigma 
            #     (you might want to clamp that to 1 if you don't want to 
            #      strengthen low alt. portion).
            f_alt = 1.0 - (dh / range_val)
            # Optionally, clamp f_alt between [0, 1]
            f_alt = np.clip(f_alt, 0.0, 1.0)

            # --- 4.4) Adjust sigma accordingly
            sigma_adjusted = np.maximum(sigma_base * f_alt, 1e-6)

            # The distance from the center in pixel units
            distance = step_idx * DS

            # Gaussian weight
            weight = np.exp(-(distance**2) / (2 * sigma_adjusted**2))

            # --- 4.5) Accumulate results for active pixels
            # Add weighted altitude
            result[active_cells] += weight[active_cells] * src[pi[active_cells], pj[active_cells]]
            # Add to weights
            weights_sum[active_cells] += weight[active_cells]

            # --- 4.6) Decrement local steps for active pixels, 
            # then update the active mask
            current_steps[active_cells] -= 1
            active_cells = (current_steps > 0)

    # --- 5) Normalize final result to get the average
    np.divide(result, np.maximum(weights_sum, 1e-12), out=result)

    return result


def LIC_iterations(
    grid: np.ndarray,
    altitude_relative: np.ndarray,
    local_range_altitude: np.ndarray,
    cellsize: float,
    profile: Dict,
    f: np.ndarray,
    num_steps: int,
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
    altitude_relative : numpy.ndarray
        2D array of normalized relative altitudes (values between 0 and 1).
    local_range_altitude : numpy ndarray
        2D array of the range of elevation in a windows centered on each point
    cellsize : float
        The spatial resolution of the DEM (in meters).
    profile : dict
        Metadata associated with the raster.
    f : numpy.ndarray
        Continuous grid used to modulate the integration length.
    num_steps : int
        Maximum length of the integration line.
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
    1. Computes a modulated LIC and a standard LIC for the DEM.
    2. Combines the two LIC results using the relative altitude grid.
    3. Calculates the maximal curvature of the combined grid.
    4. Applies Gaussian smoothing to the curvature.
    5. Updates the DEM using a weighted combination of the original DEM and the combined LIC result.
    """
    for n in range(n_iterations):
        print(f"Iteration {n+1}/{n_iterations}")

        # # Step 1: Apply modulated LIC based on the continuous grid `f`
        # filtered_grid_modulated = extended_lic_weighted_altitude_lengthModulated(
        #     grid, altitude_relative, cellsize, f, num_steps
        # )

        # # Step 2: Apply standard LIC
        # filtered_grid_cretes = extended_lic_weighted_altitude(
        #     grid, local_range_altitude, cellsize, num_steps
        # )

        # # Step 3: Combine the two LIC results
        # filtered_grid_combined = (
        #     filtered_grid_cretes * altitude_relative
        #     + filtered_grid_modulated * np.abs(altitude_relative - 1)
        # )
        
        filtered_grid = extended_lic_weighted_altitude_lengthModulated(
            grid, local_range_altitude, cellsize, f, num_steps
        )

        # Step 4: Calculate maximal curvature of the combined grid
        filtered_grid_max_c = calculate_maximal_curvature(
            filtered_grid, cellsize
        )

        # Step 5: Smooth the curvature using Gaussian filtering
        filtered_grid_max_c_blured = gaussian_filter(
            filtered_grid_max_c, sigma=sigma, mode="nearest"
        )

        # Step 6: Compute weighting factor from the smoothed curvature
        weight = k * np.abs(filtered_grid_max_c_blured)
        weight = np.nan_to_num(weight, nan=0.0, posinf=1.0, neginf=0.0)

        # Step 7: Update the DEM using the weighted combination
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


def process_geotiff_with_overlap(
    MNT_input_path: str,
    output_path: str,
    block_size: int = 2000,
    overlap: int = 20,
    sigma_max: float = 5.0,
    slope_threshold: float = 0.1,
    num_bins: int = 10,
    min_area: int = 100,
    num_steps: int = 5,
    n_iterations: int = 5,
    sigma_blur_maxcurv: float = 3.0,
    k: float = 2.5,
) -> None:
    """Processes a GeoTIFF file in overlapping blocks and applies a series of
    transformations for terrain generalization.

    Parameters
    ----------
    MNT_input_path : str
        Path to the input GeoTIFF file containing the Digital Elevation Model (DEM).
    output_path : str
        Path to the output GeoTIFF file.
    block_size : int
        Size of non-overlapping blocks (in pixels) to process. Default is 2000.
    overlap : int
        Size of the overlapping region (in pixels) between adjacent blocks. Default is 20.
    sigma_max : float
        Maximum kernel width for the adaptive Gaussian blur. Default is 5.
    slope_threshold : float
        Threshold for slope to distinguish flat and steep areas. Default is 0.1.
    num_bins : int, optional
        Number of bins for sigma approximation in the adaptive blur. Default is 10.
    min_area : int, optional
        Minimum size (in pixels) for flat areas to be preserved. Default is 100.
    num_steps : int, optional
        Maximum integration length for the LIC algorithm. Default is 5.
    n_iterations : int, optional
        Number of iterations for the LIC algorithm. Default is 5.
    sigma_blur_maxcurv : float, optional
        Gaussian blur parameter for the maximum curvature. Default is 3.
    k : float, optional
        Weighting factor for combining LIC grids. Default is 2.5.

    Returns
    -------
    None
        The processed output is written directly to the specified GeoTIFF file.

    Notes
    -----
    - The function processes the input raster in overlapping blocks to handle large datasets efficiently.
    - Each block undergoes adaptive Gaussian smoothing, binary grid generation, LIC processing,
      and a series of corrections for smooth integration.
    - The overlapping region ensures seamless transitions between adjacent blocks.
    """

    # Open the input GeoTIFF and read metadata
    with rasterio.open(MNT_input_path) as MNT_src:
        profile = MNT_src.profile.copy()
        width, height = MNT_src.width, MNT_src.height
        cellsize = MNT_src.transform[0]

        # Update output profile for BigTIFF and float32 datatype
        profile.update(
            dtype=np.float32,
            compress="lzw",
            bigtiff="YES",
        )

        print("Starting terrain generalization with adaptive Gaussian blur")

        # Calculate the total number of blocks for progress tracking
        n_blocks_x = (width + block_size - 1) // block_size
        n_blocks_y = (height + block_size - 1) // block_size
        total_blocks = n_blocks_x * n_blocks_y

        # Create the output GeoTIFF file
        with rasterio.open(output_path, "w", **profile) as dst:
            # Use tqdm for tracking progress
            with tqdm(
                total=total_blocks, desc="Processing blocks", unit="block"
            ) as pbar:
                for y in range(0, height, block_size):
                    for x in range(0, width, block_size):
                        # Define window with overlap
                        x_min = max(0, x - overlap)
                        x_max = min(width, x + block_size + overlap)
                        y_min = max(0, y - overlap)
                        y_max = min(height, y + block_size + overlap)

                        window = Window(x_min, y_min, x_max - x_min, y_max - y_min)

                        # Read the MNT block
                        MNT_block = MNT_src.read(1, window=window)

                        # Compute the maximum curvature for the block
                        curvature_block = calculate_maximal_curvature(
                            MNT_block, cellsize
                        )

                        # Apply adaptive Gaussian blur
                        MNT_blurred_block = fast_adaptive_gaussian_blur(
                            MNT_block, curvature_block, cellsize, sigma_max, num_bins
                        )

                        # Calculate relative altitude
                        relative_altitude_block = calculate_relative_altitude(
                            MNT_blurred_block, window_size=2 * 20
                        )
                        
                        # Calculate local range of altitude
                        local_range_altitude_block = calculate_local_range(
                            MNT_blurred_block, steps = num_steps
                        )


                        # Generate binary grid for flat/steep areas
                        flat_steep_grid_block = initialize_flat_steep_grid(
                            MNT_blurred_block, slope_threshold
                        )

                        # Remove small flat areas
                        cleaned_grid_block = remove_small_flat_areas(
                            flat_steep_grid_block, min_area
                        )

                        # Smooth the binary grid for continuity
                        continuous_grid_block = gaussian_filter(
                            cleaned_grid_block.astype(float), sigma=5.0, mode="nearest"
                        )

                        # Correct flat area values
                        corrected_f = correct_flat_area_values(
                            cleaned_grid_block, continuous_grid_block, sigma=2.0
                        )

                        # Update the block profile for the current region
                        block_profile = profile.copy()
                        block_profile.update(
                            {
                                "height": MNT_block.shape[0],
                                "width": MNT_block.shape[1],
                                "transform": rasterio.windows.transform(
                                    window, MNT_src.transform
                                ),
                            }
                        )

                        # Perform multi-pass LIC processing
                        processed_block = LIC_iterations(
                            MNT_blurred_block,
                            relative_altitude_block,
                            local_range_altitude_block,
                            cellsize,
                            block_profile,
                            corrected_f,
                            num_steps,
                            n_iterations,
                            sigma_blur_maxcurv,
                            k,
                        )

                        # Extract the central block without overlap
                        if x == 0 and y != 0:
                            central_block = processed_block[
                                overlap : block_size + overlap, 0:block_size
                            ]
                        elif y == 0 and x != 0:
                            central_block = processed_block[
                                0:block_size, overlap : block_size + overlap
                            ]
                        elif y == 0 and x == 0:
                            central_block = processed_block[0:block_size, 0:block_size]
                        else:
                            central_block = processed_block[
                                overlap : block_size + overlap,
                                overlap : block_size + overlap,
                            ]

                        # Define output window for the central block
                        output_window = Window(
                            x, y, central_block.shape[1], central_block.shape[0]
                        )

                        # Write the processed block to the output file
                        dst.write(
                            central_block.astype(np.float32), 1, window=output_window
                        )

                        # Update progress bar
                        pbar.update(1)

    print(f"Processing completed. Output written to {output_path}")


if __name__ == "__main__":
    """Main script for processing a GeoTIFF file with Line Integral Convolution (LIC) and
    adaptive Gaussian blur techniques for terrain generalization.

y
    Notes
    -----
    This script initializes parameters for the process and calls the `process_geotiff_with_overlap`
    function to perform LIC-based generalization. It outputs a processed raster file
    saved in the specified folder.
    """


    # File paths
    MNT_input_path = "../../../../QGIS/out/Valtellina_NASADEM.tif"
    LIC_complete_output_path = "../../../../out_scripts_test_temp/Valtellina_NASADEM_LIC_article.tif"

    # Call the processing function
    process_geotiff_with_overlap(
        MNT_input_path,
        LIC_complete_output_path,
    )
