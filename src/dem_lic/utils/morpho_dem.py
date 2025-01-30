# -*- coding: utf-8 -*-
"""
DEM Processing Utilities
"""

import numpy as np
from scipy.ndimage import sobel
from scipy.ndimage import gaussian_filter, label, maximum_filter, minimum_filter


def calculate_maximal_curvature(
        dem: np.ndarray,
        resolution: float = 1
) -> np.ndarray:
    """
    Compute the maximal curvature (k_max) of a DEM using vectorized calculations.

    Parameters
    ----------
    dem : numpy.ndarray
        A 2D array representing the Digital Elevation Model (DEM).
    resolution : float, optional
        The spatial resolution of the DEM (default is 1.0).

    Returns
    -------
    numpy.ndarray
        A 2D array containing the maximal curvature values.

    Notes
    -----
    This function computes the following terms:
    - p: dz/dx (first derivative in x direction)
    - q: dz/dy (first derivative in y direction)
    - r: d²z/dx² (second derivative in x direction)
    - t: d²z/dy² (second derivative in y direction)
    - s: d²z/dxdy (second mixed derivative)

    The maximal curvature is then computed as:
    k_max = H + sqrt(H² - K)
    where:
    H = -((1 + q²)r - 2pqs + (1 + p²)t) / (2 * (1 + p² + q²)^(3/2))
    K = (rt - s²) / (1 + p² + q²)²
    """
    # Compute first derivatives (p, q)
    dz_dx = sobel(dem, axis=1, mode="nearest") / (8 * resolution)  # dz/dx
    dz_dy = sobel(dem, axis=0, mode="nearest") / (8 * resolution)  # dz/dy

    # Compute second derivatives (r, t, s)
    d2z_dx2 = sobel(dz_dx, axis=1, mode="nearest") / (8 * resolution)  # d²z/dx²
    d2z_dy2 = sobel(dz_dy, axis=0, mode="nearest") / (8 * resolution)  # d²z/dy²
    d2z_dxdy = sobel(dz_dx, axis=0, mode="nearest") / (8 * resolution)  # d²z/dxdy

    # Intermediate terms
    p2 = dz_dx**2
    q2 = dz_dy**2
    pq2_sum = 1 + p2 + q2

    # Mean curvature (H)
    H = -((1 + q2) * d2z_dx2 - 2 * dz_dx * dz_dy * d2z_dxdy + (1 + p2) * d2z_dy2) / (
        2 * pq2_sum ** (3 / 2)
    )

    # Gaussian curvature (K)
    K = (d2z_dx2 * d2z_dy2 - d2z_dxdy**2) / pq2_sum**2

    # Maximal curvature (k_max)
    k_max = H + np.sqrt(
        np.maximum(H**2 - K, 0)
    )  # Ensure non-negative values inside sqrt

    return k_max


def fast_adaptive_gaussian_blur(
    grid: np.ndarray,
    curvature: np.ndarray,
    cellsize: float,
    sigma_max: float,
    num_bins: int = 10,
) -> np.ndarray:
    """Applies a fast adaptive Gaussian blur by grouping pixels into bins based on local curvature values.

    Parameters
    ----------
    grid : numpy.ndarray
        A 2D array representing the Digital Elevation Model (DEM).
    curvature : numpy.ndarray
        A 2D array of local curvature values corresponding to the DEM.
    cellsize : float
        The size of the raster cells in meters.
    sigma_max : float
        The maximum kernel width for the Gaussian blur.
    num_bins : int, optional
        The number of bins used to approximate the sigma values. Default is 10.

    Returns
    -------
    numpy.ndarray
        A 2D array representing the smoothed DEM after applying the adaptive Gaussian blur.

    Notes
    -----
    The function performs the following steps:
    1. Calculates the sigma value for each pixel based on the local curvature, limited by `sigma_max`.
    2. Groups pixels into bins according to their sigma values.
    3. Applies a Gaussian blur for each bin using a fixed sigma value.
    4. Interpolates results to produce a smooth output.

    The adaptive Gaussian blur ensures that flat areas are smoothed more aggressively while preserving
    sharp features like ridges and peaks.
    """

    # Step 1: Compute sigma values for each pixel based on curvature
    sigma = 1 / (np.abs(curvature) * cellsize + 1e-6)  # Avoid division by zero
    sigma = np.clip(
        sigma, 0, sigma_max
    )  # Restrict sigma values to a maximum of sigma_max

    # Step 2: Group pixels into bins based on their sigma values
    bins = np.linspace(0, sigma_max, num_bins + 1)  # Define bin edges
    bin_indices = np.digitize(sigma, bins) - 1  # Assign each pixel to a bin
    bin_indices = np.clip(
        bin_indices, 0, num_bins - 1
    )  # Ensure indices stay within valid range

    # Step 3: Apply Gaussian blur for each bin
    blurred_grids = []
    for b in range(num_bins):
        sigma_bin = bins[b]  # Sigma value for the current bin
        if sigma_bin > 0:  # Only apply Gaussian blur if sigma is greater than 0
            blurred = gaussian_filter(
                grid, sigma=sigma_bin, mode="nearest"
            )  # Apply Gaussian blur
        else:
            blurred = grid.copy()  # No smoothing for sigma = 0
        blurred_grids.append(blurred)  # Store the result for this bin

    # Step 4: Interpolate results between bins (optimized vectorized approach)

    # Create arrays for the lower and upper bin indices
    lower_idx = np.clip(bin_indices, 0, num_bins - 2)  # Lower bin index
    upper_idx = lower_idx + 1  # Upper bin index

    # Compute the interpolation factors
    t = (sigma - bins[lower_idx]) / (
        bins[upper_idx] - bins[lower_idx] + 1e-6
    )  # Avoid division by zero

    # Prepare the smoothed grid by interpolating between bins
    smoothed_grid = (1 - t) * np.choose(lower_idx, blurred_grids) + t * np.choose(
        upper_idx, blurred_grids
    )
    return smoothed_grid


# Chapter 3.6 : Preserving sharp transitions to flat areas
def initialize_flat_steep_grid(
        mnt: np.ndarray,
        slope_threshold: float
) -> np.ndarray:
    """Creates a binary grid to distinguish flat and steep areas based on slope
    calculated from the DEM.

    Parameters
    ----------
    mnt : numpy.ndarray
        2D array representing the Digital Elevation Model (DEM).
    slope_threshold : float
        Threshold value for the slope to separate flat areas (0) from steep areas (1).

    Returns
    -------
    numpy.ndarray
        2D binary array with 0 for flat areas and 1 for steep areas.

    Notes
    -----
    The slope is computed using the first-order derivatives in the x and y directions,
    combined to form the magnitude of the gradient. Areas with a slope below the
    threshold are considered flat, while others are marked as steep.

    """
    # Step 1: Compute the partial derivatives in x and y directions
    dz_dx = np.gradient(mnt, axis=1)  # Gradient in the x-direction
    dz_dy = np.gradient(mnt, axis=0)  # Gradient in the y-direction

    # Step 2: Compute the slope magnitude
    slope = np.sqrt(dz_dx**2 + dz_dy**2)

    # Step 3: Create a binary grid based on the slope threshold
    flat_steep_grid = np.where(slope < slope_threshold, 0, 1)

    return flat_steep_grid


def remove_small_flat_areas(
        flat_steep_grid: np.ndarray,
        min_area: int
) -> np.ndarray:
    """Removes small flat areas from a binary grid by replacing regions smaller
    than the specified threshold with steep areas.

    Parameters
    ----------
    flat_steep_grid : numpy.ndarray
        2D binary array representing the terrain, with 0 for flat areas and 1 for steep areas.
    min_area : int
        Minimum size of flat regions (in pixels) to retain.

    Returns
    -------
    numpy.ndarray
        2D binary array where small flat areas have been replaced with 1 (steep areas).

    Notes
    -----
    Small flat areas are identified as connected regions of 0-values that are smaller
    than `min_area`. These regions are replaced with steep areas (1) in the output grid.

    """
    # Step 1: Label connected flat regions (regions of 0)
    labeled_grid, num_features = label(flat_steep_grid == 0)

    # Step 2: Compute the size of each labeled region
    region_sizes = np.bincount(labeled_grid.ravel())

    # Step 3: Identify small regions with size less than `min_area`
    small_regions = region_sizes < min_area
    small_regions[0] = False  # Ignore background pixels (label 0)

    # Step 4: Replace small flat regions with steep areas
    cleaned_grid = flat_steep_grid.copy()
    cleaned_grid[np.isin(labeled_grid, np.where(small_regions)[0])] = 1

    return cleaned_grid


def calculate_relative_altitude(
        mnt: np.ndarray,
        window_size: int = 40
) -> np.ndarray:
    """Computes a raster of normalized relative altitude values between 0 and 1
    based on local elevation. This is intended to weight the combination
    of ridge enhancement and flat area transition treatments.

    Parameters
    ----------
    mnt : numpy.ndarray
        2D array representing the Digital Elevation Model (DEM).
    window_size : int, optional
        Size of the window (in pixels) used to compute local maximums and minimums.
        Default is 40.

    Returns
    -------
    numpy.ndarray
        2D array of normalized relative altitude values between 0 and 1.

    Notes
    -----
    The function calculates relative altitude as the difference between the current
    pixel value and the local minimum within a moving window, normalized by the range
    (local maximum - local minimum). Flat areas with no variation are assigned a neutral value of 0.5.
    """
    # Step 1: Calculate local maximum values within the specified window size
    max_local = maximum_filter(mnt, size=window_size, mode="reflect")

    # Step 2: Calculate local minimum values within the specified window size
    min_local = minimum_filter(mnt, size=window_size, mode="reflect")

    # Step 3: Compute the relative altitude as a normalized value
    denominator = max_local - min_local  # Compute the range in the window
    with np.errstate(divide="ignore", invalid="ignore"):  # Handle divisions by zero
        relative_altitude = (mnt - min_local) / denominator  # Normalize elevation
        relative_altitude[denominator == 0] = 0.5  # Assign neutral value to flat areas

    # Step 4: Ensure all values are clipped between 0 and 1
    relative_altitude = np.clip(relative_altitude, 0, 1)

    return relative_altitude


def calculate_local_range(
        mnt: np.ndarray,
        steps: int = 5
) -> np.ndarray:
    """Computes the altitude range in a window. This is intended to weight the elevation
    of each point in the line integral convolution.

    Parameters
    ----------
    mnt : numpy.ndarray
        2D array representing the Digital Elevation Model (DEM).
    steps : int, optional
        Mid size of the window (in pixels) used to compute local maximums and minimums.
        Default is 5.

    Returns
    -------
    numpy.ndarray
        2D array of local elevation range.

    Notes
    -----
    The function calculates the difference between the local
    maximum value and the local minimum within a moving window
    (local maximum - local minimum).
    """
    
    window_size = 2* steps
    max_local = maximum_filter(mnt, size=window_size, mode="reflect")
    min_local = minimum_filter(mnt, size=window_size, mode="reflect")
    
    range_alt = np.abs(max_local-min_local)


    return range_alt

