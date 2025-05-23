# -*- coding: utf-8 -*-
"""
Main code to call different stages of the processing described in the article :
"Bernhard Jenny (2021) Terrain generalization with line integral
convolution, Cartography and Geographic Information Science, 48:1, 78-92, DOI:
10.1080/15230406.2020.1833762
https://doi.org/10.1080/15230406.2020.1833762"
"""

import os


from dem_lic.utils.lic_extended import process_lic_extended, process_geotiff_in_block_with_overlap
from functools import partial



def generalization(
    MNT_input_path: str,
    output_path: str,
    block_size: int = 2000,
    overlap: int = 20,
    sigma_max: float = 5.0,
    slope_threshold: float = 6,
    num_bins: int = 10,
    min_area: int = 100,
    num_steps: int = 5,
    sigma_modulated: bool = True,
    n_iterations: int = 5,
    sigma_blur_maxcurv: float = 3.0,
    k: float = 2.5,
) -> None:
    """
    Main function to validate inputs and call the `process_geotiff_with_overlap` function.
    Processes a GeoTIFF file in overlapping blocks and applies a series of
    transformations for terrain generalization.

    Parameters
    ----------
    MNT_input_path : str
        Path to the input GeoTIFF file containing the Digital Elevation Model (DEM).
    output_path : str
        Path to the output GeoTIFF file.
    block_size : int, optional
        Size of non-overlapping blocks (in pixels) to process. Default is 2000.
    overlap : int, optional
        Size of the overlapping region (in pixels) between adjacent blocks. Default is 20.
    sigma_max : float, optional
        Maximum kernel width for the adaptive Gaussian blur. Default is 5.0.
    slope_threshold : float, optional
        Threshold for slope to distinguish flat and steep areas, in degrees. Default is 6.
    num_bins : int, optional
        Number of bins for sigma approximation in the adaptive blur. Default is 10.
    min_area : int, optional
        Minimum size (in pixels) for flat areas to be preserved. Default is 100.
    num_steps : int, optional
        Maximum integration length for the LIC algorithm. Default is 5.
     sigma_modulated : bool, optional
        If True, the sigma for the weighting in the LIC algorithm is modulated by the altitude.
        If False, a fixed Gaussian sigma is used.
    n_iterations : int, optional
        Number of iterations for the LIC algorithm. Default is 5.
    sigma_blur_maxcurv : float, optional
        Gaussian blur parameter for the maximum curvature. Default is 3.0.
    k : float, optional
        Weighting factor for combining LIC grids. Default is 2.5.

    Returns
    -------
    None
        The processed output is written directly to the specified GeoTIFF file.

    Notes
    -----
    This function ensures that the inputs are valid before executing the main processing pipeline.
    """

    # Check if input file exists
    if not os.path.isfile(MNT_input_path):
        raise FileNotFoundError(f"Input file not found: {MNT_input_path}")

    # Check that block size and overlap are positive integers
    if not isinstance(block_size, int) or block_size <= 0:
        raise ValueError("Block size must be a positive integer.")
    if not isinstance(overlap, int) or overlap < 0:
        raise ValueError("Overlap must be a non-negative integer.")

    # Check that sigma_max, slope_threshold, sigma_blur_maxcurv, and k are positive floats
    if not isinstance(sigma_max, (int, float)) or sigma_max <= 0:
        raise ValueError("Sigma max must be a positive number.")
    if not isinstance(slope_threshold, (int, float)) or slope_threshold < 0:
        raise ValueError("Slope threshold must be a positive number.")
    if not isinstance(sigma_blur_maxcurv, (int, float)) or sigma_blur_maxcurv <= 0:
        raise ValueError("Sigma blur for maximum curvature must be a positive number.")
    if not isinstance(k, (int, float)) or k <= 0:
        raise ValueError("Weighting factor k must be a positive number.")

    # Check that num_bins, min_area, num_steps, and n_iterations are positive integers
    if not isinstance(num_bins, int) or num_bins <= 0:
        raise ValueError("Number of bins must be a positive integer.")
    if not isinstance(min_area, int) or min_area <= 0:
        raise ValueError("Minimum area must be a positive integer.")
    if not isinstance(num_steps, int) or num_steps <= 0:
        raise ValueError("Number of steps must be a positive integer.")
    if not isinstance(n_iterations, int) or n_iterations <= 0:
        raise ValueError("Number of iterations must be a positive integer.")

    # Check that sigma_modulated is a boolean
    if not isinstance(sigma_modulated, bool):
        raise ValueError("sigma_modulated must be a boolean (True or False).")


    
    # Create a partially configured version of process_lic_extended
    lic_extended_partial = partial(
        process_lic_extended,
        sigma_max=sigma_max,
        slope_threshold=slope_threshold,
        num_bins=num_bins,
        min_area=min_area,
        num_steps=num_steps,
        sigma_modulated=sigma_modulated,
        n_iterations=n_iterations,
        sigma_blur_maxcurv=sigma_blur_maxcurv,
        k=k,
    )
    
    # Call the processing function with the pre-configured LIC function
    process_geotiff_in_block_with_overlap(
        MNT_input_path=MNT_input_path,
        output_path=output_path,
        processing_function=lic_extended_partial,
        block_size=block_size,
        overlap=overlap
    )

    print("Processing completed successfully.")


if __name__ == "__main__":
    MNT_input_path = "C:\Projets\cretes_talwegs\QGIS\out\elevation_1KMmd_GMTEDmd_metrique_extrait.tif"
    LIC_complete_output_path = "C:\Projets\cretes_talwegs\out_scripts_test_temp\elevation_1KMmd_GMTEDmd_metrique_extrait_LIC_st0.tif"

    
    # Call the processing function
    generalization(
        MNT_input_path,
        LIC_complete_output_path,
        sigma_max = 2.0,
        slope_threshold = 0,
        num_steps = 5,
        sigma_modulated= True,
        n_iterations = 4,
        sigma_blur_maxcurv = 2.9,

    )

    
