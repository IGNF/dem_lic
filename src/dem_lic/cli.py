# -*- coding: utf-8 -*-
"""
CLI interface to call different stages of the processing described in the article :
"Bernhard Jenny (2021) Terrain generalization with line integral
convolution, Cartography and Geographic Information Science, 48:1, 78-92, DOI:
10.1080/15230406.2020.1833762
https://doi.org/10.1080/15230406.2020.1833762"
"""



import argparse
import os


from dem_lic.utils.lic_extended import process_geotiff_with_overlap


def main():
    """Main function to parse command-line arguments, validate inputs, and call
    the `process_geotiff_with_overlap` function.
    """
    description = (
        "Processes a GeoTIFF file using Line Integral Convolution (LIC) "
        "and adaptive Gaussian blur techniques for terrain generalization.\n\n"
        "Examples:\n"
        "  dem-lic input.tif output.tif --block_size 2000 --overlap 50\n"
        "  dem-lic input.tif output.tif --n_iterations 10 --sigma_max 6.0\n"
    )
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)


    # Add arguments
    parser.add_argument("MNT_input_path", type=str, help="Path to the input GeoTIFF file (DEM).")
    parser.add_argument("output_path", type=str, help="Path to the output generalized GeoTIFF file.")
    parser.add_argument("--block_size", type=int, default=2000, help="Size of processing blocks in pixels. Default: 2000.")
    parser.add_argument("--overlap", type=int, default=20, help="Size of the overlap in pixels. Default: 20.")
    parser.add_argument("--sigma_max", type=float, default=5.0, help="Maximum Gaussian kernel width. Default: 5.0.")
    parser.add_argument("--slope_threshold", type=float, default=6, help="Slope threshold for flat/steep distinction in degrees. Default: 6.")
    parser.add_argument("--num_bins", type=int, default=10, help="Number of bins for sigma approximation. Default: 10.")
    parser.add_argument("--min_area", type=int, default=100, help="Minimum size of flat areas to preserve (in pixels). Default: 100.")
    parser.add_argument("--num_steps", type=int, default=5, help="Maximum integration length for LIC. Default: 5.")
    parser.add_argument("--sigma_modulated", type=str, default="true", choices=["true", "false"],
                        help="Whether altitude weighting is used in the LIC algorithm (true/false). Default: true.")    
    parser.add_argument("--n_iterations", type=int, default=5, help="Number of LIC iterations. Default: 5.")
    parser.add_argument("--sigma_blur_maxcurv", type=float, default=3.0, help="Gaussian blur sigma for max curvature. Default: 3.0.")
    parser.add_argument("--k", type=float, default=2.5, help="Weighting factor for combining LIC results. Default: 2.5.")
    parser.add_argument("--version", action="version", version="dem_lic 1.0.0")

    # Parse arguments
    args = parser.parse_args()

    # Validate input file
    if not os.path.isfile(args.MNT_input_path):
        raise FileNotFoundError(f"Input file not found: {args.MNT_input_path}")

    # Validate block_size and overlap
    if args.block_size <= 0:
        raise ValueError("Block size must be a positive integer.")
    if args.overlap < 0:
        raise ValueError("Overlap must be a non-negative integer.")

    # Validate numeric parameters
    if args.sigma_max <= 0:
        raise ValueError("Sigma max must be a positive number.")
    if args.slope_threshold <= 0:
        raise ValueError("Slope threshold must be a positive number.")
    if args.sigma_blur_maxcurv <= 0:
        raise ValueError("Sigma blur for max curvature must be a positive number.")
    if args.k <= 0:
        raise ValueError("Weighting factor k must be a positive number.")
    if args.num_bins <= 0:
        raise ValueError("Number of bins must be a positive integer.")
    if args.min_area <= 0:
        raise ValueError("Minimum area must be a positive integer.")
    if args.num_steps <= 0:
        raise ValueError("Number of steps must be a positive integer.")
    if args.n_iterations <= 0:
        raise ValueError("Number of iterations must be a positive integer.")
        
    
    # Convert `sigma_modulated` from string to boolean
    sigma_modulated = args.sigma_modulated.lower() == "true"

    # Call the processing function
    process_geotiff_with_overlap(
        MNT_input_path=args.MNT_input_path,
        output_path=args.output_path,
        block_size=args.block_size,
        overlap=args.overlap,
        sigma_max=args.sigma_max,
        slope_threshold=args.slope_threshold,
        num_bins=args.num_bins,
        min_area=args.min_area,
        num_steps=args.num_steps,
        sigma_modulated=sigma_modulated,
        n_iterations=args.n_iterations,
        sigma_blur_maxcurv=args.sigma_blur_maxcurv,
        k=args.k,
    )

    print(f"Processing completed successfully. Output saved to: {args.output_path}")


if __name__ == "__main__":
    main()
