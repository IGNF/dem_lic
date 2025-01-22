Command Line Interface (CLI)
============================

The command-line interface (CLI) for **dem_lic** allows users to process Digital Elevation Models (DEMs) directly from the terminal using advanced generalization techniques. This document provides detailed information on how to use the CLI, including available commands, options, and examples.

Usage
-----

To run the CLI, execute the following command:

.. code-block:: bash

    dem_lic -h

Use the `-h` option to display help and see all available arguments and their descriptions.

.. code-block:: bash

    dem_lic <MNT_input_path> <output_path> [options]


**Positional Arguments**:

- **MNT_input_path**: Path to the input GeoTIFF file containing the DEM.
- **output_path**: Path to the output GeoTIFF file where the generalized DEM will be saved.

**Optional Arguments**:

The following optional arguments can be used to customize the processing:

- `--block_size` (int): Size of processing blocks in pixels (default: 2000).
- `--overlap` (int): Size of the overlap in pixels between blocks (default: 20).
- `--sigma_max` (float): Maximum kernel width for the adaptive Gaussian blur (default: 5.0).
- `--slope_threshold` (float): Slope threshold for distinguishing flat and steep areas (default: 0.1).
- `--num_bins` (int): Number of bins for sigma approximation (default: 10).
- `--min_area` (int): Minimum size of flat areas to preserve, in pixels (default: 100).
- `--num_steps` (int): Maximum integration length for the LIC algorithm (default: 5).
- `--n_iterations` (int): Number of LIC iterations (default: 5).
- `--sigma_blur_maxcurv` (float): Gaussian blur sigma for maximal curvature (default: 3.0).
- `--k` (float): Weighting factor for combining LIC results (default: 2.5).

Examples
--------

**Basic Usage**:

Process a DEM with default settings:

.. code-block:: bash

    dem_lic input_dem.tif generalized_dem.tif

**Custom Parameters**:

Specify a custom block size and overlap:

.. code-block:: bash

    dem_lic input_dem.tif generalized_dem.tif --block_size 3000 --overlap 50

Increase the number of iterations for better generalization:

.. code-block:: bash

    dem_lic input_dem.tif generalized_dem.tif --n_iterations 10

**Error Handling**:

If the input file is missing or invalid, the CLI will raise an error:

.. code-block:: bash

    FileNotFoundError: Input file not found: input_dem.tif

If an invalid value is provided for a parameter, the CLI will report an error:

.. code-block:: bash

    ValueError: Block size must be a positive integer.

Development Notes
-----------------

The CLI is implemented in the `cli.py` module and serves as an entry point for the **dem_lic** package. It validates user inputs, parses command-line arguments, and calls the `process_geotiff_with_overlap` function from the `utils.lic_extended` module.

For additional details on the processing algorithm, refer to the code documentation or the source code of `cli.py`.


