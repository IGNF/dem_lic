# DEM Generalization with LIC

[![License: MIT](./docs/_static/License-MIT.svg)](./LICENSE.txt)
[![Python](./docs/_static/Python-3.8.svg)](https://www.python.org/)


A Python package for Digital Elevation Model (DEM) generalization, leveraging Line Integral Convolution (LIC) and adaptive Gaussian blur techniques to simplify terrain representations while preserving key features. The method processes large DEM files efficiently by dividing them into overlapping blocks.

This implementation is based on the work presented in the article:

```{note}
Bernhard Jenny (2021) Terrain generalization with line integral convolution, Cartography and Geographic Information Science, 48:1, 78-92, DOI: 10.1080/15230406.2020.1833762
```
[https://doi.org/10.1080/15230406.2020.1833762](https://doi.org/10.1080/15230406.2020.1833762)


## Documentation

Full documentation is available at:

[https://ignf.github.io/dem_lic/](https://ignf.github.io/dem_lic/)

## What does it do?

This algorithm generalizes Digital Elevation Models (DEMs) by combining adaptive Gaussian blur based on terrain curvature, flat and steep area detection for feature enhancement, and Line Integral Convolution (LIC) to emphasize ridges and key terrain features. Iterative processing ensures improved visual coherence and clarity.

![Exemple of generalization](docs/images/dem_to_generalization_25m.png)

## Installation

To install the package, clone the repository and install it in editable mode.

### HTTP method

```bash
git clone https://github.com/IGNF/dem_lic.git
cd dem_lic
pip install -e .
```

### SSH method

```bash
git clone git@github.com:IGNF/dem_lic.git
cd dem_lic
pip install -e .
```

## Command Line Quick Start

The command-line interface (CLI) for **dem_lic** provides an easy way to process Digital Elevation Models (DEMs) using the generalization algorithm. The CLI supports input and output file specifications, as well as customization of processing parameters.

```bash
dem_lic input_dem.tif output_generalized_dem.tif
```

Optional arguments can be used to adjust the process:

```bash
dem_lic input_dem.tif output_generalized_dem.tif --n_iterations 5 --overlap 20
```

For detailed usage instructions and a complete list of options, refer to the [CLI Documentation](./docs/documentation/cli.md).

## Python Quick Start

You can also use **dem_lic** directly within your Python scripts. Here is a quick example:

```python
from dem_lic.main_lic_extended import generalization

generalization(
    MNT_input_path="path_to_input_dem.tif",
    output_path="path_to_output_dem.tif",
    sigma_max = 2.0,
    n_iterations = 4,
)
```

## Input DEM requirements

The input Digital Elevation Model must:

- Be a **GeoTIFF** file.
- Use a **projected coordinate system with metric units** (e.g., UTM, Lambert 93).  
Do **not** use a geographic CRS with degrees (such as WGS84 / EPSG:4326), as the algorithm relies on metric distances for slope, curvature, and LIC integration lengths.  
Using degrees will result in **incorrect and unusable outputs**.


## License

This project is licensed under the MIT License. See the LICENSE.TXT file for details.
