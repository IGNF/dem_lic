# DEM Generalization with Line Integral Convolution

This repository contains tools for Digital Elevation Model (DEM) generalization, leveraging Line Integral Convolution (LIC) and adaptive Gaussian blur techniques to simplify terrain representations while preserving key features. The method processes large DEM files efficiently by dividing them into overlapping blocks.

---

## Table of Contents
1. [Features](#features)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Documentation](#documentation)
5. [License](#license)
6. [Contributing](#contributing)
7. [Contact](#contact)

---

## Features
- Adaptive Gaussian blur for smoothing based on terrain curvature.
- Flat and steep area detection to improve feature generalization.
- Line Integral Convolution to highlight ridges and key terrain features.
- Iterative processing for enhanced visual coherence.

---

## Installation

Clone this repository to your local drive.
You can then install this project directly using pip:

```bash
pip install .
```

Ensure you have the following dependencies installed:
- `numpy`
- `rasterio`
- `scipy`
- `tqdm`

---

## Usage

### Command-line Usage
Process a GeoTIFF DEM file by calling the main processing function directly:
```python
from src.LIC_extended_main import process_geotiff_with_overlap

process_geotiff_with_overlap(
    MNT_input_path="path_to_input_file.tif",
    output_path="path_to_output_file.tif"
)
```

### Using in Python Scripts
You can import and use individual functions directly in Python, or run the main processing function for more complex workflows:

```python
from src.LIC_extended_main import process_geotiff_with_overlap
from src.morphoDEM import calculate_maximal_curvature

# Example: Compute curvature
curvature = calculate_maximal_curvature(dem_array, resolution=1.0)

# Example: Process a DEM file
process_geotiff_with_overlap(
    MNT_input_path="path_to_input_file.tif",
    output_path="path_to_output_file.tif",
    block_size=2000,
    overlap=20,
    sigma_max=5,
    slope_threshold=0.1,
    num_bins=10,
    min_area=100,
    num_steps=5,
    n_iterations=5,
    sigma_blur_maxcurv=3,
    k=2.5
)
```

---

## Documentation
The full methodology is detailed in the document:
[Methodological Documentation DEM Generalization with Line Integral Convolution](docs/Methodological%20Documentation%20DEM%20Generalization%20with%20Line%20Integral%20Convolution.pdf)

Additionally, the HTML version of the documentation is available and can be viewed locally by opening:
`docs/build/html/index.html`

---

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

## Contributing

Contributions are welcome! If you want to report a bug, suggest improvements, or contribute code, feel free to open an issue or submit a pull request.

---

## Contacts
For further information, questions, or feedback, please contact:
**Edmond SAINT DENIS** - [edmond.saint-denis@ign.fr](mailto:edmond.saint-denis@ign.fr)