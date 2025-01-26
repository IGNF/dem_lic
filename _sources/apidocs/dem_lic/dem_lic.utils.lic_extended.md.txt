# {py:mod}`dem_lic.utils.lic_extended`

```{py:module} dem_lic.utils.lic_extended
```

```{autodoc2-docstring} dem_lic.utils.lic_extended
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`extended_lic_weighted_altitude <dem_lic.utils.lic_extended.extended_lic_weighted_altitude>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude
    :summary:
    ```
* - {py:obj}`extended_lic_weighted_altitude_lengthModulated <dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated
    :summary:
    ```
* - {py:obj}`LIC_iterations <dem_lic.utils.lic_extended.LIC_iterations>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.LIC_iterations
    :summary:
    ```
* - {py:obj}`correct_flat_area_values <dem_lic.utils.lic_extended.correct_flat_area_values>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.correct_flat_area_values
    :summary:
    ```
* - {py:obj}`process_geotiff_with_overlap <dem_lic.utils.lic_extended.process_geotiff_with_overlap>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.process_geotiff_with_overlap
    :summary:
    ```
````

### API

````{py:function} extended_lic_weighted_altitude(grid: numpy.ndarray, relative_altitude: numpy.ndarray, cellsize: float, num_steps: int) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.extended_lic_weighted_altitude

```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude
```
````

````{py:function} extended_lic_weighted_altitude_lengthModulated(grid: numpy.ndarray, relative_altitude: numpy.ndarray, cellsize: float, f: numpy.ndarray, num_steps: int) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated

```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated
```
````

````{py:function} LIC_iterations(grid: numpy.ndarray, altitude_relative: numpy.ndarray, cellsize: float, profile: typing.Dict, f: numpy.ndarray, num_steps: int, n_iterations: int, sigma: float, k: float) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.LIC_iterations

```{autodoc2-docstring} dem_lic.utils.lic_extended.LIC_iterations
```
````

````{py:function} correct_flat_area_values(b: numpy.ndarray, c: numpy.ndarray, sigma: float, epsilon: float = 1e-06) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.correct_flat_area_values

```{autodoc2-docstring} dem_lic.utils.lic_extended.correct_flat_area_values
```
````

````{py:function} process_geotiff_with_overlap(MNT_input_path: str, output_path: str, block_size: int = 2000, overlap: int = 20, sigma_max: float = 5.0, slope_threshold: float = 0.1, num_bins: int = 10, min_area: int = 100, num_steps: int = 5, n_iterations: int = 5, sigma_blur_maxcurv: float = 3.0, k: float = 2.5) -> None
:canonical: dem_lic.utils.lic_extended.process_geotiff_with_overlap

```{autodoc2-docstring} dem_lic.utils.lic_extended.process_geotiff_with_overlap
```
````
