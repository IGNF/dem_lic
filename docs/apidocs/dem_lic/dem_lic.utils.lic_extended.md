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
* - {py:obj}`process_geotiff_in_block_with_overlap <dem_lic.utils.lic_extended.process_geotiff_in_block_with_overlap>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.process_geotiff_in_block_with_overlap
    :summary:
    ```
* - {py:obj}`process_lic_extended <dem_lic.utils.lic_extended.process_lic_extended>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.process_lic_extended
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`lic_extended_partial <dem_lic.utils.lic_extended.lic_extended_partial>`
  - ```{autodoc2-docstring} dem_lic.utils.lic_extended.lic_extended_partial
    :summary:
    ```
````

### API

````{py:function} extended_lic_weighted_altitude_lengthModulated(grid: numpy.ndarray, local_range_altitude: numpy.ndarray, cellsize: float, f: numpy.ndarray, num_steps: int, sigma_modulated: bool = True) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated

```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated
```
````

````{py:function} LIC_iterations(grid: numpy.ndarray, local_range_altitude: numpy.ndarray, cellsize: float, f: numpy.ndarray, num_steps: int, sigma_modulated: bool, n_iterations: int, sigma: float, k: float) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.LIC_iterations

```{autodoc2-docstring} dem_lic.utils.lic_extended.LIC_iterations
```
````

````{py:function} correct_flat_area_values(b: numpy.ndarray, c: numpy.ndarray, sigma: float, epsilon: float = 1e-06) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.correct_flat_area_values

```{autodoc2-docstring} dem_lic.utils.lic_extended.correct_flat_area_values
```
````

````{py:function} process_geotiff_in_block_with_overlap(MNT_input_path: str, output_path: str, processing_function, block_size: int = 2000, overlap: int = 20, **kwargs) -> None
:canonical: dem_lic.utils.lic_extended.process_geotiff_in_block_with_overlap

```{autodoc2-docstring} dem_lic.utils.lic_extended.process_geotiff_in_block_with_overlap
```
````

````{py:function} process_lic_extended(MNT_block: numpy.ndarray, cellsize: float, sigma_max: float = 5.0, slope_threshold: float = 6, num_bins: int = 10, min_area: int = 100, num_steps: int = 5, sigma_modulated: bool = True, n_iterations: int = 5, sigma_blur_maxcurv: float = 3.0, k: float = 2.5) -> numpy.ndarray
:canonical: dem_lic.utils.lic_extended.process_lic_extended

```{autodoc2-docstring} dem_lic.utils.lic_extended.process_lic_extended
```
````

````{py:data} lic_extended_partial
:canonical: dem_lic.utils.lic_extended.lic_extended_partial
:value: >
   'partial(...)'

```{autodoc2-docstring} dem_lic.utils.lic_extended.lic_extended_partial
```

````
