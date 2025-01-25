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

````{py:function} extended_lic_weighted_altitude_lengthModulated(grid, relative_altitude, cellsize, f, num_steps)
:canonical: dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated

```{autodoc2-docstring} dem_lic.utils.lic_extended.extended_lic_weighted_altitude_lengthModulated
```
````

````{py:function} LIC_iterations(grid, altitude_relative, cellsize, profile, f, num_steps, n_iterations, sigma, k)
:canonical: dem_lic.utils.lic_extended.LIC_iterations

```{autodoc2-docstring} dem_lic.utils.lic_extended.LIC_iterations
```
````

````{py:function} correct_flat_area_values(b, c, sigma, epsilon=1e-06)
:canonical: dem_lic.utils.lic_extended.correct_flat_area_values

```{autodoc2-docstring} dem_lic.utils.lic_extended.correct_flat_area_values
```
````

````{py:function} process_geotiff_with_overlap(MNT_input_path, output_path, block_size=2000, overlap=20, sigma_max=5, slope_threshold=0.1, num_bins=10, min_area=100, num_steps=5, n_iterations=5, sigma_blur_maxcurv=3, k=2.5)
:canonical: dem_lic.utils.lic_extended.process_geotiff_with_overlap

```{autodoc2-docstring} dem_lic.utils.lic_extended.process_geotiff_with_overlap
```
````
