# {py:mod}`dem_lic.utils.morpho_dem`

```{py:module} dem_lic.utils.morpho_dem
```

```{autodoc2-docstring} dem_lic.utils.morpho_dem
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`calculate_maximal_curvature <dem_lic.utils.morpho_dem.calculate_maximal_curvature>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_maximal_curvature
    :summary:
    ```
* - {py:obj}`fast_adaptive_gaussian_blur <dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur
    :summary:
    ```
* - {py:obj}`initialize_flat_steep_grid <dem_lic.utils.morpho_dem.initialize_flat_steep_grid>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.initialize_flat_steep_grid
    :summary:
    ```
* - {py:obj}`remove_small_flat_areas <dem_lic.utils.morpho_dem.remove_small_flat_areas>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.remove_small_flat_areas
    :summary:
    ```
* - {py:obj}`calculate_local_range <dem_lic.utils.morpho_dem.calculate_local_range>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_local_range
    :summary:
    ```
````

### API

````{py:function} calculate_maximal_curvature(dem: numpy.ndarray, resolution: float = 1) -> numpy.ndarray
:canonical: dem_lic.utils.morpho_dem.calculate_maximal_curvature

```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_maximal_curvature
```
````

````{py:function} fast_adaptive_gaussian_blur(grid: numpy.ndarray, curvature: numpy.ndarray, cellsize: float, sigma_max: float, num_bins: int = 10) -> numpy.ndarray
:canonical: dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur

```{autodoc2-docstring} dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur
```
````

````{py:function} initialize_flat_steep_grid(mnt: numpy.ndarray, slope_threshold: float, cellsize: float) -> numpy.ndarray
:canonical: dem_lic.utils.morpho_dem.initialize_flat_steep_grid

```{autodoc2-docstring} dem_lic.utils.morpho_dem.initialize_flat_steep_grid
```
````

````{py:function} remove_small_flat_areas(flat_steep_grid: numpy.ndarray, min_area: int) -> numpy.ndarray
:canonical: dem_lic.utils.morpho_dem.remove_small_flat_areas

```{autodoc2-docstring} dem_lic.utils.morpho_dem.remove_small_flat_areas
```
````

````{py:function} calculate_local_range(mnt: numpy.ndarray, steps: int = 5) -> numpy.ndarray
:canonical: dem_lic.utils.morpho_dem.calculate_local_range

```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_local_range
```
````
