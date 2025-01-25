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
* - {py:obj}`calculate_relative_altitude <dem_lic.utils.morpho_dem.calculate_relative_altitude>`
  - ```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_relative_altitude
    :summary:
    ```
````

### API

````{py:function} calculate_maximal_curvature(dem, resolution=1)
:canonical: dem_lic.utils.morpho_dem.calculate_maximal_curvature

```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_maximal_curvature
```
````

````{py:function} fast_adaptive_gaussian_blur(grid, curvature, cellsize, sigma_max, num_bins=10)
:canonical: dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur

```{autodoc2-docstring} dem_lic.utils.morpho_dem.fast_adaptive_gaussian_blur
```
````

````{py:function} initialize_flat_steep_grid(mnt, slope_threshold)
:canonical: dem_lic.utils.morpho_dem.initialize_flat_steep_grid

```{autodoc2-docstring} dem_lic.utils.morpho_dem.initialize_flat_steep_grid
```
````

````{py:function} remove_small_flat_areas(flat_steep_grid, min_area)
:canonical: dem_lic.utils.morpho_dem.remove_small_flat_areas

```{autodoc2-docstring} dem_lic.utils.morpho_dem.remove_small_flat_areas
```
````

````{py:function} calculate_relative_altitude(mnt, window_size=40)
:canonical: dem_lic.utils.morpho_dem.calculate_relative_altitude

```{autodoc2-docstring} dem_lic.utils.morpho_dem.calculate_relative_altitude
```
````
