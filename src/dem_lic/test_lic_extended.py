# -*- coding: utf-8 -*-
"""
Test for lic_extended script modules
"""

import rasterio 
import numpy as np

import dem_lic.utils.lic_extended as lic

# File paths
MNT_input_path = "../../../QGIS/out/Valtellina_NASADEM.tif"
LIC_complete_output_path = "../../../out_scripts_test_temp/Valtellina_NASADEM_LIC_article.tif"

src_rasterio = rasterio.open(MNT_input_path)

grid = src_rasterio.read(1)
cellsize = src_rasterio.transform[0]
num_steps = 5




