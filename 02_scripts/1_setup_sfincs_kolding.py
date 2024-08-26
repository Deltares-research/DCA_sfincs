#%%
from datetime import datetime
import os
from os.path import basename, join

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from hydromt import DataCatalog
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel


sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\reduced_region\kolding_river\riv_3m"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain_reduced.geojson")

mod = SfincsModel(root=sfincs_root, mode="w+", logger=logger, **kwargs)
mod.build(region={"geom": region_fn}, opt=opt)


# mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
#                   datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
#                   datasets_riv=[{"centerlines": "river_lines_kolding_single",
#                                  "mask": "river_edges_kolding", 
#                                  "manning": 0.035, 
#                                  "gdf_zb": "riverdepth_points"}],
#                   write_dep_tif = True,
#                   nr_subgrid_pixels = 5, # 5 m
#                   nbins= 10,
#                   nrmax = 2000,
#                   )



mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
                  datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
                  datasets_riv=[{"centerlines": "river_lines_kolding_single",
                                 "mask": "river_edges_kolding", 
                                 "manning": 0.035, 
                                 "rivdph": 3}],
                  write_dep_tif = True,
                  nr_subgrid_pixels = 5, # 5 m
                  nbins= 10,
                  nrmax = 2000,
                  )


# %%
