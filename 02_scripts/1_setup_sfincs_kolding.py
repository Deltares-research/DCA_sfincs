#%%
from datetime import datetime
import os
from os.path import basename, join

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import rioxarray
import shapely
from hydromt import DataCatalog
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel
import pyflwdir 

#%%
sfincs_root = r"c:\projects\DCA\res25_riv_Q"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain_kolding.geojson")

mod = SfincsModel(root=sfincs_root, mode="r+", logger=logger, **kwargs)

#%%
mod.build(region={"geom": region_fn}, opt=opt)

#%%
zb_powlaw = hydromt.workflows.rivers.river_depth(data= mod.data_catalog.get_geodataframe("bankfull_discharge_statistik2"),
                                    method = "powlaw",
                                    qbankfull_name= "h2",
                                    hc = 0.3985,
                                    hp = 0.4500,)

gdf_powlaw = mod.data_catalog.get_geodataframe("bankfull_discharge_statistik2")
raster_regular = mod.data_catalog.get_rasterdataset("DEM_5x5m")
gdf_powlaw["zs"] = raster_regular.raster.sample(gdf_powlaw)
gdf_powlaw["rivbed"] = gdf_powlaw["zs"] - zb_powlaw
gdf_powlaw = gdf_powlaw[abs(gdf_powlaw["rivbed"])<3]

#%%
mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
                  datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
                  datasets_riv=[{"centerlines": "river_lines_kolding_single",
                                 "mask": "river_edges_kolding", 
                                 "manning": 0.035, 
                                 "point_zb": gdf_powlaw,
                                 # "rivdph": 3
                                 }],

                  write_dep_tif = True,
                  nr_subgrid_pixels = 5, # 5 m
                  nbins= 10,
                  nrmax = 2000,
                  )
mod.write()


# %%
