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
import numpy as np

#%%
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\quadtree\kolding_river\subgrid"
data_dir = r"p:\11209905-dca-sfincs-river\00_data\geometries"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding_quadtree.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain_kolding.geojson")

mod = SfincsModel(root=sfincs_root, mode="w+", logger=logger, **kwargs)

#%%

gdf_riv = mod.data_catalog.get_geodataframe("river_lines_kolding_single")
gdf_riv_buf = gdf_riv.assign(geometry=gdf_riv.buffer(100))
gdf_riv_buf["refinement_level"] = 3

gdf_refinement = gpd.GeoDataFrame(
    {"refinement_level": [3]},
    geometry=[
        gdf_riv_buf.unary_union,
    ],
    crs=gdf_riv.crs,
)

gdf_refinement.to_file(join(data_dir, "refinement_river.geojson"), driver="GeoJSON")

#%%

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

#%% Forcing

#Timeseries

timeseries = mod.data_catalog.get_dataframe('timeseries_discharge')

mod.config.update(
    {
        "tref": timeseries.index[0],
        "tstart": timeseries.index[1296],
        "tstop": timeseries.index[4896],
    }
)


#Discharge

mod.setup_discharge_forcing(timeseries= "river_discharges" ,
                            locations= "river_discharge_locations_kolding" )

# Waterlevel

t_sec = (timeseries.index - timeseries.index.min()).total_seconds()
htide = 0.7*np.cos(2*np.pi*(t_sec/3600 - 12.42/4) / 12.42)
htide= htide + 0.1

df_bzs = pd.DataFrame(index=timeseries.index, columns= [1], data= htide.values)

mod.setup_waterlevel_forcing(timeseries=df_bzs, locations="waterlevel_locations")
mod.setup_observation_points(locations = "observation_locations", merge = False )

mod.write()

#%%


mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
                  datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
                  datasets_riv=[{"centerlines": "river_lines_kolding_single",
                                 "mask": "river_edges_kolding", 
                                 "manning": 0.035, 
                                 "rivdph": 3,
                                 #"gdf_zb": "riverdepth_points"
                                 }],
                  write_dep_tif = True,
                  nr_subgrid_pixels = 5, 
                  nbins= 10,
                  nrmax = 2000,
                  max_gradient = 9999,

                  )

# NOTE: the waterlevel forcing data is now stored in the sf.forcing dictionary
mod.write()

# %%
