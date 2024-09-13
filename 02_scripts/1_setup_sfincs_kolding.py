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
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\discharge_anlysis_base\kolding\res_25_subgrid_5m_riv_zb"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain_kolding.geojson")

mod = SfincsModel(root=sfincs_root, mode="w+", logger=logger, **kwargs)
mod.build(region={"geom": region_fn}, opt=opt)

#%%
mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
                  datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
                  datasets_riv=[{"centerlines": "river_lines_kolding_single",
                                 "mask": "river_edges_kolding", 
                                 "manning": 0.035, 
                                 "point_zb": "riverdepth_points"}],
                  write_dep_tif = True,
                  nr_subgrid_pixels = 5, # 5 m
                  nbins= 10,
                  nrmax = 2000,
                  )
mod.write()














#%% Derive bottom depths based on the bankfull discharge

zb_powlaw = hydromt.workflows.rivers.river_depth(data= mod.data_catalog.get_geodataframe("bankfull_discharge_statistik"),
                                    method = "powlaw",
                                    qbankfull_name= "q01",
                                    )

#%%

gdf_riv_mask = mod.data_catalog.get_geodataframe("river_edges_kolding")
gdf_riv = mod.data_catalog.get_geodataframe("river_lines_kolding")
da_elv = mod.data_catalog.get_rasterdataset("DEM_5x5m")

#%%

gdf_mask = gpd.GeoDataFrame(
        geometry=[gdf_riv_mask.buffer(0).unary_union],
        crs=gdf_riv_mask.crs,
    )  # create single polygon to clip
gdf_riv_clip = gdf_riv.overlay(gdf_mask, how="difference")
gdf_riv_clip["buf"] = mod.res[0] / 2
gdf_riv_mask1 = gdf_riv_clip.assign(
        geometry=gdf_riv_clip.buffer(gdf_riv_clip["buf"])
    )
gdf_riv_mask = gpd.overlay(gdf_riv_mask, gdf_riv_mask1, how="union")
da_riv_mask = da_elv.raster.geometry_mask(gdf_riv_mask)

#%%
## Derive river width

riv_width = hydromt.workflows.rivers.river_width(gdf_stream= gdf_riv_clip,
                                                 da_rivmask = da_riv_mask)

gdf_riv_clip["width"] = riv_width

#%%
gdf_riv_clip = gdf_riv_clip[(gdf_riv_clip["width"] > 0) & (gdf_riv_clip["width"] < 100)]

#%% From this gdf, sample the centerline for the discahrge points

# Compute centerlines center points

gdf_riv_clip["center"] = gdf_riv_clip.geometry.centroid

# Loop trough discharge points and find the closest centerline point

gdf_points = mod.data_catalog.get_geodataframe("bankfull_discharge_statistik")

#%%
gdf_points["riv_width"] = 0
for i, row in gdf_points.iterrows():
    point = row.geometry
    nearest_point = gdf_riv_clip.center.apply(lambda geom: geom.distance(point)).idxmin()
    width = gdf_riv_clip.loc[nearest_point, "width"]
    gdf_points.loc[i, "riv_width"] = width

#%%
# Flow direction

da_flwdir = hydromt.flw.d8_from_dem(da_elv = da_elv)
#da_flwdir.rio.to_raster(r"p:\11209905-dca-sfincs-river\00_data\flwdir.tif")
# get the flow direction fro tif

da_flwdir = rioxarray.open_rasterio(r"p:\11209905-dca-sfincs-river\00_data\flwdir.tif")[0]
flwdir_class = pyflwdir.from_array(da_flwdir.values, ftype="d8")

#flwdir_class = Flwdir(idxs_ds = da_flwdir)
#%%

#assign shape and ftype

# flwdir_class.shape = da_flwdir.shape
# flwdir_class.ftype = "d8"


#%%
zb_gvf = hydromt.workflows.rivers.river_depth(data= gdf_points,
                                    method = "gvf",
                                    qbankfull_name= "q10",
                                    flwdir = flwdir_class,
                                    rivwth_name = "riv_width",
                                    rivzs_name = "h2",
                                    manning = 0.035
                                    )

#%%

mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
                  datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
                  datasets_riv=[{"centerlines": "river_lines_kolding_single",
                                 "mask": "river_edges_kolding", 
                                 "manning": 0.035, 
                                 "rivdph": 3,
                                 #"point_zb": "riverdepth_points"
                                 }],
                  write_dep_tif = True,
                  nr_subgrid_pixels = 2, # 5 m
                  nbins= 10,
                  nrmax = 2000,
                  )


# %%
