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
import numpy as np

#%%
sfincs_root = r"c:\projects\DCA\res25_riv_zb"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain_kolding.geojson")

mod = SfincsModel(root=sfincs_root, mode="r+", logger=logger, **kwargs)
# mod.build(region={"geom": region_fn}, opt=opt)

# #%%
# mod.setup_subgrid(datasets_dep=[{"elevtn": "DEM_5x5m"}, {"elevtn": "Bathymetry_50x50m"}],
#                   datasets_rgh=[{"manning": "manning_roughness_5x5m"}],
#                 #   datasets_riv=[{"centerlines": "river_lines_kolding_single",
#                 #                  "mask": "river_edges_kolding", 
#                 #                  "manning": 0.035, 
#                 #                  # "point_zb": "riverdepth_points",
#                 #                  "rivdph": 3}],

#                   write_dep_tif = True,
#                   nr_subgrid_pixels = 5, # 5 m
#                   nbins= 10,
#                   nrmax = 2000,
#                   )
# mod.write()


mod.read_config()
mod.read_geoms()
mod.read_subgrid()
mod.read_grid()


#%% Derive bottom depths based on the bankfull discharge

zb_powlaw = hydromt.workflows.rivers.river_depth(data= mod.data_catalog.get_geodataframe("bankfull_discharge_statistik2"),
                                    method = "powlaw",
                                    qbankfull_name= "h2",
                                    hc = 0.3985,
                                    hp = 0.45,
                                    )

#%% Turn depth to bed levels to compare 

# Extract location from geodataframe and remove outlier

gdf_powlaw = mod.data_catalog.get_geodataframe("bankfull_discharge_statistik2")
gdf_powlaw["depth"] = zb_powlaw
gdf_powlaw_sorted = gdf_powlaw.sort_values(by= "ogc_fid" ).reset_index(drop=True)

#%%

# Get two raster datasets

regular = mod.data_catalog.get_rasterdataset("DEM_5x5m")
burned_in = mod.data_catalog.get_rasterdataset(os.path.join(mod.root, "subgrid", "dep_subgrid.tif"))

regular_sample = regular.raster.sample(gdf_powlaw_sorted)
burned_in_sample = burned_in.raster.sample(gdf_powlaw_sorted)

difference = regular_sample - burned_in_sample
gdf_powlaw_sorted["zbed"] = regular_sample - gdf_powlaw.depth

#%%
# Compare absolute values

fig, ax = plt.subplots(1,1, figsize = (10,7))
regular_sample.plot()
burned_in_sample.plot()
gdf_powlaw_sorted.zbed.plot()

ax.set_ylim(-10,10)
#%%

# Compare differences

fig, ax = plt.subplots(1,1, figsize = (10,7))

difference.plot()
gdf_powlaw_sorted.depth.plot()

ax.set_ylim(-2,2)

#%% Linear regression

from scipy.optimize import curve_fit

# Function to remove outliers based on the IQR method
def remove_outliers(x, y):
    # Calculate Q1 (25th percentile) and Q3 (75th percentile) of y
    Q1 = np.percentile(y, 25)
    Q3 = np.percentile(y, 75)
    IQR = Q3 - Q1

    # Define the bounds for non-outliers
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Filter data to remove outliers
    mask = (y >= lower_bound) & (y <= upper_bound)
    return x[mask], y[mask]

# Sample data for predictor (Q) and target (h)
Q = gdf_powlaw["h2"].values
h = difference.values

# Remove outliers from Q and h
Q_filtered, h_filtered = remove_outliers(Q, h)

# Define the power law function
def power_law(Q, a, b):
    return a * Q**b

# Perform curve fitting
params, covariance = curve_fit(power_law, Q_filtered, h_filtered)

# Extract the parameters
a, b = params

print(f"Fitted parameters: a = {a:.4f}, b = {b:.4f}")

# Predict values using the fitted model
h_pred = power_law(Q_filtered, a, b)

# Plot the original data and the fitted curve
plt.figure(figsize=(7,7))
plt.scatter(h_filtered, h_pred, label="Original Data", color="red")
# plt.plot(Q, h_pred, label=f"Fitted Curve: h = {a:.4f} * Q^{b:.4f}", color="blue")
plt.xlabel("Predictor (h)")
plt.ylabel("Target (h)")
plt.legend()
plt.xlim([0.5,2])
plt.ylim([0.5,2])
plt.plot([0.5,2], [0.5,2])
plt.grid(True)

plt.show()


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
