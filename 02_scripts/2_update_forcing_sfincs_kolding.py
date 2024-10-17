# %%
from datetime import datetime
from os.path import basename, join

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import os
from hydromt import DataCatalog
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel

# %%
sfincs_base = r"p:\11209905-dca-sfincs-river\01_models\discharge_anlysis_base\kolding\no_tide_tributaries_base_case_large_crosssections"

names = ["res100_riv_zb",  "res25_riv_zb", "res25_riv_3m", "res25_noburned_rivers", "res25_nosubgrid"]
long_names = ["100m", "25m", "25m - RivDepth 3m", "No burned in rivers", "No subgrid"]
colors_range = np.linspace(0, 1, len(names))

start_date = "2022-02-24"
end_date = "2022-02-27"

for i, name in enumerate(names):
    sfincs_root = os.path.join(sfincs_base, name) 

    print(name)

    logger = setuplog("update", "./hydromt.log", log_level=10)
    yml_file = join("setup_sfincs_kolding.yml")
    opt = configread(yml_file, abs_path=True)  # read settings from ini file
    kwargs = opt.pop("global", {})

    mod = SfincsModel(root=sfincs_root, mode="r+", logger=logger, **kwargs)
    #mod.read_grid()

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

    df_gues = pd.read_csv(r"p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/kolding_discharge_historic_2022_kolding_river.csv")
    df_gues['date'] = pd.to_datetime(df_gues['date'], format='%Y-%m-%d')
    df_gues_point0 = df_gues[df_gues["point"] == 20610]

    # reindex on date
    df_gues_point0 = df_gues_point0.set_index('date')
    timeseries_gues = pd.DataFrame(df_gues_point0["discharge"])

    #rename discahrge column to 0

    timeseries_gues.columns = [0]

    # Read the geojson file
    geojson_path = r'p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\kolding_qpoints_v2.geojson'
    geojson_data = gpd.read_file(geojson_path, crs = 25832 )
    gues_location = geojson_data[geojson_data['qpointno'] == 20610]


    #Discharge tributaries

    df_gues = pd.read_csv(r"p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/tribut_discharge_historic_2022_kolding_river.csv")
    df_gues['date'] = pd.to_datetime(df_gues['date'], format='%Y-%m-%d')
    df_gues_point1 = df_gues[df_gues["point"] == 20729].set_index('date')
    df_gues_point1 = pd.DataFrame(df_gues_point1["discharge"])
    df_gues_point1.columns = [1]

    df_gues_point2 = df_gues[df_gues["point"] == 25500].set_index('date')
    df_gues_point2 = pd.DataFrame(df_gues_point2["discharge"])
    df_gues_point2.columns = [2]

    # reindex on date
    timeseries_tribut = pd.concat([df_gues_point1, df_gues_point2], axis = 1)


    geojson_path = r'p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\tribut_qpoints.geojson'
    geojson_data = gpd.read_file(geojson_path, crs = 25832 )
    tribut_location = geojson_data[(geojson_data['qpointno'] == 20729) | (geojson_data['qpointno'] == 25500)]


    timeseries_discharge = pd.concat([timeseries_gues, df_gues_point1, df_gues_point2], axis = 1)
    # merge the geodataframes
    locations_discharge =  pd.concat([gues_location, tribut_location]).reset_index(drop = True)


    mod.setup_discharge_forcing(timeseries= timeseries_discharge,
                                locations= locations_discharge, merge = False)


    # Waterlevel

    t_sec = (timeseries.index - timeseries.index.min()).total_seconds()
    htide = 0.7*np.cos(2*np.pi*(t_sec/3600 - 12.42/4) / 12.42)
    htide= htide + 0.1
    htide = np.zeros(len(timeseries.index))

    df_bzs = pd.DataFrame(index=timeseries.index, columns= [1], data= htide)

    mod.setup_waterlevel_forcing(timeseries=df_bzs, locations="waterlevel_locations", merge = False)
    mod.setup_observation_points(locations = "observation_locations", merge = False )


    # Read .shp file with observation crosssections

    cross_sections = gpd.read_file(r"p:\11209905-dca-sfincs-river\00_data\Discharge_observations\discharge_observations_crosssection3.geojson", crs = 32632)
    mod.setup_observation_lines(locations = cross_sections, merge = False)

    # NOTE: the waterlevel forcing data is now stored in the sf.forcing dictionary
    mod.write()
# %%
mod.plot_basemap(fn_out= join(sfincs_root, "figs", "basemap.png"))
mod.plot_forcing()
# %%
