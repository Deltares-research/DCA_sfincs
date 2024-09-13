# %%
from datetime import datetime
from os.path import basename, join

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from hydromt import DataCatalog
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel

# %%
sfincs_root =r"p:\11209905-dca-sfincs-river\01_models\discharge_anlysis_base\kolding\res_25_subgrid_5m_riv_zb"

logger = setuplog("update", "./hydromt.log", log_level=10)
yml_file = join("setup_sfincs_kolding.yml")
opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

# %%
mod = SfincsModel(root=sfincs_root, mode="r+", logger=logger, **kwargs)

#Timeseries

timeseries = mod.data_catalog.get_dataframe('timeseries_discharge')
mod.config.update(
    {
        "tref": timeseries.index[0],
        "tstart": timeseries.index[1296],
        "tstop": timeseries.index[4896],
    }
)

#%%
#Discharge

df_gues = pd.read_csv(r"p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/kolding_discharge_historic_2022_kolding_river.csv")
df_gues['date'] = pd.to_datetime(df_gues['date'], format='%Y-%m-%d')
df_gues_point0 = df_gues[df_gues["point"] == 20610]

# reindex on date
df_gues_point0 = df_gues_point0.set_index('date')
timeseries_gues = pd.DataFrame(df_gues_point0["discharge"])

#rename discahrge column to 0

timeseries_gues.columns = [0]


#%%
# Read the geojson file
geojson_path = r'p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\kolding_qpoints_v2.geojson'
geojson_data = gpd.read_file(geojson_path, crs = 25832 )
gues_location = geojson_data[geojson_data['qpointno'] == 20610]

#%%

mod.setup_discharge_forcing(timeseries= timeseries_gues ,
                            locations= gues_location, merge = False)

# Waterlevel

t_sec = (timeseries.index - timeseries.index.min()).total_seconds()
htide = 0.7*np.cos(2*np.pi*(t_sec/3600 - 12.42/4) / 12.42)
htide= htide + 0.1

df_bzs = pd.DataFrame(index=timeseries.index, columns= [1], data= htide.values)

mod.setup_waterlevel_forcing(timeseries=df_bzs, locations="waterlevel_locations")
mod.setup_observation_points(locations = "observation_locations", merge = False )

#%%
# Read .shp file with observation crosssections

cross_sections = gpd.read_file(r"p:\11209905-dca-sfincs-river\00_data\Discharge_observations\discharge_observations_crosssection.shp", crs = 25832)
mod.setup_observation_lines(locations = cross_sections, merge = False)

# NOTE: the waterlevel forcing data is now stored in the sf.forcing dictionary
mod.write()
# %%
mod.plot_basemap(fn_out= join(sfincs_root, "figs", "basemap.png"))
mod.plot_forcing()
# %%
