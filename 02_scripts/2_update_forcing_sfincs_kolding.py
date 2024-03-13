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
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA_PAK_res25_sub5_riv3m"

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


#Discharge

mod.setup_discharge_forcing(timeseries= "river_discharges" ,
                            locations= "river_discharge_locations" )

# Waterlevel

t_sec = (timeseries.index - timeseries.index.min()).total_seconds()
htide = 0.7*np.cos(2*np.pi*(t_sec/3600 - 12.42/4) / 12.42)
htide= htide + 0.1

df_bzs = pd.DataFrame(index=timeseries.index, columns= [1], data= htide.values)

mod.setup_waterlevel_forcing(timeseries=df_bzs, locations="waterlevel_locations")
mod.setup_observation_points(locations = "observation_locations", merge = False )
# NOTE: the waterlevel forcing data is now stored in the sf.forcing dictionary
mod.write()
# %%
mod.plot_basemap()
# %%
