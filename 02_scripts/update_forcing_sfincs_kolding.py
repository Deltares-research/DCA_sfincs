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
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA_05m_riv"

logger = setuplog("update", "./hydromt.log", log_level=10)
yml_file = join("setup_sfincs_kolding.yml")
opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})


# %%
mod = SfincsModel(root=sfincs_root, mode="r+", logger=logger, **kwargs)

# %%
mod.setup_river_inflow(
    rivers="river_lines_kolding",
    keep_rivers_geom=True,
)


# %%


def generate_curved_dataset(min_value, max_value, num_points):
    x = np.linspace(0, 1, num_points)
    y = (
        np.sin(x * np.pi) + 1
    ) * 0.5  # You can modify this function to create different curves

    # Scale and shift the curve to fit within the specified min and max values
    scaled_curve = y * (max_value - min_value) + min_value

    return scaled_curve


locations = mod.data_catalog.get_geodataframe(
    r"c:\Git_repos\DCA_sfincs\02_scripts\input\src_dummy.geojson", geom=mod.region
).set_crs('EPSG:25832', allow_override = True)
timeseries = pd.DataFrame(
    columns=locations.index, index=pd.date_range("20220101", "20220105", freq="D")
).fillna(0)
for column in timeseries.columns:
    timeseries[column] = generate_curved_dataset(0, locations.q01[column], 5) * 10

# %%
mod.config.update(
    {
        "tref": datetime(2022, 1, 1, 0, 0),
        "tstart": datetime(2022, 1, 1, 0, 0),
        "tstop": datetime(2022, 1, 5, 0, 0),
    }
)



#%%
mod.setup_discharge_forcing(timeseries= timeseries, locations = locations, buffer = 10000)
# %%
mod.write()
# %%
