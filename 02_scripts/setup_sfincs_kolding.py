# %%
from datetime import datetime
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

# %%
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA"

logger = setuplog("update", "./hydromt.log", log_level=10)
yml_file = join("setup_sfincs_kolding.yml")
opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join("input", "sfincs_domain.geojson")

# %%
mod = SfincsModel(root=sfincs_root, mode="w+", logger=logger, **kwargs)
mod.build(region={"geom": region_fn}, opt=opt)

# %%
opt = configread(yml_file, abs_path=True)  # read settings from ini file


mod = SfincsModel(root=sfincs_root, mode="r", logger=logger, **kwargs)
mod.read()
mod.set_root(sfincs_root + "riv", mode="w+")
mod.setup_subgrid(**opt["setup_subgrid"])
mod.write()
# %%
opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})
update_model(sfincs_root, sfincs_root + "_riv", opt=opt, kwargs=kwargs)


# %%

mod.write()
