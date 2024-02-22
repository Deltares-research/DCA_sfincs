# %%
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


# %%
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA_PAK_sbg20_waterlevel_bound"

logger = setuplog("update", "./hydromt.log", log_level=10)
working_dir = os.path.abspath(".")
yml_file = join(working_dir, "setup_sfincs_kolding.yml")
opt = configread(yml_file)  # read settings from ini file
kwargs = opt.pop("global", {})

region_fn = join(working_dir, "input", "sfincs_domain.geojson")

# %%
mod = SfincsModel(root=sfincs_root, mode="w+", logger=logger, **kwargs)
mod.build(region={"geom": region_fn}, opt=opt)


# %%
