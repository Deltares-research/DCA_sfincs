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

root = r"c:\Git_repos\DCA_sfincs"
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA"
input_file = r"p:\11209905-dca-sfincs-river\00_data\geometries\old\riverlines_kolding_area_cleaned.geojson"
output_file = r"p:\11209905-dca-sfincs-river\00_data\geometries\riverlines_kolding_area_prepared_without_width.geojson"

rivers = gpd.read_file(input_file)

# rivers = rivers.rename(columns={"MIDTBREDDE": "rivwth"})
# rivers["rivwth"].where(rivers["rivwth"] != "2.5-12", 7.25, inplace=True)
# rivers["rivwth"].where(rivers["rivwth"] != "Ukendt", 5, inplace=True)
# rivers["rivwth"].where(rivers["rivwth"] != "12-", 15, inplace=True)

rivers["rivdph"] = 1

# rivers["rivwth"] = rivers["rivwth"].astype(float)

rivers.to_file(output_file)

# %%
