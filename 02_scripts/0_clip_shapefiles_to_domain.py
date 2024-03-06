# %%
from datetime import datetime
from os.path import basename, join
from pathlib import Path

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import xarray as xr
from hydromt import DataCatalog
from hydromt.config import configread
from hydromt.log import setuplog
from hydromt_sfincs import SfincsModel

# %%
input_shape_fn = r"p:\11209905-dca-sfincs-river\00_data\vanSarah_07_12_2023\GIS\Statistik1990-2019.shp"
model_domain_fn = (
    r"c:\Git_repos\DCA-sfincs-river\02_scripts\input\sfincs_domain.geojson"
)
output_folder = r"p:\11209905-dca-sfincs-river\00_data\geometries"


# %%
def clip_shape(input_fn, domain_fn, write=False):
    input_shape = gpd.read_file(input_shape_fn)
    model_domain = gpd.read_file(model_domain_fn)
    shape_clipped = gpd.clip(input_shape, model_domain)
    if write:
        shape_clipped.to_file(
            join(output_folder, f"{Path(input_shape_fn).stem}_clipped.geojson")
        )
    return shape_clipped


# %%
input_shape_fn = r"p:\11209905-dca-sfincs-river\00_data\vanSarah_07_12_2023\GIS\Statistik1990-2019.shp"
model_domain_fn = (
    r"c:\Git_repos\DCA-sfincs-river\02_scripts\input\sfincs_domain.geojson"
)
clip_shape(input_shape_fn, model_domain_fn, write=True)
# %%
input_shape_fn = r"p:\11209905-dca-sfincs-river\00_data\vanSarah_07_12_2023\GIS\riverlines_kolding_area.shp"
clipped = clip_shape(input_shape_fn, model_domain_fn, write=False)

# %%
clipped = clipped.set_crs(epsg=25832, allow_override=True)
clipped_filtered = clipped.query(
    "MIDTBREDDE == '2.5-12' or MIDTBREDDE == 'Ukendt' or MIDTBREDDE == '12-'"
)
clipped_filtered.to_file(
    join(output_folder, f"{Path(input_shape_fn).stem}_filtered.geojson")
)
# %%
clipped_filtered = clipped.query("MIDTBREDDE == 'Ukendt'")
clipped_filtered.to_file(
    join(output_folder, f"{Path(input_shape_fn).stem}_filtered_ukendt.geojson")
)

# %%
clipped_filtered = clipped.query("MIDTBREDDE == '2.5-12'")
clipped_filtered.to_file(
    join(output_folder, f"{Path(input_shape_fn).stem}_filtered_25_12.geojson")
)
# %%
gdf = clipped_filtered.copy()
main_river_id = 1094333445

# Step 1: Identify Main River
main_river = gdf[gdf["FOT_ID"] == main_river_id]

# Step 2: Determine Connectivity using NetworkX
graph = nx.Graph()

# Add edges from river segments
for line in gdf["geometry"]:
    coordinates = list(line.coords)
    graph.add_edges_from(zip(coordinates, coordinates[1:]))

# Step 3: Find Connected Components
connected_components = list(nx.connected_components(graph))

# Filter out segments not connected to the main river
connected_to_main_river = set()
for component in connected_components:
    if main_river_id in component:
        connected_to_main_river.update(component)

filtered_gdf = gdf[
    gdf["geometry"].apply(
        lambda x: any(coord in connected_to_main_river for coord in x.coords)
    )
]


# %%
