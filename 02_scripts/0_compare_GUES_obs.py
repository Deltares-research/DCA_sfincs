#%%

import pandas as pd
import geopandas as gpd
import numpy as np
import datetime 
import matplotlib.pyplot as plt
import os
from shapely.geometry import LineString, Point
import xarray as xr
#%%


# Read the CSV file
data = pd.read_csv(r'p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/kolding_discharge_100m_historic.csv',
                   header = None,
                   sep = ':|,',
                   usecols=[3,5,7],
                   names=['point','date','discharge'],
                   # parse_dates=['date'],
                   )

#%% parse date

data['date'] = pd.to_datetime(data['date'], format='"%Y-%m-%d"')
filtered_data = data[data['date'].dt.year == 2022]
                     
# filtered_data = data[(data["points" == list_points]]
#%%Read gejson
# Read the geojson file
geojson_path = r'p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\kolding_qpoints_v2.geojson'
geojson_data = gpd.read_file(geojson_path)
data_points = geojson_data['qpointno']

# # read .shp file
# shp_path = r'p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\kolding_qpoints.shp'
# shp_data = gpd.read_file(shp_path)
# data_points = shp_data['qpointno']

# read river centerlines
geojson_path_riv = r"p:\11209905-dca-sfincs-river\00_data\geometries\centerline_kolding_merged.geojson"
centerline_data = gpd.read_file(geojson_path_riv)

#%%
# Step 1: Read the shapefile (points) and GeoJSON (river centerline)
# Make sure your shapefile and GeoJSON are in the same projection (e.g., EPSG:4326)
points_gdf = geojson_data
river_gdf = centerline_data

# Step 2: Extract the LineString (river centerline) from the GeoJSON
river_centerline = river_gdf.geometry  # Assuming the GeoJSON has a single river centerline

# # Check if the centerline is a valid LineString
# if not isinstance(river_centerline, LineString):
#     raise ValueError("The river centerline should be a LineString.")

# Step 3: Project each point onto the river centerline (find the closest point on the line)
def project_point_to_line(point, line):
    return line.interpolate(line.project(point))

# Create an empty list to store the projected points and their distances along the line
projected_points = []
distances_along_line = []

for idx, row in points_gdf.iterrows():
    point = row.geometry
    projected_point = project_point_to_line(point, river_centerline)
    
    # Step 4: Calculate the distance of the projected point along the river centerline
    distance_along_line = river_centerline.project(projected_point)

    # Append the original point, the projected point, and the distance along the line
    projected_points.append({
        'original_point': point,
        'projected_point': projected_point,
        'distance_along_line': distance_along_line,
        'attributes': row.to_dict()  # Store all attributes of the original point
    })
    distances_along_line.append(distance_along_line[0])


# Create a DataFrame from the list of projected points

projected_points_df = pd.DataFrame(projected_points)
df_distances = pd.DataFrame(distances_along_line)
df_distances["point"] = data_points


#%%
filtered_data_loc = filtered_data[filtered_data['point'].isin(data_points)]

#%%
for i in range(len(filtered_data_loc)):
    filtered_data_loc.iloc[i,2] = filtered_data_loc.iloc[i,2].removesuffix('}')
    
#%% Reset index and save to csv
filtered_data_loc.reset_index(drop=True, inplace=True)
df_combined = pd.merge(filtered_data_loc, df_distances, how='left', on='point')
df_combined.rename(columns={0: 'distance'}, inplace=True)

df_combined.to_csv(r'p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/kolding_discharge_historic_2022_kolding_river.csv', index=False)

# %%

df_gues = pd.read_csv(r"p:/11209905-dca-sfincs-river/00_data/kolding_stream_data/kolding_discharge_historic_2022_kolding_river.csv")
df_gues['date'] = pd.to_datetime(df_gues['date'], format='%Y-%m-%d')


#%% Plot the discahrge over the river length for the GUES data

# Collect the observed data for 28-02-2022

time_int = "2022-02-25"

df_gues_peak = df_gues[df_gues['date'] == time_int]
#df_gues_grouped = df_gues_peak.groupby('point').mean()

# sort by distance

df_gues_peak = df_gues_peak.sort_values(by = 'distance')
distance_tributary_1 = 4600
distance_tributary_2 = 7350

# Points in df_gues sgould be in correct order (Upstream to downstream)

fig, ax  = plt.subplots(1,1)
df_gues_peak.plot.scatter(x = 'distance', y = 'discharge', ax = ax)
ax.axvline(x = distance_tributary_1, color = 'r', linestyle = '--')
ax.axvline(x = distance_tributary_2, color = 'r', linestyle = '--')
#ax.plot(df_gues_grouped.index, df_gues_grouped['discharge'], marker = 'o', linestyle = '--')


sfincs_his = xr.open_dataset(r'p:\11209905-dca-sfincs-river\01_models\discharge_anlysis_base\kolding\res_25_subgrid_5m_riv_zb\sfincs_his.nc')
sfincs_q = sfincs_his["crosssection_discharge"] 

#sfincs_q_t = sfincs_q.sel(time = time_int)
# resample on hours

sfincs_q_t = sfincs_q.resample(time = 'D').mean()
sfincs_q_t = sfincs_q_t.sel(time = time_int)

df_sfincs_q_t = sfincs_q_t.to_dataframe()

# add column "distance" to df_sfincs_q_t
df_sfincs_q_t['distance'] = df_gues_peak['distance'].values

#df_sfincs_q_t = pd.concat([df_sfincs_q_t ,df_gues_peak['distance']], ignore_index= True, axis = 1)

ax.scatter(df_sfincs_q_t['distance'], abs(df_sfincs_q_t['crosssection_discharge']), color = "red")

# %%

base_dir =  r'p:\11209905-dca-sfincs-river\00_data\Discharge_observations'
file_names = ['Kolding aa Ejstrup 2022.txt', "Seest Molleoe 2022.txt"] 
file_names = ["Vester nebel 2022.txt"]

for file_name in file_names:
    print(file_name)

    if file_name == 'Kolding aa Ejstrup 2022.txt':

        df_obs = pd.read_csv(os.path.join(base_dir, file_name), skiprows = 4, 
                            encoding="utf-8", header =None, parse_dates=[0],
                            dayfirst=True, names=['Date', 'Waterlevel', 'Discharge'])
        
    else:
        df_obs = pd.read_csv(os.path.join(base_dir, file_name), skiprows = 4, 
                    encoding="utf-8", header =None, parse_dates=[0], sep = ';', decimal= ',',
                    dayfirst=True, names=['Date', 'Waterlevel', 'Discharge'])

# %%

df_gues = pd.read_csv(r"p:\11209905-dca-sfincs-river\00_data\kolding_stream_data\kolding_discharge_100m_historic_2022_all.csv")
df_gues['date'] = pd.to_datetime(df_gues['date'], format='%Y-%m-%d')
date_range_obs = [df_obs['Date'].min(), df_obs['Date'].max()]
df_gues_obs = df_gues[df_gues['date'].between(date_range_obs[0], date_range_obs[1])]
df_gues_obs['discharge'] = df_gues_obs['discharge']
df_gues_obs_vester_nebel = df_gues_obs[df_gues_obs['point'] == 25463] 


fig, ax  = plt.subplots(1,1)
df_obs.plot(x = 'Date', y = 'Discharge', ax = ax, label = 'Observed')
df_gues_obs_vester_nebel.plot(x = 'date', y = 'discharge', ax = ax, label = 'GUES')
# %%
