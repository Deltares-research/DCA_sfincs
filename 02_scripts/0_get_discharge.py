#%%

import numpy as np
from datetime import datetime
import os
from os.path import basename, join

import geopandas as gpd
import hydromt
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

#%% discharge timeseries

## Get water discahrge data


base_dir =  r'p:\11209905-dca-sfincs-river\00_data\Discharge_observations'

df_all = pd.DataFrame()


file_names = ['Kolding aa Ejstrup 2022.txt', "Seest Molleoe 2022.txt", "Vester nebel 2022.txt"] 

for file_name in file_names:

    print(file_name)

    if file_name == 'Kolding aa Ejstrup 2022.txt':

        df_loc = pd.read_csv(os.path.join(base_dir, file_name), skiprows = 4, 
                            encoding="utf-8", header =None, parse_dates=[0],
                            dayfirst=True, names=['Date', 'Waterlevel', 'Discharge'])
        
    else:
        df_loc = pd.read_csv(os.path.join(base_dir, file_name), skiprows = 4, 
                    encoding="utf-8", header =None, parse_dates=[0], sep = ';', decimal= ',',
                    dayfirst=True, names=['Date', 'Waterlevel', 'Discharge'])


    ## Make columsn with seconds!!

    # Calculate the seconds since the start
    df_loc['TotalSeconds'] = (df_loc['Date'] - df_loc['Date'].min()).dt.total_seconds()
    df_loc.set_index('TotalSeconds', inplace =True, drop = True)
    df_all = pd.concat([df_all, df_loc['Waterlevel']], axis = 1)



# df_all.round(4).to_csv(os.path.join(base_dir, "pre_processed", f"discharge_all.txt" ),
#                                      header =  False, index = True)

# df_loc['Date'].to_csv(os.path.join(base_dir, "pre_processed", f"discharge_timeseries.txt" ),
#                                      header =  False, index = False)




# %%

fig, ax = plt.subplots(1,1, figsize = (10,6))
#df_all = df_all.rename(columns = {df_all.columns[0]:"Kolding aa Ejstrup", df_all.columns[1]:"Seest Molleoe"})
df_all.plot(ax = ax)
ax.legend(["Kolding aa Ejstrup", "Seest Molleoe", "Vester nebel"])
ax.grid()
ax.set_ylabel("Waterlevel [m]")
ax.set_xlabel("Time [s]")

# %%
