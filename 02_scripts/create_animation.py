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
from matplotlib import animation


# %%
sfincs_root = r"p:\11209905-dca-sfincs-river\01_models\KoldingA_05m_riv"
logger = setuplog("update", "./hydromt.log", log_level=10)
yml_file = join("setup_sfincs_kolding.yml")
opt = configread(yml_file, abs_path=True)  # read settings from ini file
kwargs = opt.pop("global", {})
mod = SfincsModel(root=sfincs_root, mode="r", logger=logger, **kwargs)

# with the depfile on subgrid resolution this would be:
depfile = join(sfincs_root, "subgrid", "dep_subgrid.tif")

da_dep = mod.data_catalog.get_rasterdataset(depfile)
# secondly we are reading in the model results
mod.read_results()

# thirdly, we determine the masking of the floodmap

# we could use the GSWO dataset to mask permanent water in a similar way as above
# NOTE: this is masking the water levels on the computational grid resolution

# alternatively we could use a geodataframe (e.g. OpenStreetMap landareas) to mask water bodies
# NOTE: small rivers are not masked with this geodataframe
gdf_osm = mod.data_catalog.get_geodataframe("osm_landareas")

# and again, we can use a threshold to mask minimum flood depth
hmin = 0.05

# Fourthly, we downscale the floodmap
da_hmax = utils.downscale_floodmap(
    zsmax=da_zsmax,
    dep=da_dep,
    hmin=hmin,
    gdf_mask=gdf_osm,
    floodmap_fn=join(sfincs_root, "floodmap.tif") # uncomment to save to <mod.root>/floodmap.tif
)




#%%
# now assuming we have a subgrid model, we don't have hmax available, so we are using zsmax (maximum water levels)
# compute the maximum over all time steps
da_zsmax = mod.results["zsmax"].max(dim="timemax")

# mask water depth
hmin = 0.05
da_h = mod.results["h"].copy()
da_h = da_h.where(da_h > hmin).drop("spatial_ref")
da_h.attrs.update(long_name="flood depth", unit="m")
# create hmax plot and save to mod.root/figs/sfincs_h.mp4
# requires ffmpeg install with "conda install ffmpeg -c conda-forge"

step = 1  # one frame every <step> dtout
cbar_kwargs = {"shrink": 0.6, "anchor": (0, 0)}


def update_plot(i, da_h, cax_h):
    da_hi = da_h.isel(time=i)
    t = da_hi.time.dt.strftime("%d-%B-%Y %H:%M:%S").item()
    ax.set_title(f"SFINCS water depth {t}")
    cax_h.set_array(da_hi.values.ravel())


fig, ax = mod.plot_basemap(
    fn_out=None, variable="", bmap="sat", plot_bounds=False, figsize=(11, 7)
)
cax_h = da_h.isel(time=0).plot(
    x="x", y="y", ax=ax, vmin=0, vmax=3, cmap=plt.cm.viridis, cbar_kwargs=cbar_kwargs
)
plt.close()  # to prevent double plot

ani = animation.FuncAnimation(
    fig,
    update_plot,
    frames=np.arange(0, da_h.time.size, step),
    interval=250,  # ms between frames
    fargs=(
        da_h,
        cax_h,
    ),
)

# to save to mp4
# ani.save(join(mod.root, 'figs', 'sfincs_h.mp4'), fps=4, dpi=200)

# to show in notebook:
from IPython.display import HTML

HTML(ani.to_html5_video())
# %%
