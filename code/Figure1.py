import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cartopy
import matplotlib.ticker as mticker
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import netCDF4 as nc4

plt.rcParams["font.weight"] = "bold"

obspath = "data/observations/"

## Load in observational data ##

# hydrothermal vents
vents = pd.read_table(obspath + "hydrothermal-vent.txt", sep="\t",)
vents = vents.sort_values(by=["Latitude"], ascending=False)
EPR = vents[
    ((vents.Latitude > -15) & (vents.Latitude < 35))
    & ((vents.Longitude > -115) & (vents.Longitude < -98))
    & (vents["depth.Maximum.Or.Single.Reported"] > 100)
]

cocos = vents[
    ((vents.Latitude > 0) & (vents.Latitude < 8))
    & ((vents.Longitude > -100) & (vents.Longitude < -60))
    & (vents["depth.Maximum.Or.Single.Reported"] > 100)
]
cocos = cocos.sort_values(by=["Longitude"])

# projection of data (lat/lon)
data_proj = ccrs.PlateCarree()

# loading NetCDF file
fnm = obspath + "topo.nc"
nc = nc4.Dataset(fnm, "r")
topo = nc.variables["topo"][:]
lon = nc.variables["lon"][:]
lat = nc.variables["lat"][:]

# the 360 degree line is missing, make it cyclic
topo, lon = add_cyclic_point(topo, coord=lon)

# CO2 observations 
CO2obs = pd.read_csv(obspath + "CO2_compilation.txt", sep="\t")
CO2obs["year"] = CO2obs["year"]/1000

# ∆14C observations
d14C = pd.read_csv(obspath + "INTCAL20smoothed.txt", header=None)
D14C = d14C.rename(columns={0: "year", 3: "D14C"})
D14C["year"] = D14C["year"]/1000

# benthic ∆14C data from Rafter et al. 2018
Rafter_subsurface = pd.read_excel(obspath + "Rafter2018_benthic.xls")
Rafter_subsurface = Rafter_subsurface.loc[
    (Rafter_subsurface["species"] == "U. peregrina")
    | (Rafter_subsurface["species"] == "Planulina ariminensis")
    | (Rafter_subsurface["species"] == "U. peregrina ")
]
Rafter_subsurface = Rafter_subsurface.sort_values(by=["calendar age [kyr BP]"])
Rafter_subsurface = Rafter_subsurface[["calendar age [kyr BP]", "D14C"]]
Rafter_subsurface = Rafter_subsurface.dropna(subset=["D14C"])
Rafter_subsurface = (
    Rafter_subsurface.groupby("calendar age [kyr BP]").mean().reset_index()
)
Rafter_subsurface["calendar age [kyr BP]"] = (
    1000 * Rafter_subsurface["calendar age [kyr BP]"]
)

# planktic ∆14C data from Rafter et al 2018
Rafter_surface = pd.read_csv(obspath + "Rafter2018_planktic.tab", sep="\t", header=24)
Rafter_surface.loc[(Rafter_surface["Habitat"] == "planktic")]
Rafter_surface["Cal age [ka BP]"] = 1000 * Rafter_surface["Cal age [ka BP]"]
Rafter_surface = Rafter_surface.sort_values(by=["Cal age [ka BP]"])



# ∆14C from Chen et al 2020
chen = pd.read_csv(obspath + "Chen2020.txt", sep="\t", header=0, skiprows=110)
Chen = chen[chen["water.depth"] == 627]

# ∆14C from Marchitto et al 2007
Mar = pd.read_csv(obspath + "Marchitto2007.txt", sep="\s+")
Mar["Cal.Age"] = 1000 * Mar["Cal.Age"]

# ∆14C from Stott et al 2009
Stott = chen[chen["water.depth"] == 617]
Stott = Stott[Stott["ref."] == "Stott et al. (2009)"]

# ∆14C from Rafter et al. 2022
rafter_14C_compilation = pd.read_csv(obspath + "prafter-2022-Global-D14C-Comp-FIN.csv")
rafter_14C_deepsea = rafter_14C_compilation.loc[
    (
        rafter_14C_compilation[
            "Ocean basin and water mass (along density surfaces; see text)"
        ]
        == "PACIFIC BOTTOM"
    )
]

rafter_14C_mid = rafter_14C_compilation.loc[
    (
        rafter_14C_compilation[
            "Ocean basin and water mass (along density surfaces; see text)"
        ]
        == "PACIFIC MID"
    )
]


t = rafter_14C_deepsea["calendar age bin (years BP)"] / 1000
m_deep = rafter_14C_deepsea["loess_fit"]
upr68_deep = rafter_14C_deepsea["loess_upr68"]
lwr68_deep = rafter_14C_deepsea["loess_lwr68"]
upr95_deep = rafter_14C_deepsea["loess_upr95"]
lwr95_deep = rafter_14C_deepsea["loess_lwr95"]

m_mid = rafter_14C_mid["loess_fit"]
upr68_mid = rafter_14C_mid["loess_upr68"]
lwr68_mid = rafter_14C_mid["loess_lwr68"]
upr95_mid = rafter_14C_mid["loess_upr95"]
lwr95_mid = rafter_14C_mid["loess_lwr95"]

# Order -> 4 anoamlies (don't behave), then 3 parallel records (behave)
GoCobs = [
    Mar,
    Stott,
    Rafter_surface,
    Rafter_subsurface,
    Chen,
    rafter_14C_mid,
    rafter_14C_deepsea,
]

# make sure the columns are named the same thing
GoCobs[0] = GoCobs[0].rename(columns={"Cal.Age": "year", "D14C": "D14CintNP",})
GoCobs[1] = GoCobs[1].rename(columns={"cal.age": "year", "benthic.D14C": "D14CintNP"})
GoCobs[2] = GoCobs[2].rename(
    columns={"Cal age [ka BP]": "year", "Δ14C [‰]": "D14CintNP"}
)
GoCobs[3] = GoCobs[3].rename(
    columns={"calendar age [kyr BP]": "year", "D14C": "D14CintNP",}
)
GoCobs[4] = GoCobs[4].rename(columns={"cal.age": "year", "benthic.D14C": "D14CintNP"})
GoCobs[5] = GoCobs[5].rename(
    columns={"calendar age bin (years BP)": "year", "loess_fit": "D14CintNP",}
)
GoCobs[6] = GoCobs[6].rename(
    columns={"calendar age bin (years BP)": "year", "loess_fit": "D14CintNP",}
)


for i in range(7):
    GoCobs[i] = GoCobs[i][GoCobs[i]["year"] < 20000]
    GoCobs[i]["year"] = GoCobs[i]["year"]/1000


color_anomaly = "#BF840E"
color_parallel = "#8A0C0D"

# grabbed coordinates manualy from Chen2020.txt
Lats = [
    23.50,
    -1.22,
    22.90,
    22.90,
    -0.371,
]  # grabbed these values manually from Rafter compilation
Lons = [-111.60, -89.68, -109.50, -109.50, -90.815]


# 5 markers, none for compilations
markers = ["o", "v", "s", "<", "P"]

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(9, 11))
fig.subplots_adjust(hspace=-0.15)

# thickening axis
for axis in ["top", "left", "right", "bottom"]:
    axs[0].spines[axis].set_linewidth(3)
for axis in ["left", "right", "bottom"]:
    axs[1].spines[axis].set_linewidth(3)
# make twin axis and format figure
ax2 = axs[0].twinx()
axs[1].spines["top"].set_visible(False)
ax2.spines["bottom"].set_visible(False)
axs[1].tick_params(axis="y", labelsize="large")
axs[1].tick_params(axis="x", labelsize="large")
ax2.tick_params(axis="y", labelsize="large")
axs[0].set_yticklabels([])


# set labels
ax2.set_ylabel("CO$_2$ (ppm)", fontsize=20, fontweight="bold")
axs[1].set_ylabel("∆$^{14}$C (‰)", fontsize=20, fontweight="bold")
axs[1].set_xlabel("Calendar age (kyr BP)", fontsize=20, fontweight="bold")
axs[1].set_xlim((0, 20))
labels = [
    "Marchitto et al., 2007: Benthic",
    "Stott et al., 2009: Benthic",
    "Rafter et al., 2018: Planktic",
    "Rafter et al., 2018: Benthic",
    "Chen et al., 2020: Deep-sea coral",
    "Rafter et al., 2022: Mean mid-depth Pacific",
    "Rafter et al., 2022: Mean deep Pacific",
]
Anomalies = [GoCobs[0], GoCobs[1], GoCobs[2], GoCobs[3]]
# plot anomalies
for i in range(len(Anomalies)):
    axs[1].plot(
        Anomalies[i].year,
        Anomalies[i].D14CintNP,
        marker=markers[i],
        markeredgecolor="k",
        markerfacecolor=color_anomaly,
        color=color_anomaly,
        ls="solid",
        label=labels[i],
        lw=1,
        markersize=8,
        zorder=2,
    )

# Plot Chen et al. 2020
axs[1].plot(
    GoCobs[4].year,
    GoCobs[4].D14CintNP,
    marker=markers[4],
    markeredgecolor="k",
    markerfacecolor=color_parallel,
    color=color_parallel,
    ls="solid",
    lw=1,
    label=labels[4],
    markersize=10,
    zorder=2,
)
# Plot Rafter et al 2022
axs[1].plot(
    GoCobs[5].year,
    GoCobs[5].D14CintNP,
    linestyle="solid",
    color=color_parallel,
    label=labels[5],
    lw=2.5,
    zorder=0,
)
axs[1].fill_between(t, lwr95_mid, upr95_mid, color=color_parallel, alpha=0.25)
axs[1].fill_between(t, lwr95_deep, upr95_deep, color=color_parallel, alpha=0.25)


axs[1].plot(
    GoCobs[6].year,
    GoCobs[6].D14CintNP,
    linestyle="--",
    color=color_parallel,
    label=labels[6],
    lw=2.5,
    zorder=0,
)

# Plot atmospheric observations
axs[1].plot(
    D14C.year,
    D14C.D14C,
    color="darkgray",
    ls="-",
    label="Atmospheric observations",
    lw=4,
    zorder=0,
)
ax2.plot(
    CO2obs.year,
    CO2obs.CO2,
    color="darkgray",
    ls="-",
    label="Atmospheric CO2 observations",
    lw=4,
    zorder=0,
)


# make inset map
axins = inset_axes(
    axs[0],
    width="65%",
    height="65%",
    loc="lower left",
    axes_class=cartopy.mpl.geoaxes.GeoAxes,
    axes_kwargs=dict(map_projection=ccrs.PlateCarree()),
)

# add grey land masses
axins.add_feature(cfeature.LAND, color="lightgrey")
# and some coast lines
axins.coastlines(lw=1, color="grey")
# make colour contours with colorbar
im = axins.contourf(
    lon,
    lat,
    topo,
    levels=range(-4500, 1, 250),
    cmap="Blues_r",
    transform=data_proj,
    extend="min",
)
axins.contour(
    lon,
    lat,
    topo,
    levels=range(-4500, 1, 250),
    linewidths=0.5,
    alpha=0.5,
    colors="k",
    linestyles="-",
    transform=data_proj,
    extend="min",
)
# plot spreading centers
axins.plot(
    EPR.Longitude,
    EPR.Latitude,
    color="k",
    ls="solid",
    lw="3",
    transform=ccrs.PlateCarree(),
)
axins.plot(
    cocos.Longitude,
    cocos.Latitude,
    color="k",
    ls="solid",
    lw="3",
    transform=ccrs.PlateCarree(),
)


# limit extent of map 
axins.set_extent([-120, -80, -5, 35], crs=data_proj)

# plot where observations fall on map
for i in range(len(Anomalies)):
    axins.plot(
        Lons[i],
        Lats[i],
        marker=markers[i],
        markeredgecolor="k",
        markerfacecolor=color_anomaly,
        markersize=12,
        color=color_anomaly,
        transform=ccrs.PlateCarree(),
    )
axins.plot(
    Lons[4],
    Lats[4],
    marker=markers[4],
    markeredgecolor="k",
    markerfacecolor=color_parallel,
    markersize=12,
    color=color_parallel,
    transform=ccrs.PlateCarree(),
)
axins.plot(
    Lons[4],
    Lats[4],
    marker=markers[4],
    markeredgecolor="k",
    markerfacecolor=color_parallel,
    markersize=12,
    color=color_parallel,
    transform=ccrs.PlateCarree(),
)

# style inset map
for axis in ["top", "left", "right", "bottom"]:
    axins.spines[axis].set_linewidth(5)
g1 = axins.gridlines(
    draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.0, zorder=0
)
g1.xlocator = mticker.FixedLocator([-110, -90, -70])
g1.ylocator = mticker.FixedLocator([0, 10, 20, 30])
g1.ylabel_style = {"size": 12, "color": "black", "weight": "bold"}
g1.xlabel_style = {"size": 12, "color": "black", "weight": "bold"}
g1.left_labels = False
g1.top_labels = False

# style full plot
ax2.yaxis.set_ticks(np.arange(200, 300, 20))
axs[1].yaxis.set_ticks(np.arange(-600, 600, 200))
axs[0].set_yticks([])
axs[0].set_yticklabels([])
axs[1].tick_params(bottom=True, top=False, left=True, right=True)
axs[1].tick_params(axis="both", direction="in", length=7, width=3, color="black")
axs[0].tick_params(bottom=False, top=True, left=True, right=False)
axs[0].tick_params(axis="both", direction="in", length=7, width=3, color="black")
ax2.tick_params(bottom=True, top=False, left=True, right=True)
ax2.tick_params(axis="both", direction="in", length=7, width=3, color="black")

# set axes
ax2.set_ylim(175, 290)
axs[1].set_ylim(-700, 650)
axs[1].set_xlim(0, 20)

# make legend
axs[1].legend(
    loc="lower left", fontsize=9, ncol=1, frameon=False, borderpad=0.6
) 



plt.savefig("figures/Figure1.pdf", bbox_inches="tight")
print("The plotting for Figure 1 is complete!")
# plt.show()

