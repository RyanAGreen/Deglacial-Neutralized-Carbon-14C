import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as f
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys


# read in data
obspath = "data/observations/"
modelpath = "data/model/"


Rafter_subsurface = pd.read_excel(obspath + "prafter-2019-Gulf-CA-Data-for-Ryan.xls")
Rafter_subsurface = Rafter_subsurface.loc[
    (Rafter_subsurface["species"] == "U. peregrina")
    | (Rafter_subsurface["species"] == "Planulina ariminensis")
    | (Rafter_subsurface["species"] == "U. peregrina ")
]
Rafter_subsurface = Rafter_subsurface.sort_values(by=["calendar age [kyr BP]"])
Rafter_subsurface = Rafter_subsurface[["calendar age [kyr BP]", "D14C"]]
Rafter_subsurface = Rafter_subsurface.dropna(subset=["D14C"])
Rafter = Rafter_subsurface.groupby("calendar age [kyr BP]").mean().reset_index()
Rafter_anomalies = Rafter.rename(
    columns={"calendar age [kyr BP]": "year", "D14C": "D14CintNP"}
)

chen = pd.read_csv(obspath + "Chen2020.txt", sep="\t", header=0, skiprows=110)
Chen = chen[chen["water.depth"] == 627]
Chen = Chen.rename(columns={"cal.age": "year", "benthic.D14C": "D14CintNP"})
Chen = Chen[Chen["year"] < 20000]
Chen["year"] = Chen["year"].apply(f.fix)

# CO3 observational data
indian_obs = pd.read_table(obspath + "WIND28K_Yuetal2010.txt", sep="\t", header=None)
indian_obs = indian_obs.drop(labels=0, axis=1)
indian_obs = indian_obs.rename(columns={1: "time", 2: "CO3"})

centralpac = pd.read_table(obspath + "TTNO13PC61_Yuetal2013.txt", sep="\t", skiprows=3)
southatl = pd.read_table(obspath + "TNO57-21_Yuetal2014.txt", sep="\t", skiprows=2)


atlantic_obs = pd.read_table(obspath + "BOFS8K_Yuetal2008.txt", sep="\t", header=None)
atlantic_obs = atlantic_obs.rename(columns={0: "time", 1: "CO3"})
pacific_obs = pd.read_fwf(obspath + "GGC48_Yuetal2010.txt", sep="\t", header=None)
pacific_obs = pacific_obs.drop(labels=[0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12], axis=1)
pacific_obs = pacific_obs.rename(columns={3: "time", 13: "CO3"})

# read in model data
control = pd.read_fwf(
    NoISpath + "ForwardRun/ControlRun.txt", header=None, infer_nrows=1000
)
control = f.organizedata(control)
IScontrol = pd.read_fwf(
    ISpath + "ForwardRun/ControlRun.txt", header=None, infer_nrows=1000
)
IScontrol = f.organizedata(IScontrol)
ex1 = pd.read_table(
    NoISpath + "2Dinversion/Experiment1/Powell2Dinversion.txt", header=None, sep="\s+"
)
ex2 = pd.read_table(
    NoISpath + "2Dinversion/Experiment2/Powell2Dinversion.txt", header=None, sep="\s+"
)
ex3 = pd.read_table(
    NoISpath + "2Dinversion/Experiment3/Powell2Dinversion.txt", header=None, sep="\s+"
)
ex4 = pd.read_table(
    ISpath + "2Dinversion/Experiment3/Powell2Dinversion.txt", header=None, sep="\s+"
)
# ex4 = pd.read_table(
#     ISpath + "2Dinversion/Experiment3/Control.txt", header=None, sep="\s+"
# )

# ex3CO2 = pd.read_table(
#     NoISpath + "ForwardRun/ex3CO2only.txt", header=None, sep="\s+"
# )

ex3 = f.organizedata(ex3)
ex4 = f.organizedata(ex4)

color_anomaly = "#44AA99"
color_parallel = "#8A0C0D"
color_model = "#6699A2"

atlantic_color = "#1763F9"
indo_pac_color = "#BF840E"

# d14C observations
d14C = pd.read_csv(obspath + "IntCalSmoothed.txt", header=None)
D14C = d14C.rename(columns={0: "year", 3: "D14C"})
D14C["year"] = D14C["year"].apply(f.fix)

# Tree Ring Data
def calc_D14C(df):
    df["∆14C"] = (
        np.exp(-1 * df["r_date"] / 8033) * np.exp((1950.5 - df["t"]) / 8267) - 1
    ) * 1000
    return df


# Read the CSV file
df = pd.read_csv(
    "/home/rygreen/CY3/observations/TreeRing_IntCal20.csv",
    header=None,
    skiprows=lambda x: x in range(1, x),
)

# Identify the row indices where the separators occur
separator_rows = df.index[df.isnull().all(axis=1)].tolist()

# Initialize an empty list to store the section dataframes
section_dfs = []

# Loop through the separator rows and create separate dataframes for each section
for i in range(len(separator_rows)):
    if i == 0:
        start_idx = 0
        end_idx = separator_rows[i]
        section_df = df.iloc[start_idx:end_idx]
    else:
        start_idx = separator_rows[i - 1] + 1
        end_idx = separator_rows[i]
        section_df = df.loc[start_idx:end_idx]
    section_df.columns = section_df.iloc[0]
    section_df = section_df[1:]  # Remove first row from DataFrame
    section_df = section_df.dropna(axis=1, how="all")  # Remove empty columns
    section_df = section_df.loc[:, ["r_date", "t"]]  # Extract columns "z" and "calage"
    section_df = section_df.apply(
        pd.to_numeric, errors="coerce"
    )  # Convert "calage" column to numeric
    section_df = section_df.reset_index(drop=True)  # Reset index
    section_dfs.append(section_df)

# Concatenate all section dataframes into a single dataframe
merged_df = pd.concat(section_dfs, ignore_index=True)
merged_df = calc_D14C(merged_df)
merged_df["kyrBP"] = (1950 - merged_df["t"]) / 1000

### Hulu Cave observations
Hulu = pd.read_csv(
    "/home/rygreen/CY3/observations/HuluCaveD14C.tab", sep="\t", skiprows=15
)
Hulu = Hulu[(Hulu["Age [ka BP]"] >= 15) & (Hulu["Age [ka BP]"] <= 26)]


intcal20 = pd.read_csv("/home/rygreen/CY3/observations/INTCAL20.txt", skiprows=10)
intcal13 = pd.read_csv(
    "/home/rygreen/CY3/observations/INTCAL13.txt", skiprows=10, encoding="latin-1"
)
intcal09 = pd.read_csv("/home/rygreen/CY3/observations/INTCAL09.txt", skiprows=10)

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

# rafter_14C_deepsea.loess_se = rafter_14C_deepsea.loess_se.fillna(
#     (rafter_14C_deepsea.loess_se.shift() + rafter_14C_deepsea.loess_se.shift(-1)) / 2
# )
# rafter_14C_deepsea.loess_se = rafter_14C_deepsea.loess_se.interpolate()
t = rafter_14C_deepsea["calendar age bin (years BP)"] / 1000
# m_deep = rafter_14C_deepsea["binned mean (no outliers)"]
m_deep = rafter_14C_deepsea["loess_fit"]
upr68_deep = rafter_14C_deepsea["loess_upr68"]
lwr68_deep = rafter_14C_deepsea["loess_lwr68"]
upr95_deep = rafter_14C_deepsea["loess_upr95"]
lwr95_deep = rafter_14C_deepsea["loess_lwr95"]

# m_mid = rafter_14C_mid["binned mean (no outliers)"]
m_mid = rafter_14C_mid["loess_fit"]
upr68_mid = rafter_14C_mid["loess_upr68"]
lwr68_mid = rafter_14C_mid["loess_lwr68"]
upr95_mid = rafter_14C_mid["loess_upr95"]
lwr95_mid = rafter_14C_mid["loess_lwr95"]

fig, ax = plt.subplots(2, 1, sharex=True, dpi=300)
plt.rcParams["font.weight"] = "bold"

# thickening axis
for i in range(2):
    ax[i].tick_params(bottom=True, top=True, left=True, right=True)
    ax[i].tick_params(axis="both", direction="in", length=7, width=2, color="black")
    for axis in ["top", "left", "right", "bottom"]:
        ax[i].spines[axis].set_linewidth(3)
    ax[i].axvspan(11.6, 12.9, alpha=0.4, color="darkgray", zorder=0)
    ax[i].axvspan(14.5, 18, alpha=0.4, color="darkgray", zorder=0)


# Rafter Compilation
ax[0].plot(
    t, m_mid, linestyle="solid", color=color_parallel, lw=2, zorder=0,
)
# ax[0].fill_between(t, lwr95_mid, upr95_mid, color=color_parallel, alpha=0.25)
# ax[0].fill_between(t, lwr95_deep, upr95_deep, color=color_parallel, alpha=0.25)


ax[0].plot(
    t, m_deep, linestyle="--", color=color_parallel, lw=2, zorder=0,
)


ax[0].plot(
    intcal09["CAL BP"] / 1000,
    intcal09["Delta 14C"],
    color="#1763F9",
    ls=":",
    label="Atmospheric observations",
    lw=1.5,
    zorder=0,
)

ax[0].plot(
    intcal13["CAL BP"] / 1000,
    intcal13["Delta 14C"],
    color="#BF840E",
    ls=":",
    lw=1.5,
    zorder=0.1,
)

ax[0].plot(
    intcal20["CAL BP"] / 1000,
    intcal20["Delta 14C"],
    color="darkgray",
    ls="-",
    lw=1.5,
    zorder=0.2,
)

# ax[0].plot(
#     D14C.year,
#     D14C.D14C,
#     color="darkgray",
#     ls="-",
#     label="Atmospheric observations",
#     lw=2,
#     zorder=0,
# )

# ax[0].plot(
#     Hulu["Age [ka BP]"],
#     Hulu["Δ14C [‰]"],
#     color="#BF840E",
#     ls=":",
#     label="Hulu Cave",
#     lw=1,
#     zorder=1,
# )

ax[0].plot(
    merged_df["kyrBP"],
    merged_df["∆14C"],
    color="#627A10",
    ls="-.",
    label="tree ring data",
    lw=1,
    zorder=2,
)

# ax[0].fill_between(t, m - s, m + s, color=color_parallel, alpha=0.5)
# ax[0].plot(
#     Chen.year,
#     Chen.D14CintNP,
#     linewidth=1,
#     color=color_parallel,
#     marker="P",
#     markeredgecolor=color_parallel,
#     markerfacecolor=color_parallel,
#     ls="solid",
#     markersize=3,
#     alpha=0.5,
#     label="Galápagos",
#     zorder=2,
# )
ax[0].plot(
    Chen.year,
    Chen.D14CintNP,
    color=color_parallel,
    marker="P",
    markeredgecolor="k",
    markerfacecolor=color_parallel,
    ls=" ",
    alpha=1,
    markersize=9,
    zorder=2,
)


# ax[0].plot(ex3.year,ex3.D14CintNP,linewidth=3.5,ls='-',color = 'red',label= 'GeoCLand-Inv-P')
ax[0].plot(
    ex4.year,
    ex4.D14CintNP,
    linewidth=2.5,
    ls="-",
    color=color_model,
    label="NP+LC+PF+RC",
)
ax[0].plot(
    IScontrol.year,
    IScontrol.D14CintNP,
    linewidth=2,
    ls=":",
    color=color_model,
    label="Control run",
    zorder=2,
)

ax[0].fill_between(
    IScontrol.year,
    IScontrol.D14CintNP,
    ex4.D14CintNP,
    color=color_model,
    alpha=0.5,
    zorder=0,
)


# ax[0].legend(loc="upper left", frameon=False, fontsize="small")

# ax[1].plot(
#     ex4.year,
#     ex4.IndCO3,
#     color=color_model,
#     alpha=1,
#     linestyle=":",
#     lw=1,
#     label="Indian",
# )
# ax[1].plot(
#     IScontrol.year,
#     IScontrol.IndCO3,
#     color="black",
#     alpha=1,
#     linestyle=":",
#     lw=1,
#     label="Indian",
#     zorder=0,
# )

# ax[1].fill_between(
#     IScontrol.year, IScontrol.IndCO3, ex4.IndCO3, color="darkgray", alpha=0.7, zorder=0,
# )

ax[1].plot(
    ex4.year,
    ex4.AtlCO3,
    color=atlantic_color,
    alpha=1,
    linestyle="-",
    lw=2,
    label="Atlantic",
)
ax[1].plot(
    IScontrol.year,
    IScontrol.AtlCO3,
    color=atlantic_color,
    alpha=1,
    linestyle=":",
    lw=2,
    label="Atlantic",
    zorder=0,
)
ax[1].fill_between(
    IScontrol.year,
    IScontrol.AtlCO3,
    ex4.AtlCO3,
    color=atlantic_color,
    alpha=0.25,
    zorder=0,
)

ax[1].plot(
    ex4.year,
    (ex4.NPacCO3 + ex4.SPacCO3 + ex4.IndCO3) / 3,
    color=indo_pac_color,
    linestyle="-",
    lw=2,
    label="Indo-pacific",
)
ax[1].plot(
    IScontrol.year,
    (IScontrol.NPacCO3 + IScontrol.SPacCO3 + IScontrol.IndCO3) / 3,
    color=indo_pac_color,
    linestyle=":",
    lw=2,
    label="Indo-pacific",
    zorder=0,
)

ax[1].fill_between(
    IScontrol.year,
    (IScontrol.NPacCO3 + IScontrol.SPacCO3 + IScontrol.IndCO3) / 3,
    (ex4.NPacCO3 + ex4.SPacCO3 + ex4.IndCO3) / 3,
    color=indo_pac_color,
    alpha=0.25,
    zorder=0,
)
ax[1].plot(
    centralpac["Age-ka"],
    centralpac["CO3"],
    color=color_parallel,
    alpha=0.5,
    lw=1,
    marker="x",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    ls="None",
    markersize=5,
    zorder=0,
    label="TTNO13 PC61 Eq. Pacific",
)
# ax[1].plot(
#     southatl["Age-kyr"],
#     southatl["CO3-umol/kg"],
#     color=color_parallel,
#     alpha=0.5,
#     lw=1,
#     marker="<",
#     markeredgecolor=atlantic_color,
#     markerfacecolor=atlantic_color,
#     ls="None",
#     markersize=5,
#     zorder=0,
#     label="TNO57-21 S. Atlantic",
# )

ax[1].plot(
    atlantic_obs["time"],
    atlantic_obs["CO3"],
    color=color_parallel,
    alpha=0.5,
    lw=1,
    marker="o",
    markeredgecolor=atlantic_color,
    markerfacecolor=atlantic_color,
    ls="None",
    markersize=5,
    zorder=0,
    label="BOFS 8K N. Atlantic",
)
# ax[1].plot(
#     atlantic_obs["time"],
#     atlantic_obs["CO3"],
#     color=color_parallel,
#     alpha=0.4,
#     lw=1,
#     marker="o",
#     markeredgecolor=color_parallel,
#     markerfacecolor=color_parallel,
#     ls="None",
#     markersize=3,
#     zorder=0,
# )
ax[1].plot(
    indian_obs["time"],
    indian_obs["CO3"],
    color=color_parallel,
    alpha=0.5,
    lw=1,
    marker="^",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    ls="None",
    markersize=5,
    zorder=0,
    label="WIND 28K Indian",
)
# ax[1].plot(
#     indian_obs["time"],
#     indian_obs["CO3"],
#     color=color_parallel,
#     alpha=0.4,
#     lw=1,
#     marker="^",
#     markeredgecolor=color_parallel,
#     markerfacecolor=color_parallel,
#     ls="None",
#     markersize=5,
#     zorder=0,
# )
ax[1].plot(
    pacific_obs["time"],
    pacific_obs["CO3"],
    color=indo_pac_color,
    alpha=0.5,
    lw=1,
    marker="s",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    ls="None",
    markersize=5,
    zorder=0,
    label="GGC48 EQ Pacific",
)
# ax[1].plot(
#     pacific_obs["time"],
#     pacific_obs["CO3"],
#     color=indo_pac_color,
#     alpha=0.4,
#     lw=1,
#     marker="s",
#     markeredgecolor=indo_pac_color,
#     markerfacecolor=indo_pac_color,
#     ls="None",
#     markersize=3,
# )
# ax[1].legend(loc="upper left", frameon=False, fontsize="x-small", ncol=2)


ax[1].invert_yaxis()
ax[1].invert_yaxis()


ax[0].text(2, 330, "(a)", fontweight="bold", fontsize=10)
ax[1].text(2, 140, "(b)", fontweight="bold", fontsize=10)


# ax[0].set_title("Pacific ∆$^{14}$C", fontweight="bold", fontsize=10)
ax[1].set_title("Deep Ocean [CO$_{3}$$^{2-}$]", fontweight="bold", fontsize=10)
ax[0].set_ylabel("∆$^{14}$C (‰)", fontweight="bold", fontsize=10)
ax[1].set_ylabel("CO$_{3}$$^{2-}$ (µmol kg$^{-1}$)", fontweight="bold", fontsize=10)
ax[1].set_xlabel("Calendar age (kyr BP)", fontweight="bold", fontsize=10)
ax[1].set_xlim(0, 20)

# ax[0].legend(frameon=False,loc='upper left')

rafter_deep_legend = mlines.Line2D(
    [],
    [],
    color=color_parallel,
    linestyle="--",
    lw=2,
    label="Deep Pacific\n(Rafter et al., 2019)",
)
rafter_mid_legend = mlines.Line2D(
    [],
    [],
    color=color_parallel,
    linestyle="solid",
    lw=2,
    label="Mid-depth Pacific\n(Rafter et al., 2022)",
)


D14C_model_legend = mlines.Line2D(
    [],
    [],
    color=color_model,
    linestyle="solid",
    lw=2.5,
    label="Intermediate North Pacific\nNP+LC+PF+RC",
)

D14C_control_legend = mlines.Line2D(
    [],
    [],
    color=color_model,
    linestyle="--",
    lw=2,
    label="Intermediate North Pacific\nControl",
)

intcal20_legend = mlines.Line2D(
    [],
    [],
    color="darkgray",
    linestyle="solid",
    lw=2,
    label="IntCal20\n(Reimer et al. 2020)",
)

intcal13_legend = mlines.Line2D(
    [],
    [],
    color="#BF840E",
    linestyle=":",
    lw=2,
    label="IntCal13\n(Reimer et al. 2013)",
)

intcal09_legend = mlines.Line2D(
    [],
    [],
    color="#1763F9",
    linestyle=":",
    lw=2,
    label="IntCal09\n(Reimer et al. 2009)",
)

tree_ring_legend = mlines.Line2D(
    [],
    [],
    color="#627A10",
    linestyle="-",
    lw=2,
    label="IntCal20 tree-ring\n(Reimer et al. 2020)",
)

chen_D14C_legend = mlines.Line2D(
    [],
    [],
    color=color_parallel,
    linestyle=" ",
    marker="P",
    markeredgecolor="k",
    markerfacecolor=color_parallel,
    markersize=9,
    lw=1,
    label="Galápagos\n(Chen et al., 2020)",
)

ax[0].legend(
    handles=[
        intcal20_legend,
        tree_ring_legend,
        intcal13_legend,
        intcal09_legend,
        chen_D14C_legend,
        rafter_deep_legend,
        rafter_mid_legend,
        D14C_control_legend,
        D14C_model_legend,
    ],
    bbox_to_anchor=(1.01, 0.9),
    frameon=False,
    ncol=2,
    fontsize="x-small",
)

pacific_co3_model_legend = mlines.Line2D(
    [],
    [],
    color=indo_pac_color,
    linestyle="solid",
    lw=2,
    label="Deep Indo-Pacific\nNP+LC+PF+RC",
)

# indian_co3_model_legend = mlines.Line2D(
#     [], [], color=color_model, linestyle=":", lw=1, label="Indian",
# )

atlantic_co3_model_legend = mlines.Line2D(
    [],
    [],
    color=atlantic_color,
    linestyle="-",
    lw=2,
    label="Deep Atlantic\nNP+LC+PF+RC",
)
pacific_co3_control_legend = mlines.Line2D(
    [],
    [],
    color=indo_pac_color,
    linestyle=":",
    lw=2,
    label="Deep Indo-Pacific\nControl",
)
centralpac_legend = mlines.Line2D(
    [],
    [],
    color=indo_pac_color,
    linestyle=" ",
    marker="x",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    markersize=5,
    lw=1,
    label="Central Eq. Pacific, 4.3 km\nYTTNO13 PC61\n(Yu et al., 2013)",
)


atlantic_co3_control_legend = mlines.Line2D(
    [], [], color=atlantic_color, linestyle=":", lw=2, label="Deep Atlantic\nControl",
)


indian_co3_legend = mlines.Line2D(
    [],
    [],
    color=indo_pac_color,
    linestyle=" ",
    marker="^",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    markersize=5,
    lw=1,
    label="Indian, 4.1 km\nWIND 28K\n(Yu et al., 2010)",
)

pacific_co3_legend = mlines.Line2D(
    [],
    [],
    color=indo_pac_color,
    linestyle=" ",
    marker="s",
    markeredgecolor=indo_pac_color,
    markerfacecolor=indo_pac_color,
    markersize=5,
    lw=1,
    label="West Eq. Pacific, 3.4 km\nGGC48\n(Yu et al., 2010)",
)

atlantic_co3_legend = mlines.Line2D(
    [],
    [],
    color=atlantic_color,
    linestyle=" ",
    marker="o",
    markeredgecolor=atlantic_color,
    markerfacecolor=atlantic_color,
    markersize=5,
    lw=1,
    label="N. Atlantic, 4 km\nBOFS 8K\n(Yu et al., 2008)",
)

ax[1].legend(
    handles=[
        atlantic_co3_legend,
        indian_co3_legend,
        centralpac_legend,
        pacific_co3_legend,
        atlantic_co3_model_legend,
        atlantic_co3_control_legend,
        pacific_co3_model_legend,
        pacific_co3_control_legend,
    ],
    bbox_to_anchor=(1.5, 0.34),
    ncol=2,
    frameon=False,
    fontsize="x-small",
)
ax[0].set_ylim(-250, 450)

ax[0].text(21, 415, "observations")
ax[0].text(32, 415, "model")

# plt.tight_layout()
plt.savefig("/home/rygreen/CY3/PLOTTING/Figure3.pdf", bbox_inches="tight")
