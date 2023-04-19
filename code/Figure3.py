import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys

def calc_D14C(df):
    """calculate ∆14C from radiocarbon age"""
    df["∆14C"] = (
        np.exp(-1 * df["r_date"] / 8033) * np.exp((1950.5 - df["t"]) / 8267) - 1
    ) * 1000
    return df

# read in data
obspath = "data/observations/"
modelpath = "data/model/"

# read in ∆14C observational data 

# Chen et al. 2020
chen = pd.read_csv(obspath + "Chen2020.txt", sep="\t", header=0, skiprows=110)
Chen = chen[chen["water.depth"] == 627]
Chen = Chen.rename(columns={"cal.age": "year", "benthic.D14C": "D14CintNP"})
Chen = Chen[Chen["year"] < 20000]
Chen["year_kyrBP"] = Chen["year"]/1000

# Rafter et al 2022
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

# ∆14C from IntCal 
intcal20 = pd.read_csv("/home/rygreen/CY3/observations/INTCAL20.txt", skiprows=10)
intcal13 = pd.read_csv(
    "/home/rygreen/CY3/observations/INTCAL13.txt", skiprows=10, encoding="latin-1"
)
intcal09 = pd.read_csv("/home/rygreen/CY3/observations/INTCAL09.txt", skiprows=10)

# tree ring data from IntCal20
# needed extra code to parse this dataset
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
tree_ring = pd.concat(section_dfs, ignore_index=True)
tree_ring = calc_D14C(tree_ring)
tree_ring["year_kyrBP"] = (1950 - tree_ring["t"]) / 1000


# read in CO3 observational data
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

# load in model output 
control = pd.read_table(modelpath + "Control.txt", sep="\s+")
control_RC = pd.read_table(modelpath + "Control+RC.txt", sep="\s+")
NP = pd.read_table(modelpath + "NP.txt", sep="\s+")
NP_LC = pd.read_table(modelpath + "NP+LC.txt", sep="\s+")
NP_LC_PF = pd.read_table(modelpath + "NP+LC+PF.txt", sep="\s+")
NP_LC_PF_RC = pd.read_table(modelpath + "NP+LC+PF+RC.txt", sep="\s+")

color_anomaly = "#44AA99"
color_parallel = "#8A0C0D"
color_model = "#6699A2"

atlantic_color = "#1763F9"
indo_pac_color = "#BF840E"

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
ax[0].plot(
    tree_ring["year_kyrBP"],
    tree_ring["∆14C"],
    color="#627A10",
    ls="-.",
    label="tree ring data",
    lw=1,
    zorder=2,
)
ax[0].plot(
    Chen.year_kyrBP,
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
ax[0].plot(
    NP_LC_PF_RC["year_kyrBP"],
    NP_LC_PF_RC["intNP_∆14C_permil"],
    linewidth=2.5,
    ls="-",
    color=color_model,
    label="NP+LC+PF+RC",
)
ax[0].plot(
    control_RC["year_kyrBP"],
    control_RC["intNP_∆14C_permil"],
    linewidth=2,
    ls=":",
    color=color_model,
    label="Control run",
    zorder=2,
)
ax[0].fill_between(
    control_RC["year_kyrBP"],
    control_RC["intNP_∆14C_permil"],
    NP_LC_PF_RC["intNP_∆14C_permil"],
    color=color_model,
    alpha=0.5,
    zorder=0,
)
ax[1].plot(
    NP_LC_PF_RC["year_kyrBP"],
    NP_LC_PF_RC["deep_atlantic_CO3_umolkg"],
    color=atlantic_color,
    alpha=1,
    linestyle="-",
    lw=2,
    label="Atlantic",
)
ax[1].plot(
    NP_LC_PF_RC["year_kyrBP"],
    NP_LC_PF_RC["deep_atlantic_CO3_umolkg"],
    color=atlantic_color,
    alpha=1,
    linestyle=":",
    lw=2,
    label="Atlantic",
    zorder=0,
)
ax[1].fill_between(
    control_RC["year_kyrBP"],
    control_RC["deep_atlantic_CO3_umolkg"],
    NP_LC_PF_RC["deep_atlantic_CO3_umolkg"],
    color=atlantic_color,
    alpha=0.25,
    zorder=0,
)
ax[1].plot(
    NP_LC_PF_RC["year_kyrBP"],
    (NP_LC_PF_RC["deep_north_pacific_CO3_umolkg"] + NP_LC_PF_RC["deep_south_pacific_CO3_umolkg"] + NP_LC_PF_RC["deep_indian_CO3_umolkg"]) / 3,
    color=indo_pac_color,
    linestyle="-",
    lw=2,
    label="Indo-pacific",
)
ax[1].plot(
    control_RC["year_kyrBP"],
    (control_RC["deep_north_pacific_CO3_umolkg"] + control_RC["deep_south_pacific_CO3_umolkg"] + control_RC["deep_indian_CO3_umolkg"]) / 3,
    color=indo_pac_color,
    linestyle=":",
    lw=2,
    label="Indo-pacific",
    zorder=0,
)
ax[1].fill_between(
    control_RC["year_kyrBP"],
    (control_RC["deep_north_pacific_CO3_umolkg"] + control_RC["deep_south_pacific_CO3_umolkg"] + control_RC["deep_indian_CO3_umolkg"]) / 3,
    (NP_LC_PF_RC["deep_north_pacific_CO3_umolkg"] + NP_LC_PF_RC["deep_south_pacific_CO3_umolkg"] + NP_LC_PF_RC["deep_indian_CO3_umolkg"]) / 3,
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
ax[1].invert_yaxis()
ax[1].invert_yaxis()


# set axes and text
ax[1].set_title("Deep Ocean [CO$_{3}$$^{2-}$]", fontweight="bold", fontsize=10)
ax[0].set_ylabel("∆$^{14}$C (‰)", fontweight="bold", fontsize=10)
ax[1].set_ylabel("CO$_{3}$$^{2-}$ (µmol kg$^{-1}$)", fontweight="bold", fontsize=10)
ax[1].set_xlabel("Calendar age (kyr BP)", fontweight="bold", fontsize=10)
ax[1].set_xlim(0, 20)
ax[0].set_ylim(-250, 450)
ax[0].text(2, 330, "(a)", fontweight="bold", fontsize=10)
ax[1].text(2, 140, "(b)", fontweight="bold", fontsize=10)
ax[0].text(21, 415, "observations")
ax[0].text(32, 415, "model")

# make proxy legends
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

plt.savefig("figures/Figure3.pdf", bbox_inches="tight")
print("The plotting for Figure 3 is complete!")
# plt.show()
