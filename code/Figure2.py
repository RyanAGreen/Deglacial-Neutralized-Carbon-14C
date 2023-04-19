import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import functions as f
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import sys
import matplotlib.patheffects as pe
    
def decompose(df):
    """decomposes the geologic carbon addition into carbon species based on ALK-to-DIC ratios"""
    # set the total amount of each
    df["PgCO2"] = 0
    df["PgHCO3"] = 0
    df["PgCO3"] = 0
    df.loc[df.alk_to_dic_ratio <= 1, "PgCO2"] = (1 - df.alk_to_dic_ratio) * (df.geologic_carbon_rate_PgCyr * 100)
    df.loc[df.alk_to_dic_ratio <= 1, "PgHCO3"] = df.alk_to_dic_ratio * (df.geologic_carbon_rate_PgCyr * 100)
    df.loc[df.alk_to_dic_ratio > 1, "PgHCO3"] = (2 - df.alk_to_dic_ratio) * (df.geologic_carbon_rate_PgCyr * 100)
    df.loc[df.alk_to_dic_ratio > 1, "PgCO3"] = (df.alk_to_dic_ratio - 1) * (df.geologic_carbon_rate_PgCyr * 100)
    df["PgCO2rate"] = df["PgCO2"] / 100
    df["PgHCO3rate"] = df["PgHCO3"] / 100
    df["PgCO3rate"] = df["PgCO3"] / 100
    # get the cum sum through time
    df["PgCO2cumsum"] = (df["PgCO2"].cumsum()).round(0)
    df["PgHCO3cumsum"] = (df["PgHCO3"].cumsum()).round(0)
    df["PgCO3cumsum"] = (df["PgCO3"].cumsum()).round()
    return df

plt.rcParams["font.weight"] = "bold"

obspath = "data/observations/"
modelpath = "data/model/"

# load in observations
# CO2 observations 
CO2obs = pd.read_csv(obspath + "CO2_compilation.txt", sep="\t")
CO2obs = CO2obs.rename(columns={"year": "year_kyrBP", "CO2": "atmospheric_CO2_ppm"})
CO2obs["year_kyrBP"] = CO2obs["year_kyrBP"]/1000
CO2obs["geologic_carbon_rate_PgCyr"] = 0
CO2obs["geologic_carbon_cumulative_PgC"] = 0
CO2obs = CO2obs.iloc[:-2]

# ∆14C observations
d14Cobs = pd.read_csv(obspath + "D14Cdata.txt", header=None, sep="\t")
D14Cobs = d14Cobs.rename(columns={0: "year_kyrBP", 1: "atmospheric_∆14C_permil"})
D14Cobs["year_kyrBP"] = D14Cobs.year_kyrBP/1000
D14Cobs["geologic_carbon_rate_PgCyr"] = 0
D14Cobs["geologic_carbon_cumulative_PgC"] = 0

# load in model output 
control = pd.read_table(modelpath + "Control.txt", sep="\s+")
control_RC = pd.read_table(modelpath + "Control+RC.txt", sep="\s+")
NP = pd.read_table(modelpath + "NP.txt", sep="\s+")
NP_LC = pd.read_table(modelpath + "NP+LC.txt", sep="\s+")
NP_LC_PF = pd.read_table(modelpath + "NP+LC+PF.txt", sep="\s+")
NP_LC_PF_RC = pd.read_table(modelpath + "NP+LC+PF+RC.txt", sep="\s+")

experiments = [NP, NP_LC, NP_LC_PF, NP_LC_PF_RC]
# Crate_colors = []

for i in range(4):
    experiments[i] = decompose(experiments[i])
 

# colors
forest_color = "#117733"
permafrost_co2_color = "#8a0c0d"
optimized_color = "#6699A2"
imposed_color = "#E69F00"
model_color = "#6699A2"
bicarb_color = "#BF840E"
carb_color = "#1763F9"


colors = ["black", "black", "#6699A2"]

obs_color = "black"
control_color = "#520120"

D14C_atm_ex1 = [D14Cobs, control, experiments[0]]
CO2_atm_ex1 = [CO2obs, control, experiments[0]]

D14C_atm_ex2 = [D14Cobs, control, experiments[1]]
CO2_atm_ex2 = [CO2obs, control, experiments[1]]

D14C_atm_ex3 = [D14Cobs, control, experiments[2]]
CO2_atm_ex3 = [CO2obs, control, experiments[2]]

D14C_atm_ex4 = [D14Cobs, control, experiments[3]]
CO2_atm_ex4 = [CO2obs, control, experiments[3]]

D14C_all = [D14C_atm_ex1, D14C_atm_ex2, D14C_atm_ex3, D14C_atm_ex4]
CO2_all = [CO2_atm_ex1, CO2_atm_ex2, CO2_atm_ex3, CO2_atm_ex4]

linestyles = [":", "-", "-"]
linewidths = ["2.5", "2.5", "2.5"]
labels = ["Observations", "Control", "2D inversion"]

# make figure
fig, ax = plt.subplots(4, 4, sharey="row", figsize=(13, 10.5))
ax = ax.flatten()
plt.rcParams["font.weight"] = "bold"

# style figure 
for i in range(4, 16):
    for axis in ["top", "left", "right", "bottom"]:
        ax[i].spines[axis].set_linewidth(2)
    ax[i].tick_params(bottom=True, top=True, left=True, right=True)
    ax[i].tick_params(axis="both", direction="in", length=7, width=1.5, color="black")
    ax[i].axvspan(11.6, 12.9, alpha=0.4, color="darkgray", zorder=0)
    ax[i].axvspan(14.5, 18, alpha=0.4, color="darkgray", zorder=0)

x_values = [0, 5, 10, 15, 20]

# make cartoons

# NP+LC
# weathering
ax[1].arrow(
    -0.2,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# CaCO3 burial
ax[1].arrow(
    0.8,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# DIC and ALK addition
ax[1].arrow(
    -0.2,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# Land C uptake
ax[1].arrow(
    0.8,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[1].text(
    0.6,
    0.81,
    "land C uptake",
    color="#8A0C0D",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.6,
)


# NP+LC
# weathering
ax[2].arrow(
    -0.2,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# CaCO3 burial
ax[2].arrow(
    0.8,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# DIC and ALK addition
ax[2].arrow(
    -0.2,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# Land C uptake
ax[2].arrow(
    0.8,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[2].text(
    0.6,
    0.81,
    "land C uptake",
    color="#8A0C0D",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.6,
)
# Permafrost CO2
ax[2].arrow(
    -0.2,
    0.5,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[2].text(
    -0.25,
    0.56,
    "permafrost CO$_2$",
    color="#8A0C0D",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.6,
)


# NP+LC+PF+RC
# weathering
ax[3].arrow(
    -0.2,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# CaCO3 burial
ax[3].arrow(
    0.8,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# DIC and ALK addition
ax[3].arrow(
    -0.2,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
# Land C uptake
ax[3].arrow(
    0.8,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[3].text(
    0.6,
    0.81,
    "land C uptake",
    color="#8A0C0D",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.6,
)
# Permafrost CO2
ax[3].arrow(
    -0.2,
    0.5,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[3].text(
    -0.25,
    0.56,
    "permafrost CO$_2$",
    color="#8A0C0D",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.6,
)

# IS change
ax[3].arrow(
    0.5,
    1.1,
    0.0,
    -0.35,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#1763F9",
    ec="#1763F9",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[3].text(
    0.55,
    1.04,
    "Increase in ∆$^{14}$C",
    color="#1763F9",
    bbox=dict(facecolor="white", edgecolor="none"),
    zorder=2.4,
)

# NP
ax[0].arrow(
    -0.2,
    0.75,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[0].arrow(
    0.8,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="k",
    ec="k",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)
ax[0].arrow(
    -0.2,
    0.25,
    0.4,
    0.0,
    width=0.02,
    head_width=0.08,
    head_length=0.1,
    fc="#8A0C0D",
    ec="#8A0C0D",
    clip_on=False,
    length_includes_head=True,
    zorder=3,
)


for i in range(4):
    ax[i].text(
        -0.25,
        0.81,
        "weathering",
        color="k",
        bbox=dict(facecolor="white", edgecolor="none"),
        zorder=2.6,
    )
    ax[i].text(
        0.65,
        0.31,
        "CaCO$_3$ burial",
        color="k",
        bbox=dict(facecolor="white", edgecolor="none"),
        zorder=2.6,
    )
    ax[i].text(
        -0.25,
        0.31,
        "DIC and ALK",
        color="#8A0C0D",
        bbox=dict(facecolor="white", edgecolor="none"),
        zorder=2.6,
    )
    for axis in ["top", "left", "right", "bottom"]:
        ax[i].spines[axis].set_linewidth(2)
    ax[i].tick_params(bottom=False, top=False, left=False, right=False)
    ax[i].set_yticklabels([])
    ax[i].set_xticklabels([])
    ax[i].set_xlim(0, 1)
    ax[i].set_ylim(0, 1)

for i in range(3):
    ax[0 + 4].plot(
        D14C_atm_ex1[i].year_kyrBP,
        D14C_atm_ex1[i]["atmospheric_∆14C_permil"],
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[1 + 4].plot(
        D14C_atm_ex2[i].year_kyrBP,
        D14C_atm_ex2[i]["atmospheric_∆14C_permil"],
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[2 + 4].plot(
        D14C_atm_ex3[i].year_kyrBP,
        D14C_atm_ex3[i]["atmospheric_∆14C_permil"],
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )

    ax[3 + 4].plot(
        D14C_atm_ex4[i].year_kyrBP,
        D14C_atm_ex4[i]["atmospheric_∆14C_permil"],
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[4 + 4].plot(
        CO2_atm_ex1[i].year_kyrBP,
        CO2_atm_ex1[i].atmospheric_CO2_ppm,
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[5 + 4].plot(
        CO2_atm_ex2[i].year_kyrBP,
        CO2_atm_ex2[i].atmospheric_CO2_ppm,
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[6 + 4].plot(
        CO2_atm_ex3[i].year_kyrBP,
        CO2_atm_ex3[i].atmospheric_CO2_ppm,
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )
    ax[7 + 4].plot(
        CO2_atm_ex4[i].year_kyrBP,
        CO2_atm_ex4[i].atmospheric_CO2_ppm,
        linewidth="2.5",
        ls=linestyles[i],
        color=colors[i],
    )


for i in range(4):

    # CO3
    ax[i + 8 + 4].bar(
        D14C_all[i][2].year_kyrBP,
        D14C_all[i][2].PgCO3rate,
        width=0.15,
        color="#1763F9",
        zorder=0.25,
    )

    # HCO3
    ax[i + 8 + 4].bar(
        D14C_all[i][2].year_kyrBP,
        D14C_all[i][2].PgHCO3rate,
        bottom=D14C_all[i][2].PgCO3rate,
        width=0.15,
        color="#BF840E",
        zorder=0.275,
    )

    # CO2
    ax[i + 8 + 4].bar(
        D14C_all[i][2].year_kyrBP,
        D14C_all[i][2].PgCO2rate,
        bottom=D14C_all[i][2].PgCO3rate + D14C_all[i][2].PgHCO3rate,
        width=0.15,
        color="#8a0c0d",
        zorder=0.25,
    )
    ax[i + 8 + 4].set_xlabel("Calendar age (kyr BP)", fontweight="bold", fontsize=10)
    ax[i + 8 + 4].set_xlim(0, 20)


for i in range(4):
    ax[i + 4].fill_between(
        D14C_all[i][0].year_kyrBP,
        D14C_all[i][0]["atmospheric_∆14C_permil"],
        D14C_all[i][2]["atmospheric_∆14C_permil"],
        color="darkgray",
        alpha=0.7,
        zorder=2,
        label="Misfit and algorithm opportunity",
    )
    ax[i + 4 + 4].fill_between(
        CO2_all[i][0].year_kyrBP,
        CO2_all[i][0].atmospheric_CO2_ppm,
        CO2_all[i][2].atmospheric_CO2_ppm,
        color="darkgray",
        alpha=0.7,
        zorder=2,
        label="Misfit/algorithm opportunity",
    )


# Forest uptake plotting
ax[9 + 4].plot(
    D14C_all[1][2].year_kyrBP,
    D14C_all[1][2].terrestrial_carbon_uptake_rate_PgCyr,
    alpha=1,
    color="#3C7A2D",
    lw=1,
    zorder=0.5,
)

ax[9 + 4].bar(
    D14C_all[1][2].year_kyrBP,
    D14C_all[1][2].terrestrial_carbon_uptake_rate_PgCyr,
    width=0.1,
    alpha=1,
    color="#3C7A2D",
    zorder=1,
    hatch="///",
)

ax[10 + 4].plot(
    D14C_all[2][2].year_kyrBP,
    D14C_all[2][2].terrestrial_carbon_uptake_rate_PgCyr,
    alpha=1,
    color="#3C7A2D",
    lw=1,
    zorder=0.5,
)

ax[10 + 4].bar(
    D14C_all[2][2].year_kyrBP,
    D14C_all[2][2].terrestrial_carbon_uptake_rate_PgCyr,
    width=0.1,
    alpha=1,
    color="#3C7A2D",
    zorder=1,
    hatch="///",
)


ax[11 + 4].bar(
    D14C_all[3][2].year_kyrBP,
    D14C_all[3][2].terrestrial_carbon_uptake_rate_PgCyr,
    width=0.1,
    alpha=1,
    color="#3C7A2D",
    zorder=1,
    hatch="///",
)
ax[11 + 4].plot(
    D14C_all[3][2].year_kyrBP,
    D14C_all[3][2].terrestrial_carbon_uptake_rate_PgCyr,
    alpha=1,
    color="#3C7A2D",
    lw=1,
    zorder=0.5,
)


# permafrost plotting
ax[10 + 4].bar(
    D14C_all[2][2].year_kyrBP,
    D14C_all[2][2].terrestrial_carbon_release_rate_PgCyr,
    width=0.1,
    alpha=1,
    color="#8a0c0d",
    zorder=1,
    hatch="///",
)
ax[10 + 4].plot(
    D14C_all[2][2].year_kyrBP,
    D14C_all[2][2].terrestrial_carbon_release_rate_PgCyr,
    alpha=1,
    color="#8a0c0d",
    lw=1,
    zorder=0.5,
)

ax[11 + 4].bar(
    D14C_all[3][2].year_kyrBP,
    D14C_all[3][2].terrestrial_carbon_release_rate_PgCyr,
    width=0.1,
    alpha=1,
    color="#8a0c0d",
    zorder=1,
    hatch="///",
)
ax[11 + 4].plot(
    D14C_all[3][2].year_kyrBP,
    D14C_all[3][2].terrestrial_carbon_release_rate_PgCyr,
    alpha=1,
    color="#8a0c0d",
    lw=1,
    zorder=0.5,
)

labels = [
    "(a)",
    "(b)",
    "(c)",
    "(d)",
    "(e)",
    "(f)",
    "(g)",
    "(h)",
    "(i)",
    "(j)",
    "(k)",
    "(l)",
]

pulses_CO2 = []
pulses_HCO3 = []
pulses_CO3 = []
for i in range(4):

    pulses_CO2.append(
        [
            int(experiments[i].PgCO2cumsum[55] - experiments[i].PgCO2cumsum[20]),
            int(experiments[i].PgCO2cumsum[100] - experiments[i].PgCO2cumsum[55]),
            int(experiments[i].PgCO2cumsum[200] - experiments[i].PgCO2cumsum[100]),
        ]
    )
    pulses_HCO3.append(
        [
            int(experiments[i].PgHCO3cumsum[55] - experiments[i].PgHCO3cumsum[20]),
            int(experiments[i].PgHCO3cumsum[100] - experiments[i].PgHCO3cumsum[55]),
            int(experiments[i].PgHCO3cumsum[200] - experiments[i].PgHCO3cumsum[100]),
        ]
    )
    pulses_CO3.append(
        [
            int(experiments[i].PgCO3cumsum[55] - experiments[i].PgCO3cumsum[20]),
            int(experiments[i].PgCO3cumsum[100] - experiments[i].PgCO3cumsum[55]),
            int(experiments[i].PgCO3cumsum[200] - experiments[i].PgCO3cumsum[100]),
        ]
    )

# pulse order is CO2,HCO3,CO3-
x_location = [15.5, 8.75, 3]
y_location = [0.8, 0.55, 0.46]
y_location = [0.85, 0.85, 0.85]
pulse_names = ["1$^{st}$", "2$^{nd}$", "3$^{rd}$"]
for i in range(4):
    for j in range(3):
        ax[i + 8 + 4].text(
            x=x_location[j],
            y=y_location[j],
            s=str(pulses_CO2[i][j]),
            color="#8a0c0d",
            path_effects=[pe.withStroke(linewidth=1, foreground="black")],
            fontsize=12,
        )
        ax[i + 8 + 4].text(
            x=x_location[j],
            y=y_location[j] - 0.142,
            s=str(pulses_HCO3[i][j]),
            color="#BF840E",
            path_effects=[pe.withStroke(linewidth=1, foreground="black")],
            fontsize=12,
        )
        ax[i + 8 + 4].text(
            x=x_location[j],
            y=y_location[j] - 0.3,
            s=str(pulses_CO3[i][j]),
            color="#1763F9",
            path_effects=[pe.withStroke(linewidth=1, foreground="black")],
            fontsize=12,
        )
        ax[i + 8 + 4].text(
            x=x_location[j], y=1, s=pulse_names[j], color="black", fontsize=12,
        )


for i in range(4):
    ax[i + 4].text(x=17.5, y=450, s=labels[i], fontweight="bold", fontsize=10)
    ax[i + 4 + 4].text(x=17.5, y=290, s=labels[i + 4], fontweight="bold", fontsize=10)
    ax[i + 8 + 4].text(x=17.5, y=1.2, s=labels[i + 8], fontweight="bold", fontsize=10)

for i in range(1, 4):
    ax[i + 8 + 4].text(
        x=6,
        y=0.185,
        s=str(int(experiments[i]["terrestrial_carbon_uptake_cumulative_PgC"][200])),
        color="#3C7A2D",
        path_effects=[pe.withStroke(linewidth=1, foreground="black")],
        fontsize=13,
    )

for i in range(2, 4):
    ax[i + 8 + 4].text(
        x=15.5,
        y=0.185,
        s=str(int(experiments[i]["terrestrial_carbon_release_cumulative_PgC"][200])),
        color="#8a0c0d",
        path_effects=[pe.withStroke(linewidth=1, foreground="black")],
        fontsize=13,
    )

# making legends
CO2 = mpatches.Patch(color="#8a0c0d", label="CO$_2$ (ALK:DIC < 0.5)")
bicarb = mpatches.Patch(color="#BF840E", label="HCO$_3$$^{-}$ (0.5 < ALK:DIC < 1.5)")
carb = mpatches.Patch(color="#1763F9", label="CO$_{3}$$^{2-}$ (ALK:DIC > 1.5)")
forestation = mpatches.Patch(
    facecolor="#3C7A2D", alpha=1, hatch="///", label="Terrestrial carbon uptake"
)
permafrost = mpatches.Patch(
    facecolor="#8a0c0d", hatch="///", label="Terrestrial carbon release"
)

control_legend = mlines.Line2D(
    [], [], color="black", linestyle="solid", lw=2.5, label="Control run"
)
observations_legend = mlines.Line2D(
    [], [], color="black", linestyle=":", lw=2.5, label="Observations"
)
model_legend = mlines.Line2D(
    [], [], color="#6699A2", linestyle="solid", lw=2.5, label="Model simulation"
)
misfit_legend = mpatches.Patch(color="darkgray", alpha=0.7, label="Model-data misfit")

ax[7 + 4].legend(
    handles=[
        observations_legend,
        control_legend,
        model_legend,
        misfit_legend,
        CO2,
        bicarb,
        carb,
        forestation,
        permafrost,
    ],
    ncol=2,
    bbox_to_anchor=(1.0, -1.3),
    frameon=False,
)
# Arrow(x, y, dx, dy, *[, width])
simulated_legend = mpatches.Patch(color="black", label="Simulated by CYCLOPS")
prescribed_legend = mpatches.Patch(color="#1763F9", label="Imposed open-system flux")
calculated_legend = mpatches.Patch(color="#8A0C0D", label="Optimized open-system flux")

ax[7 + 2].legend(
    handles=[simulated_legend, prescribed_legend, calculated_legend,],
    ncol=1,
    bbox_to_anchor=(0.7, -1.3),
    frameon=False,
)

for i in range(4, 12):
    ax[i].set_xticklabels([])
# labeling
ax[4 + 4].set_ylim(175, 305)
ax[0 + 4].set_ylim(-75, 500)
ax[8 + 4].set_xticks(x_values)
ax[4 + 4].set_ylabel("Atmospheric CO$_2$ \n (ppm)", fontweight="bold", fontsize=10)
ax[0 + 4].set_ylabel("Atmospheric ∆$^{14}$C \n (‰)", fontweight="bold", fontsize=10)
ax[8 + 4].set_ylabel("Rate of release/uptake \n (PgC)", fontweight="bold", fontsize=10)
ax[0].set_title("NP", fontweight="bold", fontsize=15, pad=16)
ax[1].set_title("NP+LC", fontweight="bold", fontsize=15, pad=16)
ax[2].set_title("NP+LC+PF", fontweight="bold", fontsize=15, pad=16)
ax[3].set_title("NP+LC+PF+RC", fontweight="bold", fontsize=15, pad=16)


fig.subplots_adjust(hspace=0.1, wspace=0.15)

for i in range(4):
    pos1 = ax[i].get_position().bounds
    newpos = [pos1[0] + 0.025, pos1[1], pos1[2] * 0.7, pos1[3] * 0.8]
    ax[i].set_position(newpos)

plt.savefig("figures/Figure2.pdf", bbox_inches="tight")
print("The plotting for Figure 2 is complete!")

