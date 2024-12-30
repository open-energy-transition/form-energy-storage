#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import pypsa
import seaborn as sns
from _helpers import configure_logging, set_scenario_config
from pypsa.statistics import get_bus_and_carrier

sns.set_theme("paper", style="whitegrid")

carrier_groups = {
    "Offshore Wind (AC)": "Offshore Wind",
    "Offshore Wind (DC)": "Offshore Wind",
    "Open-Cycle Gas": "Gas",
    "Combined-Cycle Gas": "Gas",
    "Reservoir & Dam": "Hydro",
    "Pumped Hydro Storage": "Hydro",
    "Solar":"Solar",
    "solar-hsat":"Solar",
    "solar rooftop":"Solar",
}

elec_array = [["nuclear", "None", "None", "None", "None", "None", "None"],
            ["Onshore Wind", "Offshore Wind", "Hydro", "Run of River", "Solar", "biomass", "geothermal"],
            ["coal", "lignite", "Gas", "oil","Other", "None", "None"]]

def autopct_format_inner(value, threshold=1):
    return f'{value:.1f}%' if value >= threshold else ''


def autopct_format_outer(value, threshold=3):
    return f'{value:.1f}%' if value >= threshold else ''


def filter_labels(labels, autopct_values):
    return [label if autopct_values[i] else None for i, label in enumerate(labels)]


def plot_pies(ax, elec_mix_array, colors):
    size = 0.4
    vals = np.array(elec_mix_array)

    cmap = plt.colormaps["tab20c"]
    inner_colors = cmap([1, 2, 5, 6, 9, 10])
    inner_colors = [colors[i] if i != "None" else "black" for i in list(np.array(elec_array).flatten())]

    outer_labels = ["nuclear", "VRES", "Fossil"]
    outer_colors = [colors[i] for i in outer_labels]

    # threshold to show outer label
    threshold = 2.0
    # get percentage of outer layer  
    autopct_values = [100*value/vals.sum(axis=1).sum() for value in vals.sum(axis=1)]
    # filter out percentages lower than threshold
    autopct_values = [value if value >= threshold else None for value in autopct_values]
    # outer labels that has value >= threshold 
    valid_outer_labels = filter_labels(outer_labels, autopct_values)

    _, texts, autotexts = ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
        wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2), 
        autopct=autopct_format_outer,  pctdistance=1.5, #textprops={'fontsize': 15},
        labels=valid_outer_labels)
    
    # Adjust the position of autopct labels
    for autotext, label in zip(autotexts, texts):
        x, y = label.get_position()  # Get position of corresponding wedge label
        autotext.set_position((x , y - 0.11))  # Set position of autopct label below the wedge label
        #autotext.set_fontsize(14)
        align = "right" if x < 0 else "left"
        autotext.set_horizontalalignment(align)

    all_vals = vals.flatten()
    ax.pie(
        all_vals[all_vals!=0], radius=1.15-size,
        colors=[inner_colors[i] for i in range(len(inner_colors)) if all_vals[i] != 0],
        wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2),
        autopct=autopct_format_inner, pctdistance=0.7 #textprops={'fontsize': 14}
    )
    
    # Calculate total generated electricity
    total_electricity = np.sum(vals)
    # Add total generated electricity to the center of the pie chart
    ax.text(0, 0, f"{total_electricity:.2f}"+"\nTWh", ha='center', va='center') #fontsize=14.5
    ax.set(aspect="equal")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_validation_electricity_production",
            opts="Ept",
            clusters="37c",
            ll="v1.0",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    print(snakemake.input.network)

    n = pypsa.Network(str(snakemake.input.network))
    #n.loads.carrier = "load"

    historic = pd.read_csv(
        snakemake.input.electricity_production,
        index_col=0,
        header=[0, 1],
        parse_dates=True,
    )
    subset_technologies = ["Geothermal", "Nuclear", "Biomass", "Lignite", "Oil", "Coal"]
    lowercase_technologies = [
        technology.lower() if technology in subset_technologies else technology
        for technology in historic.columns.levels[1]
    ]
    historic.columns = historic.columns.set_levels(lowercase_technologies, level=1)

    colors = n.carriers.set_index("nice_name").color.where(
        lambda s: s != "", "lightgrey"
    )
    colors = colors.combine_first(pd.Series(snakemake.config["plotting"]["tech_colors"]))
    colors["Offshore Wind"] = colors["Offshore Wind (AC)"]
    colors["Gas"] = '#a85522'
    colors["Hydro"] = colors["Reservoir & Dam"]
    colors["Other"] = "lightgray"
    colors["Fossil"] = "#A18181"
    colors["VRES"] = "#0fa101"

    if len(historic.index) > len(n.snapshots):
        if n.snapshots.inferred_freq:
            historic = historic.resample(n.snapshots.inferred_freq).sum().loc[n.snapshots]
        else:
            # if frequency is inconsistent due to segmentation
            snapshots_interval = pd.Series({n.snapshots[x]:(n.snapshots[x+1] - n.snapshots[x])/np.timedelta64(1,'h') for x in range(0,len(n.snapshots)-1)})
            snapshots_interval_last = pd.Series({n.snapshots[-1]:(pd.Timestamp(snakemake.config["snapshots"]["end"]) - n.snapshots[-1])/np.timedelta64(1,'h')})
            snapshots_interval = pd.concat([snapshots_interval, snapshots_interval_last])

            first_bin = [n.snapshots[0]]
            mid_bins = [x + pd.DateOffset(hours=snapshots_interval[x]/2) for x in n.snapshots[:-1]]
            last_bin = [n.snapshots[-1] + pd.DateOffset(hours=snapshots_interval[-1])]
            bins = pd.DatetimeIndex(first_bin + mid_bins + last_bin)

            interval_labels = pd.cut(historic.index, bins=bins, right=False)
            historic = historic.groupby(interval_labels).sum()
            historic.index = pd.to_datetime(n.snapshots)

    # Set the historic date based on the snapshot year
    historic_year = historic.index.year.unique()[0]
    if historic_year != n.snapshots.year.unique()[0]:
        historic.index = historic.index.map(lambda x: x.replace(year=n.snapshots.year.unique()[0]))

    optimized = n.statistics.energy_balance(
        groupby=get_bus_and_carrier, aggregate_time=False
    ).T
    optimized = optimized[["Generator","StorageUnit","Link"]].droplevel(0, axis=1)
    optimized = optimized.rename(columns=n.buses.country, level=0)
    optimized = optimized.rename(columns=carrier_groups, level=1)
    optimized = optimized.T.groupby(level=[0, 1]).sum().T

    # Remove all carriers originated not from a country
    optimized = optimized.loc[:, optimized.columns.get_level_values(0) != '']

    # Compare carrier where historical data are available
    optimized = optimized.loc[:, optimized.columns.get_level_values(1).isin(historic.columns.get_level_values(1))]
    optimized = optimized.mul(n.snapshot_weightings.generators, axis=0)

    data = pd.concat([historic, optimized], keys=["Historic", "Optimized"], axis=1)
    data.columns.names = ["Kind", "Country", "Carrier"]

    # revert back datetime according to historical year
    data.index = data.index.map(lambda x: x.replace(year=historic_year))

    # total production per carrier
    fig, ax = plt.subplots(figsize=(6, 6))

    df = data.T.groupby(level=["Kind", "Carrier"]).sum().T.sum().unstack().T
    df = df / 1e6  # TWh
    df.plot.barh(ax=ax, xlabel="Electricity Production [TWh]", ylabel="")
    ax.grid(axis="y")
    fig.savefig(snakemake.output.production_bar, bbox_inches="tight")
    plt.close(fig)

    # total production per carrier in pie form
    fig, ax = plt.subplots(1, 2, figsize=(10,5))

    df = data.T.groupby(level=["Kind", "Carrier"]).sum().T.sum().unstack().T
    df = df / 1e6  # TWh
    df = df.fillna(0)
    df.loc["None",:] = 0

    row = 0
    for col in ["Historic", "Optimized"]:
        elec_mix_array = [[df.loc[i,col] for i in elec_array[r]] for r in [0,1,2]]
        if col == "Optimized":
            row = 1
        plot_pies(ax[row], elec_mix_array, colors)

        ax[row].set(aspect="equal")
        ax[row].set_title(col)
        
    elec_array_flat = list(np.unique(np.array(elec_array).flatten()))
    elec_array_flat.remove("None")

    fig.legend(
        handles=[mpatches.Patch(color=colors[c], label=c) for c in ["VRES","Fossil"] + elec_array_flat],
        loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.19)
    )
    fig.savefig(snakemake.output.production_pie, bbox_inches="tight")
    plt.close(fig)

    # highest diffs

    fig, ax = plt.subplots(figsize=(6, 10))

    df = data.sum() / 1e6  # TWh
    df = df["Optimized"] - df["Historic"]
    df = df.dropna().sort_values()
    df = pd.concat([df.iloc[:5], df.iloc[-5:]])
    c = colors[df.index.get_level_values(1)]
    df.plot.barh(
        xlabel="Optimized Production - Historic Production [TWh]", ax=ax, color=c.values
    )
    ax.set_title("Strongest Deviations")
    ax.grid(axis="y")
    fig.savefig(snakemake.output.production_deviation_bar, bbox_inches="tight")
    plt.close(fig)

    # seasonal operation

    fig, axes = plt.subplots(3, 1, figsize=(9, 9))

    df = (
        data.div(n.snapshot_weightings.generators, axis=0)
        .groupby(level=["Kind", "Carrier"], axis=1)
        .sum()
        .resample("1D")
        .mean()
        .clip(lower=0)
    )
    df = df / 1e3

    order = (
        (df["Historic"].diff().abs().sum() / df["Historic"].sum()).sort_values().index
    )
    c = colors[order]
    optimized = df["Optimized"].reindex(order, axis=1, level=1)
    historical = df["Historic"].reindex(order, axis=1, level=1)

    kwargs = dict(color=c, legend=False, ylabel="Production [GW]", xlabel="")

    optimized.plot.area(ax=axes[0], **kwargs, title="Optimized")
    historical.plot.area(ax=axes[1], **kwargs, title="Historic")

    diff = optimized - historical
    diff.clip(lower=0).plot.area(
        ax=axes[2], **kwargs, title="$\Delta$ (Optimized - Historic)"
    )
    diff.clip(upper=0).plot.area(ax=axes[2], **kwargs)

    ylim = max(axes[0].get_ylim()[1],axes[1].get_ylim()[1])
    axes[0].set_ylim(top=ylim)
    axes[1].set_ylim(top=ylim)
    axes[2].set_ylim(bottom=-ylim/2, top=ylim/2)

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(
        h[::-1],
        l[::-1],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        ncol=1,
        frameon=False,
        labelspacing=1,
    )
    fig.savefig(snakemake.output.seasonal_operation_area, bbox_inches="tight")
    plt.close(fig)

    # touch file
    with open(snakemake.output.plots_touch, "a"):
        pass
