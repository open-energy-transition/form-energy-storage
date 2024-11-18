# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates plots for KPIs for the Project "The Role of Energy Storage in Germany's Path to a 100% Clean Power System"
with Form Energy and Open Energy Transition.
"""

import logging
import os

import numpy as np

logger = logging.getLogger(__name__)

from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from matplotlib import rc
from plot_power_network import load_projection
from matplotlib.colors import Normalize
from pypsa.statistics import get_bus_and_carrier
from plot_summary import preferred_order
from pypsa.statistics import get_bus_and_carrier_and_bus_carrier

# activate latex text rendering
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'sans-serif': ['Computer Modern Sans serif']})
plt.style.use(["ggplot"])


def plot_curtailment(network, regions, path, show_fig=True, focus_de=True):
    n = network.copy()
    map_opts = map_opts_params.copy()

    # curtailment = n.statistics.curtailment(groupby=n.statistics.groupers.get_bus_and_carrier).droplevel(0).loc[regions.index,:,:]
    elec_generation = n.generators.carrier.replace(pretty_gen).unique()
    curtailment = n.statistics.curtailment(groupby=n.statistics.groupers.get_bus_and_carrier).droplevel(0)
    elec_i = pd.Index(set(curtailment.index.get_level_values(1)).intersection(elec_generation))
    curtailment_elec = curtailment.loc[:, elec_i, :].reset_index()
    curtailment_elec = curtailment_elec.rename(index=curtailment_elec.bus.map(n.buses.location))
    curtailment_elec["bus"] = curtailment_elec.index
    curtailment_elec = curtailment_elec.groupby(["bus", "carrier"]).sum()
    curtailment_elec = curtailment_elec.loc[regions.index, 0].div(1e3)  # GW
    electricity_price = n.buses_t.marginal_price.loc[:, regions.index].mean(axis=0)
    regions["elec_price"] = (
        electricity_price
    )

    bus_size_factor = 3e5

    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    if focus_de:
        curtailment_elec = curtailment_elec.filter(like="DE", axis=0)
        regions = regions.filter(like="DE", axis=0)
        for c in n.iterate_components(n.branch_components):
            c.df.drop(c.df.index[~((c.df.bus0.str.startswith("DE")) | (c.df.bus1.str.startswith("DE")))], inplace=True)
        map_opts["boundaries"] = [4, 17, 46, 56]
        bus_size_factor = 2e5

    bus_sizes = curtailment_elec / bus_size_factor

    proj = ccrs.EqualEarth()
    regions = regions.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})

    n.plot(
        geomap=True,
        bus_sizes=bus_sizes,
        bus_colors=tech_colors,
        link_widths=0,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    regions.plot(
        ax=ax,
        column="elec_price",
        cmap="Greens",
        linewidths=0,
        legend=True,
        vmax=regions.elec_price.max(),
        vmin=regions.elec_price.min(),
        legend_kwds={
            "label": "Avg. Electricity Price [EUR/MWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    legend_x = -0.37
    legend_y = 0.57
    sizes = [6, 3, 1]  # GW

    labels = [f"{s} GWh" for s in sizes]
    sizes = [s * 1e3 / bus_size_factor for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(legend_x, 1.01),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    carriers = curtailment_elec.index.get_level_values(1).unique()
    colors = [tech_colors[c] for c in carriers]
    labels = list(carriers)

    labels = list(pd.Series(labels).replace(pretty_names))

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(legend_x, legend_y + 0.21),
        frameon=False,
        title=r"\textbf{Curtailment}",
        alignment="left",
    )

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw
    )

    ax.set_facecolor("white")

    if show_fig:
        fig.show()
    fig.savefig(path, bbox_inches="tight")

def plot_line_loading(network, regions, path, focus_de=True, value="mean", show_fig=True):

    n = network.copy()
    map_opts = map_opts_params.copy()

    def filter_connections(l, l_t, opt_name, query_func=False):
        if query_func:
            hide_connection = l.query(query_func).index
        else:
            hide_connection = []
    
        p0 = l_t.p0.copy()
        opt = l[opt_name].copy()
        alpha = pd.Series(index=l.index)
        
        p0.loc[:, hide_connection] = 0
        opt.loc[hide_connection] = 0
        alpha.loc[hide_connection] = 0

        alpha = alpha.fillna(1)
    
        return p0, opt, alpha

    if focus_de:
        n.lines["country0"] = [n.buses.country[busname] for busname in n.lines.bus0]
        n.lines["country1"] = [n.buses.country[busname] for busname in n.lines.bus1]
        n.links["country0"] = [n.buses.country[busname] for busname in n.links.bus0]
        n.links["country1"] = [n.buses.country[busname] for busname in n.links.bus1]
        
        query_lines = "country0 != 'DE' or country1 != 'DE'"
        query_links = "country0 != 'DE' or country1 != 'DE' or carrier != 'DC'"
    else:
        query_lines = False
        query_links = "carrier != 'DC'"
        
    lines_p0, lines_opt, lines_alpha = filter_connections(n.lines, n.lines_t,"s_nom_opt", query_lines)
    links_p0, links_opt, links_alpha = filter_connections(n.links, n.links_t,"p_nom_opt", query_links)

    lines_loading = abs(lines_p0)/lines_opt
    links_loading = abs(links_p0)/links_opt
    
    if value == "mean":
        line_colors = lines_loading.mean() * 100
        link_colors = links_loading.mean() * 100
        title = "Average"
        
    elif value == "max":
        line_colors = lines_loading.max() * 100
        link_colors = links_loading.max() * 100
        title = "Average max"

    line_weight = pd.DataFrame(index=n.lines.index)
    line_weight["w"] = [lines_opt[i] * n.lines.length[i] / (lines_opt * n.lines.length).sum() for i in n.lines.index]
    line_ave = round((line_weight["w"] * line_colors).sum(),1)

    link_weight = pd.DataFrame(index=n.links.index)
    link_weight["w"] = [links_opt[i] * n.links.length[i] / (links_opt * n.links.length).sum() for i in n.links.index]
    link_ave = round((link_weight["w"] * link_colors).sum(),1)

    proj = load_projection(snakemake.params.plotting)

    map_opts = snakemake.params.plotting["map"]

    if focus_de:
        map_opts["boundaries"] = [4, 17, 46, 56]

    if map_opts["boundaries"] is None:
        map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

    fig, ax = plt.subplots(subplot_kw={"projection": proj}, figsize=(7, 6))
    norm = Normalize(vmin=0, vmax=100, clip=True)

    collection = n.plot(
        ax=ax,
        line_colors=line_colors,
        line_cmap=plt.cm.jet,
        line_norm=norm,
        line_alpha=lines_alpha,
        link_colors=link_colors,
        link_cmap=plt.cm.jet,
        link_alpha=links_alpha,
        link_norm=norm,
        title="Line loading",
        bus_sizes=1e-3,
        bus_alpha=0.7,
        boundaries=map_opts["boundaries"]
    )

    ax.set_title(f'{title} AC line loading: {line_ave} %\n{title} DC link loading: {link_ave} %', 
                loc='left', 
                x=0.02, 
                y=0.90,
                fontsize=12
                )

    plt.colorbar(collection[1], label="line loading [%]", ax=ax)

    if show_fig:
        fig.show()
    fig.savefig(path, bbox_inches="tight")

def plot_power_time_series(network, start_date, end_date, path, focus_component=["Generator","StorageUnit","Link"], focus_de=True, show_fig=True):
    n = network.copy()
    
    optimized = n.statistics.energy_balance(groupby=get_bus_and_carrier, aggregate_time=False).T
    optimized = optimized[focus_component].droplevel(0, axis=1).T
    optimized = optimized.loc[optimized.index.get_level_values('bus').isin(n.buses.query("carrier == 'AC'").index)]
    optimized = optimized.loc[optimized.index.get_level_values('carrier') != 'DC'].T
    optimized = optimized.rename(columns=n.buses.country, level=0)
    optimized = optimized.rename(columns=pretty_gen, level=1)
    optimized = optimized.T.groupby(level=[0, 1]).sum().T

    if focus_de:
        df_carrier = optimized["DE"]
    else:
        df_carrier = optimized.T.groupby(optimized.columns.get_level_values(1)).sum().T

    #Set color
    colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "green")
    
    for c in df_carrier.columns:
        if (df_carrier[c] < 0).values.any():
            df_carrier[c + " Discharge"] = df_carrier[c].clip(lower=0)
            colors[c + " Discharge"] = colors[c]
            df_carrier[c + " Charge"] = df_carrier[c].clip(upper=0)
            colors[c + " Charge"] = colors[c]
            df_carrier = df_carrier.drop(columns=c)

    order = ((df_carrier.diff().abs().sum() / df_carrier.sum()).sort_values().index)
    df_carrier = df_carrier[order]
    
    #cut into the specific timesteps
    df_carrier = df_carrier[pd.Timestamp(start_date):pd.Timestamp(end_date)]

    #kW to MW
    df_carrier = df_carrier/1e3

    #remove if smaller than 1 MWh
    df_carrier = df_carrier.loc[:,abs(df_carrier.sum()) > 1]
    
    fig, axes = plt.subplots(figsize=(12,5))
    
    if "Link" in focus_component:
        df_plot = df_carrier.drop('electricity distribution grid Charge', axis=1)
        df_plot.plot.area(ax = axes, legend=False, color = [colors[c] for c in df_plot.columns])
        abs(df_carrier['electricity distribution grid Charge']).plot(ax = axes, color = "black", linestyle='dashed')
    else:
        df_carrier.plot.area(ax = axes, legend=False, color = [colors[c] for c in df_carrier.columns])

    df_legend = pd.DataFrame()
    df_legend["handle"], df_legend["label"] = axes.get_legend_handles_labels()
    df_legend["label"] = [w.replace(' Charge', '').replace(' Discharge', '') for w in df_legend["label"]]
    df_legend["label"] = [w.replace(' &', '').replace(' and', '') for w in df_legend["label"]] #NOTE: Latex hates '&' strings because its their seperator
    df_legend = df_legend.drop_duplicates(subset = ["label"],keep = 'first').iloc[::-1]

    axes.legend(df_legend["handle"], df_legend["label"], loc = "upper center", bbox_to_anchor = (0.5, -0.15), frameon = False, ncol = 4, 
                title = "Components", title_fontproperties = {'weight':'bold'})

    axes.grid(axis="y")
    axes.set_ylabel("Energy Balance [MW]")
    axes.set_xlabel("")
    axes.set_facecolor("white")

    if show_fig:
        fig.show()
    fig.savefig(path, bbox_inches="tight")

def find_electricity_carrier(n):
    df = pd.DataFrame(n.statistics.energy_balance(groupby=get_bus_and_carrier_and_bus_carrier))
    df = df.groupby(["carrier","bus_carrier"]).sum().reset_index().drop(columns=0)
    df = df[df["bus_carrier"].isin(["AC"])]
    return df["carrier"].unique()

def plot_system_cost(network, nodal_costs, path, focus_component=["Generator","StorageUnit","Link","Store"], focus_de=True, electricity_only=True, show_fig=False):
    n = network.copy()

    component_nice_name = {"Generator":'generators',
                           "StorageUnit":'storage_units',
                           "Link":'links',
                           "Store":'stores'
                          }

    cost_df = pd.read_csv(nodal_costs, index_col=list(range(3)), header=list(range(4)))
    cost_df = cost_df.set_axis(['carrier', 'costs'], axis=1)
    cost_df = cost_df.rename(index=n.buses.country, level=2)

    cost_df = cost_df[cost_df.index.get_level_values(0).isin(component_nice_name[comp] for comp in focus_component)]

    if focus_de:
        cost_df = cost_df[cost_df.index.get_level_values(2) == 'DE']

    df = cost_df.groupby(["carrier"]).sum()

    # convert to billions
    df = df / 1e9

    df["nice_name"] = list(pd.Series(df.index).replace(n.carriers.nice_name))
    df = df.set_index("nice_name")

    if electricity_only:
        elec_carrier = find_electricity_carrier(n)
        df = df[df.index.isin(elec_carrier)]

        electricity_title = "electricity only "
    else:
        electricity_title = ""

    to_drop = df.index[df.max(axis=1) < snakemake.config["plotting"]["costs_threshold"]] #snakemake.params.plotting["costs_threshold"]]

    logger.info(
        f"Dropping technology with costs below {snakemake.params['plotting']['costs_threshold']} EUR billion per year"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    logger.info(f"Total {electricity_title}system cost of {round(df.sum().iloc[0])} EUR billion per year")

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    #Set color
    colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "green")

    fig, ax = plt.subplots(figsize=(8, 8))

    df.loc[new_index].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[colors[i] for i in new_index],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.set_ylim([0, snakemake.params.plotting["costs_max"]])

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="x")

    labels = [label + ": \n " + str(round(df.loc[label,"costs"],2)) for label in labels]
    labels = [label.replace(' &', '').replace(' and', '') for label in labels] #NOTE: Latex hates '&' strings because its their seperator

    ax.legend(
        handles, labels, ncol=1, loc="center left", bbox_to_anchor=[1, 0.5], frameon=False, 
        title = f"{electricity_title}system cost: \n{round(df.sum().iloc[0],2)}", alignment="left"
    )

    ax.set_facecolor("white")

    if show_fig:
        fig.show()
    fig.savefig(path, bbox_inches="tight")

if __name__ == "__main__":

    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_KPIs",
            opts="",
            clusters="52",
            ll="v1.11",
            sector_opts="",
            planning_horizons="2035",
            configfiles="config/config.form.yaml",
        )

    tech_colors_custom = {
        "H2 Electrolysis": "#ff29d9",
        "SMR": "#f073da",
        "SMR CC": "#c251ae",
        "Fischer-Tropsch": "#25c49a",
        "H2 Fuel Cell": "#2d8077",
        "H2 for industry": "#cd4f41",
        "H2 for shipping": "#238fc4",
        "OCGT H2 retrofitted": "#1c404c",
        "Sabatier": "#de9e46",
        "methanation": "#de9e46",
        "land transport fuel cell": "#c6dfa2",
        "methanolisation": "#edba1c",
        "Iron-Air Battery Storage": "#edba1c",
        "Li-Ion Battery Storage": "#1c404c",
        "Reservoir \& Dam": '#298c81',
        "Solar Thermal": '#d7b24c',
        "Solar (HSAT)": "#fdb915",
    }

    pretty_gen = {
        "onwind": "Onshore Wind",
        "offwind-ac": "Offshore Wind (AC)",
        "offwind-dc": "Offshore Wind (DC)",
        "offwind-float": "Offshore Wind (Floating)",
        "solar": "Solar",
    }

    # choose pretty names
    pretty_names = {
        "H2 Electrolysis": "H2 electrolysis",
        "H2 pipeline": "H2 pipeline constructed",
        "H2 pipeline retrofitted": "H2 pipeline retrofitted",
        "SMR": "SMR",
        "SMR CC": "SMR CC",
        "Fischer-Tropsch": "Fischer-Tropsch process",
        "H2 Fuel Cell": "H2 fuel cell",
        "H2 for industry": "Industry H2 demand",
        "H2 for shipping": "Shipping H2 demand",
        "OCGT H2 retrofitted": "OCGT H2 retrofitted",
        "Sabatier": "Methanation (Sabatier)",
        "methanation": "Methanation (Sabatier)",
        "land transport fuel cell": "Land transport H2 demand",
        "methanolisation": "Methanol synthesis",
        "offshore wind (AC)": "Offshore Wind (AC)",
        "offshore wind (DC)": "Offshore Wind (DC)",
        "offwind-ac": "Offshore Wind (AC)",
        "offwind-dc": "Offshore Wind (DC)",
        "offshore wind": "Offshore Wnd",
        "onwind": "Onshore Wind",
        "onshore wind": "Onshore Wind",
        "solar PV": "Solar PV (utility)",
        "solar": "Solar PV (utility)",
        "solar rooftop": "Solar PV (rooftop)",
        "hydroelectricity": "Hydroelectricity",
        "uranium": "Uranium",
        "solid biomass": "Solid biomass",
        "solid biomass for industry": "Industry biomass demand",
        "solid biomass for industry CC": "Industry biomass demand CC",
        "gas for industry": "Industry methane demand",
        "gas for industry CC": "Industry methane demand CC",
        "shipping oil": "Shipping oil demand",
        "shipping methanol": "Shipping methanol demand",
        "oil emissions": "Oil emissions",
        "shipping oil emissions": "Shipping oil emissions",
        "shipping methanol emissions": "Shipping methanol emissions",
        "process emissions": "Process emissions",
        "process emissions CC": "Process emissions CC",
        "agriculture machinery oil emissions": "Agriculture machinery oil emissions",
        "gas": "Methane",
        "oil": "Oil",
        "oil primary": "Oil",
        "coal": "Coal",
        "oil boiler": "Oil boiler",
        "gas boiler": "Gas boiler",
        "nuclear": "Nuclear",
        "lignite": "Lignite",
        "land transport oil": "Land transport oil demand",
        "land transport EV": "Land transport EV",
        "naphtha for industry": "Industry naphtha demand",
        "low-temperature heat for industry": "Industry low-temperature heat demand",
        "kerosene for aviation": "Aviation kerosene demand",
        "industry electricity": "Industry electricity demand",
        "heat": "Heat demand",
        "electricity": "Electricity demand",
        "co2": "CO2 emissions",
        "biomass boiler": "Biomass boiler",
        "agriculture machinery oil": "Agriculture machinery oil demand",
        "agriculture heat": "Agriculture heat demand",
        "agriculture electricity": "Agriculture electricity demand",
        "hot water storage": "Hot water storage",
        "resistive heater": "Resistive heater",
        "air heat pump": "Heat pump (air)",
        "ground heat pump": "Heat pump (ground)",
        "electricity distribution grid": "Electricity distribution grid",
        "transmission lines": "Transmission lines",
        "Reservoir & Dam": "Reservoir \& Dam",
        "ror": "Run of River",
        "offwind-float": "Offshore Wind (Floating)",
        "solar-hsat": "Solar (HSAT)",
        "residential rural solar thermal": "Solar Thermal",
        "urban central solar thermal": "Solar Thermal",
        "rural solar thermal": "Solar Thermal",
        "residential urban decentral solar thermal": "Solar Thermal",
        "urban decentral solar thermal": "Solar Thermal",
    }

    # define tech colors and update with custom
    tech_colors = snakemake.config["plotting"]["tech_colors"]
    tech_colors.update(tech_colors_custom)
    tech_colors = {(pretty_gen[k] if k in pretty_gen else k):v  for (k,v) in tech_colors.items()}

    # get threshold capacity
    threshold = snakemake.config["existing_capacities"]["threshold_capacity"]

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    map_opts_params = snakemake.params.plotting["map"]
    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    time_series_params = snakemake.params.plotting["time_series"]
    start_date = time_series_params["start_date"]
    end_date = time_series_params["end_date"]

    n = pypsa.Network(snakemake.input.network)

    plot_curtailment(n, regions, path=snakemake.output.curtailment_map, focus_de=True, show_fig=False)

    plot_line_loading(n, regions, path=snakemake.output.line_loading_map, focus_de=True, value="mean", show_fig=False)

    plot_power_time_series(n, start_date, end_date, path=snakemake.output.energy_balance, focus_component=["Generator","StorageUnit","Link"], focus_de=True, show_fig=False)
    
    plot_power_time_series(n, start_date, end_date, path=snakemake.output.storage_energy_balance, focus_component=["StorageUnit"], focus_de=True, show_fig=False)

    plot_system_cost(n, snakemake.input.nodal_costs, path=snakemake.output.system_cost, focus_component=["Generator","StorageUnit","Link","Store"], focus_de=True, electricity_only=True, show_fig=False)

    plot_system_cost(n, snakemake.input.nodal_costs, path=snakemake.output.storage_system_cost, focus_component=["StorageUnit"], focus_de=True, electricity_only=True, show_fig=False)