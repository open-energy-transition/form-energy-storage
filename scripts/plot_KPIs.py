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
        vmax=60,
        vmin=40,
        legend_kwds={
            "label": "Avg. Electricity Price [EUR/MWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    legend_x = -0.37
    legend_y = 0.57
    sizes = [6, 3, 1]  # GW

    labels = [f"{s} GW" for s in sizes]
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

    # get threshold capacity
    threshold = snakemake.config["existing_capacities"]["threshold_capacity"]

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    map_opts_params = snakemake.params.plotting["map"]
    regions = gpd.read_file(snakemake.input.regions).set_index("name")

    n = pypsa.Network(snakemake.input.network)

    # create save path for figures
    # save_path = f"data_exchange_pypsa/{'_'.join(run_name.split('_')[:2])}/results/{scenario}"
    # os.makedirs(save_path, exist_ok=True)

    plot_curtailment(n, regions, path=snakemake.output.curtailment_map, focus_de=True, show_fig=False)

