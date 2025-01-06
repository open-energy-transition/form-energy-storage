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
from pypsa.plot import add_legend_circles, add_legend_patches
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from matplotlib import rc
from plot_power_network import load_projection
from matplotlib.colors import Normalize
from plot_summary import preferred_order
from pypsa.statistics import get_bus_and_carrier_and_bus_carrier, get_country_and_carrier
from plot_validation_cross_border_flows import sort_one_country, color_country
import country_converter as coco

cc = coco.CountryConverter()

# activate latex text rendering
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'sans-serif': ['Computer Modern Sans serif']})
rc('axes', **{'edgecolor': 'None', 'titlesize': 18,'titleweight': 'bold', "labelsize": 14})
rc('figure', **{'edgecolor': 'None'})
rc('patch', **{'edgecolor': 'None'})
rc('savefig', **{'edgecolor': 'None'})
rc('legend', **{'fontsize': 12, 'title_fontsize': 12})
rc('xtick', **{'labelsize': 12})
rc('ytick', **{'labelsize': 12})
# plt.style.use(["ggplot"])

def plot_curtailment(network, regions, path, show_fig=True, focus_de=True, legend_circles=[5, 3, 1], bus_size_factor_ = 1e4, vmax_price=105, vmin_price=45):
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
    curtailment_elec = curtailment_elec.loc[regions.index, 0].div(1e3)  # GWh
    electricity_price = n.buses_t.marginal_price.loc[:, regions.index].mean(axis=0)
    regions["elec_price"] = (
        electricity_price
    )

    bus_size_factor = bus_size_factor_

    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    if focus_de:
        curtailment_elec = curtailment_elec.filter(like="DE", axis=0)
        regions = regions.filter(like="DE", axis=0)
        for c in n.iterate_components(n.branch_components):
            c.df.drop(c.df.index[~((c.df.bus0.str.startswith("DE")) | (c.df.bus1.str.startswith("DE")))], inplace=True)
        map_opts["boundaries"] = [4, 17, 46, 56]
        bus_size_factor = bus_size_factor_

    # calculate total curtailment
    total_curtailment = curtailment_elec.groupby(level=0).sum().sum()
    logger.info(f"Total annual curtailment of: {total_curtailment} GWh")

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
        vmax=vmax_price,
        vmin=vmin_price,
        legend_kwds={
            "label": "Avg. Electricity Price [EUR/MWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    legend_x = -0.47
    legend_y = 0.57
    sizes = legend_circles  # TWh

    labels = [f"{s} TWh" for s in sizes]
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

    ax.set_title(f'Total annual curtailment: \n {round(total_curtailment/1e3, 3)} TWh',
                loc='left',
                x=legend_x + 0.03,
                y=legend_y + 0.05,
                fontsize=12
                )

    carriers = curtailment_elec.index.get_level_values(1).unique()
    colors = [tech_colors[c] for c in carriers]
    labels = list(carriers)

    labels = list(pd.Series(labels).replace(pretty_names))

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(legend_x, legend_y),
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

    ax.set_title(f'{title} AC line loading: {line_ave} \%\n{title} DC link loading: {link_ave} \%',
                loc='left',
                x=0.28,
                y=-0.12,
                fontsize=12
                )

    plt.colorbar(collection[1], label="line loading [\%]", ax=ax)

    if show_fig:
        fig.show()
    fig.savefig(path, bbox_inches="tight")

def plot_energy_trade(network, countries, path):
    n = network.copy()

    # Remove EU bus before starting
    to_drop_links=n.links[(n.links.bus0.str[:2] == "EU") | (n.links.bus1.str[:2] == "EU")].index
    n.remove("Link",to_drop_links)

    # Preparing network data to be shaped similar to ENTSOE datastructure
    df_links = n.links_t.p0.rename(
        columns=dict(n.links.bus0.str[:2] + " - " + n.links.bus1.str[:2])
    )
    df_lines = n.lines_t.p0.rename(
        columns=dict(n.lines.bus0.str[:2] + " - " + n.lines.bus1.str[:2])
    )
    df = pd.concat([df_links, df_lines], axis=1)

    # Drop internal country connection
    df.drop(
        [c for c in df.columns if c[:2] == c[5:]], axis=1, inplace=True
    )

    # align columns name
    for c1 in df.columns:
        for c2 in df.columns:
            if c1[:2] == c2[5:] and c2[:2] == c1[5:]:
                df = df.rename(columns={c1: c2})

    df = df.T.groupby(lambda x: x).sum().T

    df_new = pd.DataFrame()

    for country in countries:
        df_country = sort_one_country(country, df)
        df_country_new = pd.DataFrame(data= n.snapshot_weightings.objective @ df_country)
        df_country_new = df_country_new.T.rename(index={"objective":cc.convert(country, to="name_short")})
        df_new = pd.concat([df_country_new, df_new])

    # MWh to GWh
    df_new = df_new/1e3
    
    df_new = df_new.T
    df_new = df_new.rename(index={link:link[5:] for link in df_new.index})
    df_new = df_new.groupby(lambda x: x).sum()
    df_new["color"] = [color_country[country] for country in df_new.index]
    df_new = df_new.rename(index={country:cc.convert(country, to="name_short") for country in df_new.index})

    plot_kw = {"title": "Total Energy Trade", "ylabel": "Import (-)/Export (+) [GWh]"}

    plot_by_country(df_new, plot_kw, path)

def prepare_energy_balance(n):
    df = n.statistics.energy_balance(groupby=get_bus_and_carrier_and_bus_carrier, aggregate_time=False)
    df = df.rename(index=n.buses.country, level="bus")

    return df

def filter_plot_energy_balance(network, dataframe, kpi_param, path):

    carrier_filter = kpi_param.get("carrier_filter", "electricity")
    group_carrier = kpi_param.get("group_carrier", None)
    plot_kw = kpi_param.get("plot_kw", {})
    
    n = network.copy()
    df = dataframe.copy()

    # Model 1: Electricity system with generators and batteries out
    if carrier_filter == "electricity+":
        extract_carrier = ["solar rooftop",'BEV charger', 'V2G']
        df = df.query('index.get_level_values("bus_carrier").isin(["AC","DC"]) or index.get_level_values("carrier").isin(@extract_carrier)', engine="python")
        df = df.rename(index={i:'Electricity trade' for i in ["AC","DC"]}, level="carrier")
        df_low = df[df.index.get_level_values("carrier").isin(["solar rooftop",'BEV charger', 'V2G'])]
        df_low = df_low[df_low.index.get_level_values("bus_carrier").isin(["low voltage"])] * -1
        df_low = df_low.rename(index={i:'electricity distribution grid' for i in extract_carrier}, level="carrier")
        
        df = pd.concat([df,df_low])
        df = df.groupby(["bus","carrier"]).sum()
        df.loc[df.index.get_level_values("carrier").isin(['electricity distribution grid']),:] *= -1
    
        line_carrier = "electricity distribution grid"
    
    # Model 2: Just electricty system
    if carrier_filter == "electricity":
        df = df[df.index.get_level_values("bus_carrier").isin(["AC","DC"])]
        df = df.rename(index={i:'Electricity trade' for i in ["AC","DC"]}, level="carrier")
        df = df.groupby(["bus","carrier"]).sum()
        df.loc[df.index.get_level_values("carrier").isin(['electricity distribution grid']),:] *= -1
    
        line_carrier = "electricity distribution grid"
    
    # Model 3: Just low voltage system
    if carrier_filter == "low voltage":
        df = df[df.index.get_level_values("bus_carrier").isin(["low voltage"])]
        df = df.groupby(["bus","carrier"]).sum()
        df.loc[df.index.get_level_values("carrier").isin(['electricity distribution grid']),:] *= -1
    
        line_carrier = "electricity distribution grid"
    
    # Model 4: Just the heating system
    if carrier_filter == "heat":
        extract_carrier = ['residential rural heat', 'urban central heat','residential urban decentral heat']
        df = df[df.index.get_level_values("bus_carrier").isin(extract_carrier)]
        df = df.rename(index={i:'heat' for i in extract_carrier}, level="carrier")
        df = df.groupby(["bus","carrier"]).sum()
        df.loc[df.index.get_level_values("carrier").isin(['heat']),:] *= -1
        
        line_carrier = "heat"

    # Model 5: Hydrogen input and output
    if carrier_filter == "hydrogen":
        df = df[df.index.get_level_values("bus_carrier").isin(['Hydrogen Storage'])]
        df = df.groupby(["bus","carrier"]).sum()
        df = df[~df.index.get_level_values("carrier").isin(['H2 Store'])]

        line_carrier = 'H2 Store'
    
    df = df.groupby(["carrier"]).sum()
    
    #Set color
    colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "green")
    if group_carrier == "pretty":
        line_carrier = pretty_names[line_carrier] if line_carrier in pretty_names else line_carrier
        df = df.rename(index=pretty_names).groupby(["carrier"]).sum()
        colors = colors.rename(index=pretty_names)
        colors = colors[~colors.index.duplicated(keep='first')]
        
    elif group_carrier == "sector": 
        line_carrier = sector_colors[line_carrier] if line_carrier in sector_colors else line_carrier
        df = df.rename(index=sector_names).groupby(["carrier"]).sum()
        colors = sector_colors
        
    df = df.T

    colors['Electricity trade'] = 'lightgrey'
    
    for c in df.columns:
        if c in line_carrier:
            continue
        
        if (df[c] < 0).values.any():
            df[c + " Discharge"] = df[c].clip(lower=0)
            colors[c + " Discharge"] = colors[c]
            df[c + " Charge"] = df[c].clip(upper=0)
            colors[c + " Charge"] = colors[c]
            df = df.drop(columns=c)

    #kW to MW
    df = df/1e3
    
    #remove if smaller than 1 MWh
    to_drop = abs(n.snapshot_weightings.objective @ df) < 1
    df = df.loc[:,~to_drop]
    
    #cut into the specific timesteps
    start_date = kpi_param.get("start_date", df.index[0])
    end_date = kpi_param.get("end_date", df.index[-1])
    df = df[pd.Timestamp(start_date):pd.Timestamp(end_date)]
    
    fig, ax = plt.subplots(figsize=(12,5))
    
    df_plot = df[df.columns.difference([line_carrier])]
    
    order = ((df_plot.diff().abs().sum() / df_plot.sum()).sort_values().index)
    df_plot = df_plot[order]
    df_plot.plot.area(ax = ax, legend=False, color = [colors[c] for c in df_plot.columns], **plot_kw)
    
    df_demand = df.columns.intersection([line_carrier])
    if df_demand.any():
        df[df_demand].plot(ax = ax, color = "black", linestyle='dashed', **plot_kw)
    
    ax.grid(axis="y")
    ax.set_xlabel("")
    ax.set_facecolor("white")
    
    df_legend = pd.DataFrame()
    df_legend["handle"], df_legend["label"] = ax.get_legend_handles_labels()
    df_legend["label"] = [w.replace(' Charge', '').replace(' Discharge', '') for w in df_legend["label"]]
    #df_legend["label"] = [w.replace('Electricity distribution grid', 'Electricity demand') for w in df_legend["label"]]
    df_legend = df_legend.drop_duplicates(subset = ["label"],keep = 'first').iloc[::-1]
    
    ax.legend(df_legend["handle"], df_legend["label"], loc = "upper center", bbox_to_anchor = (0.5, -0.15), frameon = False, ncol = 3, 
                title = "Components", title_fontproperties = {'weight':'bold'})

    fig.savefig(path, bbox_inches="tight")

def filter_plot_SOC(network, dataframe, kpi_param, path):

    carrier_filter = kpi_param.get("carrier_filter", "Li-Ion Battery Storage")
    group_carrier = kpi_param.get("group_carrier", None)
    plot = kpi_param.get("plot", None)
    plot_kw = kpi_param.get("plot_kw", {})
    
    n = network.copy()
    df = dataframe.copy()

    if isinstance(carrier_filter, list):
        df = df[df.index.get_level_values("carrier").isin(carrier_filter)]
    elif isinstance(carrier_filter, str):
        df = df[df.index.get_level_values("carrier").isin([carrier_filter])]

    # readjust charge and discharge efficiency in storage units
    n_storage_units = n.storage_units
    n_storage_units["carrier"] = n_storage_units["carrier"].map(n.carriers.nice_name)
    efficiency_store = n_storage_units.groupby(["carrier"]).efficiency_store.mean()
    efficiency_dispatch = n_storage_units.groupby(["carrier"]).efficiency_dispatch.mean()
    roundtrip_efficiency = efficiency_store * efficiency_dispatch

    df_charge = df.copy(deep=True).clip(upper=0)
    df_discharge = df.copy(deep=True).clip(lower=0)

    for c in efficiency_store.index:
        df_charge[df_charge.index.get_level_values("carrier").isin([c])] *= roundtrip_efficiency[c]

    df = df_charge + df_discharge

    df = df.groupby(["carrier"]).sum()

    # Factor in snapshot weightings
    df = df.mul(n.snapshot_weightings.stores, axis=1)

    # Accumulate energy based on the value of the previous energy level
    for i in range(len(df.columns)):
        if i == 0:
            continue
        df.iloc[:,i] = df.iloc[:,i] + df.iloc[:,i-1]

    #Set color
    colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "green")
    if group_carrier == "pretty":
        df = df.rename(index=pretty_names).groupby(["carrier"]).sum()
        colors = colors.rename(index=pretty_names)
        colors = colors[~colors.index.duplicated(keep='first')]

    df = df.T

    #MWh to GWh
    df = df/1e3

    #remove if smaller than 1 MWh
    to_drop = abs(df.max()) < 1
    df = df.loc[:,~to_drop]

    # Set the minimum value to zero
    df = (df - df.min())

    if plot == "share":
        # Set the maximum value to 100
        df = (df/df.max() * 100)

        if not plot_kw.get("ylim", False):
            plot_kw["ylim"] = [0, 100]
    else:
        if not plot_kw.get("ylim", False):
            plot_kw["ylim"] = [0, None]

    #cut into the specific timesteps
    start_date = kpi_param.get("start_date", df.index[0])
    end_date = kpi_param.get("end_date", df.index[-1])
    df = df[pd.Timestamp(start_date):pd.Timestamp(end_date)]

    fig, ax = plt.subplots(figsize=(12,5))

    df.plot(ax = ax, color = [colors[c] for c in df.columns], **plot_kw)

    ax.grid(axis="y")
    ax.set_xlabel("")
    ax.set_facecolor("white")
    
    df_legend = pd.DataFrame()
    df_legend["handle"], df_legend["label"] = ax.get_legend_handles_labels()
    ax.legend(df_legend["handle"], df_legend["label"], loc = "upper center", bbox_to_anchor = (0.5, -0.15), frameon = False, ncol = 3, 
                title = "Components", title_fontproperties = {'weight':'bold'})

    fig.savefig(path, bbox_inches="tight")

def calculate_emission(network, countries):
    n = network.copy()
    
    index_with_emission = (n.links == 'co2 atmosphere').T.any()
    df_links = n.links
    
    for num in [0,1,2,3,4]:
        index_with_co2 = (df_links[f"bus{num}"] == 'co2 atmosphere')
        df_links.loc[index_with_co2,"co2_num"] = str(num)
    
        df_links[f"p{num}"] =  n.snapshot_weightings.stores @ n.links_t[f"p{num}"]
        df_links[f"country{num}"] = df_links[f"bus{num}"].map(n.buses.country)
    
    df_links["p_co2"] = [-df_links.loc[i,f"p{df_links.loc[i,"co2_num"]}"] if not pd.isna(df_links.loc[i,"co2_num"]) else np.nan for i in df_links.index]
    #[df_links[f"p{num}"] if not pd.isna(num) else np.nan for num in df_links["co2_num"]]
    
    for country in n.buses.country.unique():
        if country == "":
            continue
        index_with_country = (df_links[[f"country{num}" for num in [0,1,2,3,4]]] == country).T.any()
        df_links.loc[index_with_country,"country"] = country
    
    df_links.loc[~df_links.country.isin(countries),"country"] = "EU"
    
    df = pd.DataFrame(df_links.groupby(["country","carrier"])["p_co2"].sum()).unstack("country")
    df.columns = df.columns.get_level_values("country")
    df = df.replace(0, np.nan)
    df = df.dropna(axis = 0, how = "all")

    # convert tCO2 to MtCO2
    df = df/1e6

    return df

def calculate_system_cost_csv(network, csv_path, countries):
    n = network.copy()
    
    df = pd.read_csv(csv_path, index_col=list(range(3)), header=list(range(4)))
    df = df.set_axis(['carrier', 'cost'], axis=1)
    df.index = df.index.rename(["country"], level=[2])
    df = df.rename(index=n.buses.country, level=2)
    df = df.reset_index(["country"])
    df.loc[~df.country.isin(countries),"country"] = "EU"
    
    df = pd.DataFrame(df.groupby(["country","carrier"])["cost"].sum()).unstack("country")
    df.columns = df.columns.get_level_values("country")
    
    # convert to billions
    df = df / 1e9

    return df

def calculate_capacity_csv(network, csv_path, countries):
    n = network.copy()
    
    df = pd.read_csv(csv_path, index_col=list(range(2)), header=list(range(3)))
    df = df.set_axis(['carrier', 'capacity'], axis=1)
    df.index = df.index.rename(["component","country"], level=[0,1])
    df = df.rename(index=n.buses.country, level=1)
    df = df.reset_index(["country","component"])
    df.loc[~df.country.isin(countries),"country"] = "EU"
    
    # convert links power gen techs capacity from MWth to MWel by multiplying by efficiency
    df = df.set_index(["carrier"])
    links_i = df.query("component == 'links' and carrier in @power_generation_tech").index
    power_ge_links_eff = n.links.query("carrier in @links_i").groupby("carrier").first().efficiency
    df["efficiency"] = df.index.map(power_ge_links_eff)
    df.loc[links_i, "capacity"] *= df.loc[links_i].efficiency

    df = pd.DataFrame(df.groupby(["country","carrier"])["capacity"].sum()).unstack("country")
    df.columns = df.columns.get_level_values("country")

    # convert from MW to GW
    df = df / 1e3
    
    return df

def calculate_capacity(network, countries, stats = "optimal", storage = False):
    
    n = network.copy()
    
    if stats == "install":
        df = pd.DataFrame(n.statistics.installed_capacity(groupby=get_country_and_carrier, storage = storage))
    elif stats == "optimal":
        df = pd.DataFrame(n.statistics.optimal_capacity(groupby=get_country_and_carrier, storage = storage))
    elif stats == "expand" and storage == True:
        df = pd.DataFrame(n.statistics.optimal_capacity(groupby=get_country_and_carrier, storage = True)).subtract(pd.DataFrame(n.statistics.installed_capacity(groupby=get_country_and_carrier, storage = True)), fill_value=0)
    elif stats == "expand" and storage == False:
        df = pd.DataFrame(n.statistics.expanded_capacity(groupby=get_country_and_carrier))
    df = df.reset_index(["country","component"])
        
    # convert links power gen techs capacity from MWth to MWel by multiplying by efficiency
    links_i = df.query("component == 'links' and carrier in @power_generation_tech").index
    power_ge_links_eff = n.links.query("carrier in @links_i").groupby("carrier").first().efficiency
    df["efficiency"] = df.index.map(power_ge_links_eff)
    df.loc[links_i, 0] *= df.loc[links_i].efficiency
    
    df.loc[~df.country.isin(countries),"country"] = "EU"
    
    df = pd.DataFrame(df.groupby(["country","carrier"])[0].sum()).unstack("country")
    df.columns = df.columns.get_level_values("country")

    # convert from MWh to GWh
    df = df / 1e3

    return df

def calculate_generation(network, countries):
    n = network.copy()
    
    df = pd.DataFrame(n.statistics.energy_balance(groupby=get_country_and_carrier))
    df = df.reset_index(["country","component"])
    df.loc[~df.country.isin(countries),"country"] = "EU"

    df = pd.DataFrame(df.groupby(["country","carrier"])[0].sum()).unstack("country")
    df.columns = df.columns.get_level_values("country")

    # convert from MW to GW
    df = df / 1e3
    
    return df

def filter_and_rename(network, df, carrier_filter = None, group_carrier = None):
    n = network.copy()

    # Replace countries with short country names
    short_name = cc.convert(names = df.columns, src = 'ISO2', to = 'name_short', not_found=None)
    if len(df.columns) == 1:
        df.columns = [short_name]
    else:
        df.columns = short_name

    # convert to nice name
    df["nice_name"] = list(pd.Series(df.index).replace(n.carriers.nice_name))
    df = df.set_index("nice_name")

    # Create your own filter
    stats = pd.DataFrame(n.statistics.energy_balance(groupby=get_bus_and_carrier_and_bus_carrier))

    if carrier_filter == "electricity":
        carrier = stats[stats.index.get_level_values("bus_carrier").isin(["AC"])].index.get_level_values("carrier").unique()
        df = df.loc[df.index.isin(carrier),:]

    elif carrier_filter == "electricity+":
        carrier = stats[stats.index.get_level_values("bus_carrier").isin(["AC"])].index.get_level_values("carrier").unique()
        carrier = carrier.append(stats.filter(like="water tanks",axis=0).index.get_level_values("carrier").unique())
        carrier = carrier.append(stats.filter(like="EV battery",axis=0).index.get_level_values("carrier").unique())
        df = df.loc[df.index.isin(carrier),:]

    elif carrier_filter == "storage":
        carrier = pd.Index([
            "Iron-Air Battery Storage","Li-Ion Battery Storage","lair","vanadium","pair",
            "H2 Fuel Cell","H2 Store","H2 Electrolysis","Pumped Hydro Storage",
            "BEV charger","EV battery","V2G",'residential rural water tanks',
            'urban central water tanks','residential urban decentral water tanks',
        ])
        df = df.loc[df.index.isin(carrier),:]

    elif carrier_filter == "storage-cap":
        carrier = pd.Index([
            "Iron-Air Battery Storage","Li-Ion Battery Storage","lair","vanadium","pair",
            "Adiabatic CAES","H2 Fuel Cell","H2 Electrolysis","Pumped Hydro Storage",
            "V2G",'residential rural water tanks',
            'urban central water tanks','residential urban decentral water tanks',
        ])
        df = df.loc[df.index.isin(carrier),:]

    elif carrier_filter == "storage-energy":
        carrier = pd.Index([
            "Iron-Air Battery Storage","Li-Ion Battery Storage","lair","vanadium","pair",
            "Adiabatic CAES","H2 Store","Pumped Hydro Storage","EV battery",'residential rural water tanks',
            'urban central water tanks','residential urban decentral water tanks',
        ])
        df = df.loc[df.index.isin(carrier),:]

    elif carrier_filter == "power":
        carrier = pd.Index([
            "solar rooftop","Solar","solar-hsat","Onshore Wind","Offshore Wind (DC)",
            "Offshore Wind (AC)","Offshore Wind (Floating)","Run of River","Reservoir & Dam",
            "Open-Cycle Gas","Combined-Cycle Gas","nuclear","oil","lignite","coal"
        ])
        df = df.loc[df.index.isin(carrier),:]

    # Group carrier
    if group_carrier == "pretty":
        df = df.rename(index=pretty_names).groupby(["nice_name"]).sum()
    elif group_carrier == "sector": 
        df = df.rename(index=sector_names).groupby(["nice_name"]).sum()

    to_drop = df.loc[abs(df.T.sum()) < 1,:].index
    # An exception if Iron Air Battery Storage is still in
    if 'Iron-Air Battery Storage' in to_drop:
        to_drop = to_drop.drop('Iron-Air Battery Storage')
    df = df.drop(to_drop)
    
    # Sort the index by preferred order and set color
    colors = n.carriers.set_index("nice_name").color.where(lambda s: s != "", "green")
    
    if group_carrier == "pretty":
        new_index = preferred_order_pretty_names_reversed.intersection(df.index).append(
            df.index.difference(preferred_order_pretty_names)
        )
        colors = colors.rename(index=pretty_names)
        colors = colors[~colors.index.duplicated(keep='first')]
        
    elif group_carrier == "sector": 
        new_index = preferred_order_sector_names_reversed.intersection(df.index).append(
            df.index.difference(preferred_order_sector_names)
        )
        colors = sector_colors
    else:
        new_index = preferred_order_reversed.intersection(df.index).append(
            df.index.difference(preferred_order)
        )

    df["color"] = [colors[i] for i in df.index]

    return df.loc[new_index,:]
        
def plot_by_country(df, plot_kw, path):

    fig, ax = plt.subplots(figsize=(12, 7))
    
    df.drop(columns={"color"}).T.plot(
        kind="bar",
        width = 0.9,
        ax=ax,
        stacked=True,
        color=df["color"],
        rot=90,
        legend=False,
        **plot_kw,
    )

    ax.grid(axis="x")
    ax.set_xlabel("")
    ax.set_facecolor("white")

    if not plot_kw.get("ylim", False):
        ymin, ymax = ax.get_ylim()
        margin = abs(ymax - ymin) * 0.1
        if ymin == 0:
            ymin = margin
        ax.set_ylim(ymin - margin, ymax + margin)

    handles, labels = ax.get_legend_handles_labels()
    
    handles.reverse()
    labels.reverse()
    labels = [label.replace(' &', '').replace(' and', '') for label in labels] #NOTE: Latex hates '&' strings because its their seperator
    
    ax.legend(handles, labels, loc = "upper center", bbox_to_anchor = (0.5, -0.20), frameon = False, ncol = 3, 
                title = "Components", title_fontproperties = {'weight':'bold'})

    fig.savefig(path, bbox_inches="tight")

def plot_in_detail(df, plot_kw, path):
    
    c = df.columns.difference(["color"])
    df = df[["color"]].assign(Total=df[c].sum(axis=1))

    fig, ax = plt.subplots(figsize=(6, 8))

    df.drop(columns={"color"}).T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=df["color"],
        **plot_kw,
    )

    if not plot_kw.get("ylim", False):
        ymin, ymax = ax.get_ylim()
        margin = abs(ymax - ymin) * 0.1
        if ymin == 0:
            ymin = margin
        ax.set_ylim(ymin - margin, ymax + margin)
    
    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    ax.grid(axis="x")

    labels = [label + ": \n " + str(round(df.loc[label,"Total"],2)) for label in labels]
    labels = [label.replace(' &', '').replace(' and', '') for label in labels] #NOTE: Latex hates '&' strings because its their seperator

    ax.legend(
        handles, labels, ncol=1, loc="center left", bbox_to_anchor=[1, 0.5], frameon=False, 
        title = f"Total: \n{round(df.drop(columns={"color"}).sum().iloc[0],2)}", alignment="left"
    )

    ax.set_facecolor("white")
    fig.tight_layout()
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

    power_generation_tech = pd.Index([
        "solar rooftop","Solar","solar-hsat","Onshore Wind","Offshore Wind (DC)",
        "Offshore Wind (AC)","Offshore Wind (Floating)","Run of River","Reservoir & Dam",
        "Open-Cycle Gas","Combined-Cycle Gas","nuclear","oil","lignite","coal"
    ])

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
        "lair": "Liquid Air Battery Storage",
        "pair": "Compressed Air Battery Storage",
        "vanadium": "Vanadium-Redox Battery Storage",
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
        "Reservoir & Dam": "Reservoir & Dam", # Note, revert this back to /&
        "ror": "Run of River",
        "offwind-float": "Offshore Wind (Floating)",
        "solar-hsat": "Solar (HSAT)",
        "residential rural solar thermal": "Solar Thermal",
        "urban central solar thermal": "Solar Thermal",
        "rural solar thermal": "Solar Thermal",
        "residential urban decentral solar thermal": "Solar Thermal",
        "urban decentral solar thermal": "Solar Thermal",
        "urban central water tanks": "Water Tanks",
        "residential rural water tanks": "Water Tanks",
        "residential urban decentral water tanks": "Water Tanks",
        "V2G":"Vehicle-to-Grid",
    }

    sector_colors = {
        "Power Sector": '#ffbf2b',
        "Transport Sector": '#baf238',
        "Heating Sector": '#db6a25',
        "CCUS": '#d2d2d2',
        "Primary Fuel": '#692e0a'
    }

    sector_names = {
        "Iron-Air Battery Storage": "Power Sector",
        "Li-Ion Battery Storage": "Power Sector",
        "lair": "Power Sector",
        "vanadium": "Power Sector",
        "pair": "Power Sector",
        "H2 Fuel Cell": "Power Sector",
        "H2 Store": "Power Sector",
        "H2 Electrolysis": "Power Sector",
        "Pumped Hydro Storage": "Power Sector",
        "solar rooftop": "Power Sector",
        "Solar": "Power Sector",
        "solar-hsat": "Power Sector",
        "Onshore Wind": "Power Sector",
        "Offshore Wind (DC)": "Power Sector",
        "Offshore Wind (AC)": "Power Sector",
        "Offshore Wind (Floating)": "Power Sector",
        "Run of River": "Power Sector",
        "Reservoir & Dam": "Power Sector",
        "Open-Cycle Gas": "Power Sector",
        "Combined-Cycle Gas": "Power Sector",
        "nuclear": "Power Sector",
        "AC": "Power Sector",
        "DC": "Power Sector",
        "electricity distribution grid": "Power Sector",
        "BEV charger": "Transport Sector",
        "EV battery": "Transport Sector",
        "V2G": "Transport Sector",
        "land transport oil": "Transport Sector",
        "residential rural gas boiler": "Heating Sector",
        'residential rural ground heat pump': "Heating Sector",
        "residential rural resistive heater": "Heating Sector",
        "residential rural solar thermal": "Heating Sector",
        "residential rural water tanks": "Heating Sector",
        "residential rural water tanks discharger": "Heating Sector",
        "residential rural water tanks discharger": "Heating Sector",
        "residential urban decentral air heat pump": "Heating Sector",
        'residential urban decentral gas boiler': "Heating Sector",
        "residential urban decentral resistive heater": "Heating Sector",
        "residential urban decentral solar thermal": "Heating Sector",
        "residential urban decentral water tanks": "Heating Sector",
        "residential urban decentral water tanks charger": "Heating Sector",
        "rural oil boiler": "Heating Sector",
        "rural gas boiler": "Heating Sector",
        "rural ground heat pump": "Heating Sector",
        "services rural gas boiler": "Heating Sector",
        "services rural air heat pump": "Heating Sector",
        "services rural ground heat pump": "Heating Sector",
        "services rural resistive heater": "Heating Sector",
        "services rural water tanks charger": "Heating Sector",
        "services urban decentral gas boiler": "Heating Sector",
        "services urban decentral water tanks discharger": "Heating Sector",
        "services urban decentral resistive heater": "Heating Sector",
        "urban central CHP": "Heating Sector",
        "urban central CHP CC": "Heating Sector",
        "urban central air heat pump": "Heating Sector",
        "urban central gas boiler": "Heating Sector",
        "urban central resistive heater": "Heating Sector",
        "urban central solar thermal": "Heating Sector",
        "urban central solid biomass CHP": "Heating Sector",
        "urban central water tanks": "Heating Sector",
        "urban central water tanks charger": "Heating Sector",
        "urban central water tanks discharger": "Heating Sector",
        "urban decentral oil boiler": "Heating Sector",
        "urban decentral air heat pump": "Heating Sector",
        "urban decentral gas boiler": "Heating Sector",
        "co2": "Primary Fuel",
        "co2 stored": "CCUS",
        "co2 sequestered": "CCUS",
        "DAC": "CCUS",
        "SMR CC": "Primary Fuel",
        "SMR": "Primary Fuel",
        "Sabatier": "Primary Fuel",
        "oil refining": "Primary Fuel",
        "oil primary": "Primary Fuel",
        "oil": "Primary Fuel",
        "methanol": "Primary Fuel",
        "gas": "Primary Fuel",
        "lignite": "Primary Fuel",
        "coal": "Primary Fuel",
        "uranium": "Primary Fuel",
    }

    preferred_order = pd.Index([
        "Iron-Air Battery Storage",
        "Li-Ion Battery Storage",
        "lair",
        "vanadium",
        "pair",
        "H2 Fuel Cell",
        "H2 Store",
        "H2 Electrolysis",
        "Pumped Hydro Storage",
        "solar rooftop",
        "Solar",
        "solar-hsat",
        "Onshore Wind",
        "Offshore Wind (DC)",
        "Offshore Wind (AC)",
        "Offshore Wind (Floating)",
        "Run of River",
        "Reservoir & Dam",
        "Open-Cycle Gas",
        "Combined-Cycle Gas",
        "nuclear",
        "AC",
        "DC",
        "electricity distribution grid",
        "BEV charger",
        "EV battery",
        "V2G",
        "land transport oil",
        "residential rural gas boiler",
        "residential rural resistive heater",
        "residential rural solar thermal",
        "residential rural water tanks",
        "residential rural water tanks discharger",
        "residential urban decentral air heat pump",
        "residential urban decentral resistive heater",
        "residential urban decentral solar thermal",
        "residential urban decentral water tanks",
        "residential urban decentral water tanks charger",
        "rural oil boiler",
        "services rural gas boiler",
        "services rural air heat pump",
        "services rural ground heat pump",
        "services rural resistive heater",
        "services rural water tanks charger",
        "services urban decentral gas boiler",
        "services urban decentral water tanks discharger",
        "urban central CHP",
        "urban central CHP CC",
        "urban central air heat pump",
        "urban central gas boiler",
        "urban central resistive heater",
        "urban central solar thermal",
        "urban central solid biomass CHP",
        "urban central water tanks",
        "urban central water tanks charger",
        "urban central water tanks discharger",
        "urban decentral oil boiler",
        "co2",
        "co2 stored",
        "co2 sequestered",
        "DAC",
        "SMR CC",
        "SMR",
        "Sabatier",
        "oil refining",
        "oil primary",
        "oil",
        "methanol",
        "gas",
        "lignite",
        "coal",
        "uranium"
    ])

    # defined preferred order under all naming schemes
    preferred_order_reversed = preferred_order[::-1]
    preferred_order_pretty_names = pd.Index(pd.Series(preferred_order).apply(lambda x: pretty_names[x] if x in pretty_names.keys() else x).drop_duplicates(keep = 'first'))
    preferred_order_pretty_names_reversed = preferred_order_pretty_names[::-1]

    preferred_order_sector_names = pd.Index(pd.Series(preferred_order).apply(lambda x: sector_names[x] if x in sector_names.keys() else x).drop_duplicates(keep = 'first'))
    preferred_order_sector_names_reversed = preferred_order_sector_names[::-1]

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
    # calculate and print emissions
    co2_emissions = n.stores_t.e.filter(like="co2 atmosphere", axis=1).iloc[-1].div(1e6)[0]  # in MtCO2
    logger.info(f"Total annual CO2 emissions of {co2_emissions} MtCO2.")

    plot_curtailment(n, regions, path=snakemake.output.curtailment_map, focus_de=True,
                     show_fig=False, legend_circles=[12, 6, 3], bus_size_factor_=8e4, vmax_price=105, vmin_price=45)
    plot_curtailment(n, regions, path=snakemake.output.curtailment_map_EU, focus_de=False,
                     show_fig=False, legend_circles=[20, 10, 5], bus_size_factor_=4e4, vmax_price=105, vmin_price=45)

    plot_line_loading(n, regions, path=snakemake.output.line_loading_map, focus_de=True, value="mean", show_fig=False)

    # extract all the nessesary statistics
    countries = snakemake.params.countries

    plot_energy_trade(n, countries, path=snakemake.output.energy_trade)

    df_cost = calculate_system_cost_csv(n, snakemake.input.nodal_costs, countries)
    df_capacity_csv = calculate_capacity_csv(n, snakemake.input.nodal_capacity, countries)
    df_gen = calculate_generation(n, countries)
    df_co2 = calculate_emission(n, countries)
    df_eql = prepare_energy_balance(n)

    config_kpi = snakemake.params.kpi
    if config_kpi:
        for fn, kpi_param in config_kpi.items():
            try:
                logger.info(f"Preparing plot {fn}")

                extract_param = kpi_param.get("extract",None)
                include = kpi_param.get("include",False)
                exclude = kpi_param.get("exclude",False)

                if extract_param == "system cost":
                    logger.info("extracting system cost")
                    df = df_cost.copy(deep=True)
                elif extract_param == "capacity":
                    logger.info("extracting capacity")
                    df = df_capacity_csv.copy(deep=True)
                elif extract_param == "capacity stats":
                    stats = kpi_param.get("stats","optimal")
                    storage = kpi_param.get("storage",None)
                    df = calculate_capacity(n, countries, stats = stats, storage = storage)
                elif extract_param == "generation":
                    logger.info("extracting generation")
                    df = df_gen.copy(deep=True)
                elif extract_param == "emission":
                    logger.info("extracting emission")
                    df = df_co2.copy(deep=True)

                #time series plots have their own route
                elif extract_param == "energy balance" or extract_param == "SOC":
                    df = df_eql.copy(deep=True)

                    if include:
                        print(f"include only {include}")
                        df = df[df.index.get_level_values("bus").isin(include)]
                    if exclude:
                        print(f"exclude {exclude}")
                        df = df[~df.index.get_level_values("bus").isin(exclude)]

                    if extract_param == "energy balance":
                        filter_plot_energy_balance(n, df, kpi_param, snakemake.output[fn])
                    elif extract_param == "SOC":
                        filter_plot_SOC(n, df, kpi_param, snakemake.output[fn])
                    continue

                if include:
                    logger.info(f"include only {include}")
                    df = df[include]
                if exclude:
                    logger.info(f"exclude {exclude}")
                    df = df.loc[:,df.columns.difference(exclude)]
                
                carrier_filter = kpi_param.get("carrier_filter",None)
                group_carrier = kpi_param.get("group_carrier",None)
                df = filter_and_rename(n, df, carrier_filter = carrier_filter, group_carrier = group_carrier)
                
                plot_param = kpi_param.get("plot",None)
                plot_kw = kpi_param.get("plot_kw",{})
                if plot_param == "detail":
                    plot_in_detail(df, plot_kw, snakemake.output[fn])
                elif plot_param == "overview":
                    plot_by_country(df, plot_kw, snakemake.output[fn])

            except (TypeError):
                logger.warning(f"Unable to plot {fn} because the datasets is empty after applying the filters")

                fig, ax = plt.subplots(figsize=(8,8))
                ax.set_axis_off()
                fig.text(4,4,f"Unable to plot {fn} because the datasets is empty after applying the filters")
                fig.savefig(snakemake.output[fn])
                continue

