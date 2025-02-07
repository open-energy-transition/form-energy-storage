# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
import os
import logging

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

def connection_limit_ntc(n,fn,year="2035"):
    logger.info(f"Adjusting connection limit between EU countries base on TYNDP_NTC.csv for the year {year}")
    df_ntc = pd.read_csv(fn, index_col=["country1","country2"])

    n.lines["country0"] = n.lines.bus0.map(n.buses.country)
    n.lines["country1"] = n.lines.bus1.map(n.buses.country)
    df_line = pd.DataFrame(round(n.lines.groupby(["country0","country1"]).s_nom.sum()))
    
    n.links = n.links
    n.links["country0"] = n.links.bus0.map(n.buses.country)
    n.links["country1"] = n.links.bus1.map(n.buses.country)
    df_link = pd.DataFrame(round(n.links.query("carrier == 'DC'").groupby(["country0","country1"]).p_nom.sum()))
    
    df = pd.concat([df_ntc[year],df_line, df_link],axis=1).loc[df_ntc.index,:]
    df = df.fillna(0)
    
    df["nom_sum"] = df["s_nom"] + df["p_nom"]
    df["s_nom_share"] = df["s_nom"] * df[year] / df["nom_sum"]
    df["p_nom_share"] = df["p_nom"] * df[year] / df["nom_sum"]
    
    # adjusting the lines
    for country0, country1 in df.index:
        original = df.loc[(country0,country1),"s_nom"]
        value = df.loc[(country0,country1),"s_nom_share"]
        if value == 0:
            continue
        logger.info(f"Set AC transmission limit between {country0} - {country1} from {original} MW to NTC of {value} MW")
        lines = n.lines.query("country0 == @country0 & country1 == @country1")
        proportion = (n.lines.loc[lines.index,"s_nom"]/n.lines.loc[lines.index,"s_nom"].sum())
        n.lines.loc[lines.index,"s_nom_max"] = value * proportion / n.lines.loc[lines.index,"s_max_pu"]
        n.lines.loc[lines.index,["s_nom_min","s_nom"]] = n.lines.loc[lines.index,["s_nom","s_nom_max"]].T.min()
    
    # adjusting the links
    for country0, country1 in df.index:
        original = df.loc[(country0,country1),"p_nom"]
        value = df.loc[(country0,country1),"p_nom_share"]
        if value == 0:
            continue
        logger.info(f"Set DC transmission limit between {country0} - {country1} from {original} to {value}")
        links = n.links.query("country0 == @country0 & country1 == @country1 & carrier == 'DC'")
        proportion = (n.links.loc[links.index,"p_nom"]/n.links.loc[links.index,"p_nom"].sum())
        n.links.loc[links.index,"p_nom_max"] = value * proportion / n.links.loc[links.index,"p_max_pu"]
        n.links.loc[links.index,["p_nom_min","p_nom"]] = n.links.loc[links.index,["p_nom","p_nom_max"]].T.min()
    
        links_rev = n.links.query("country0 == @country1 & country1 == @country0 & carrier == 'DC'")
        proportion = (n.links.loc[links_rev.index,"p_nom"]/n.links.loc[links_rev.index,"p_nom"].sum())
        n.links.loc[links_rev.index,"p_nom_max"] = value * proportion / n.links.loc[links_rev.index,"p_max_pu"]
        n.links.loc[links_rev.index,["p_nom_min","p_nom"]] = n.links.loc[links_rev.index,["p_nom","p_nom_max"]].T.min()
    
    return n

def disable_links_expansion(n, country, carrier):
    logger.info(f"Disable links expansion for {carrier} in {country}")
    n.links.loc[(n.links.bus1.map(n.buses.country).isin(country)) & (n.links.carrier.isin(carrier)),"p_nom_extendable"] = False

    return n

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "final_adjustment_perfect",
            configfiles="../config/test/config.perfect.yaml",
            opts="",
            clusters="5",
            ll="v1.0",
            sector_opts="",
            # planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    n = pypsa.Network(snakemake.input.network)

    n = connection_limit_ntc(n,snakemake.input.ntc,year="2035")

    if not snakemake.params.chp_extendable_DE:
        n = disable_links_expansion(n, ["DE"], ["urban central CHP"])

    n.export_to_netcdf(snakemake.output.network)