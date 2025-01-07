#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config
from entsoe import EntsoePandasClient
from entsoe.exceptions import NoMatchingDataError

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_electricity_loads")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    api_key = snakemake.config["private"]["keys"]["entsoe_api"]
    client = EntsoePandasClient(api_key=api_key)

    start = pd.Timestamp(snakemake.params.snapshots["start"], tz="Europe/Brussels")
    end = pd.Timestamp(snakemake.params.snapshots["end"], tz="Europe/Brussels")

    countries = snakemake.params.countries

    loads = []
    unavailable_countries = []

    for country in countries:
        country_code = country

        try:
            gen = client.query_load(country, start=start, end=end)
            gen = gen.tz_localize(None).resample("1h").mean()
            gen = gen.loc[start.tz_localize(None) : end.tz_localize(None)]
            loads.append(gen)
        except NoMatchingDataError:
            unavailable_countries.append(country)

    if unavailable_countries:
        logger.warning(
            f"Historical electricity loads for countries {', '.join(unavailable_countries)} not available."
        )

    keys = [c for c in countries if c not in unavailable_countries]
    loads = pd.concat(loads, keys=keys, axis=1)
    loads.to_csv(snakemake.output[0])