# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve cost data from ``technology-data``.
"""

import logging
from pathlib import Path
import pandas as pd
import re

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_cost_data", year=2030)
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    version = snakemake.params.version
    if "/" in version:
        baseurl = f"https://raw.githubusercontent.com/{version}/outputs/"
    else:
        baseurl = f"https://raw.githubusercontent.com/PyPSA/technology-data/{version}/outputs/"
    filepath = Path(snakemake.output[0])
    url = baseurl + filepath.name

    print(url)

    to_fn = Path(rootpath) / filepath

    print(to_fn)

    logger.info(f"Downloading technology data from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, to_fn, disable=disable_progress)

    df_costs = pd.read_csv(to_fn)

    if df_costs.columns[0] == '404: Not Found':
        logger.warning(f"Unable to download technology data from '{url}'.")

        # Possible reason for missing file: the years are atypical
        typical_years = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
        year = int(re.findall(r'\d+', filepath.name)[0])
        if year not in typical_years:
            closest_year = min(typical_years, key=lambda x: abs(x-year))
            logger.info(f"Technology data for '{year}' is not available, using '{closest_year}' instead.")

            url = baseurl + f"costs_{closest_year}.csv"
            logger.info(f"Downloading technology data from '{url}'.")
            progress_retrieve(url, to_fn, disable=disable_progress)

    logger.info(f"Technology data available at at {to_fn}")
