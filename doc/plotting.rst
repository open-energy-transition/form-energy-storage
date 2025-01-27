..
  SPDX-FileCopyrightText: 2019-2024 The PyPSA-Eur Authors

  SPDX-License-Identifier: CC-BY-4.0

##########################################
Plotting and Summaries
##########################################


Rule ``make_summary``
========================

.. automodule:: make_summary

.. _summary_plot:

Create summary CSV files for all scenario runs including costs, capacities,
capacity factors, curtailment, energy balances, prices and other metrics.

Rule ``plot_summary``
========================

.. automodule:: plot_summary

.. _map_plot:

Creates plots for optimised power network topologies and regional generation,
storage and conversion capacities built.

Rule ``plot_power_network``
===========================

.. automodule:: plot_power_network

Creates plots for optimised power network topologies and regional generation,
storage and conversion capacities built.

Rule ``plot_power_network_perfect``
===================================

.. automodule:: plot_power_network_perfect

Creates plots for optimised power network topologies and regional generation,
storage and conversion capacities built for the perfect foresight scenario.

Rule ``plot_hydrogen_network``
==============================

.. automodule:: plot_hydrogen_network

Creates map of optimised hydrogen network, storage and selected other
infrastructure.

Rule ``plot_gas_network``
=========================

.. automodule:: plot_gas_network

Creates map of optimised gas network, storage and selected other
infrastructure.

Rule ``plot_KPIs``
=========================

.. automodule:: plot_KPIs

Plot_KPIs creates both predefined and configurable figures. The predefined figures are:
- A map of curtailed energy for Germany and the wider region.
- A map of transmission line loading for Germany and the wider region.
- An overview plot of the energy trade.
- A summary CSV for storage capacities.

The configurable figures are defined based on ``config/config.kpi.yaml``.

Key components of each configuration:

1. `extract`: Defines the network statistics for the figures (e.g., system cost, generation, storage, emissions).
2. `include` and/or `exclude`: Specifies the regions or entities to include or exclude in the calculation (e.g., "DE" for Germany, "EU" for the European Union).
3. `carrier_filter`: Specifies the energy carriers or technologies (e.g., power, electricity, storage, heat) relevant to the figure.
4. `group_carrier`: Specifies the names used for each carrier. Same names are aggregated.
5. `plot`: Determines the type of plot or visualization to be used (e.g., "detail" for detailed data, "overview" for broader data representation).
6. `figsize`: Determines the figsize of the plots. If not defined, the default size are chosen.
6. `plot_kw`: Additional keyword arguments for the plot (e.g., title, labels, and axis).

Breakdown of the components:

1. **extract**:

  - `system cost`: Extract data from `csvs/nodal_costs.csv`.
  - `capacity`: Extract data from `csvs/nodal_capacities.csv`.
  - `capacity stats`: Extract capacity data from `n.statistics`. See notes for more detail.
  - `generation`: Extract power generation data from `n.statistics.energy_balance`.
  - `emission`: Extract emissions data based on emission links to the atmosphere.
  - `energy balance`: Extract energy balance data from `n.statistics.energy_balance`.
  - `SOC`: Extract state of charge from the `storage_units` component.

3. **carrier_filter**

  a. for system cost, capacity, capacity stats, generation, and emission:
    
    - `electricity`: Filter carrier with AC bus carrier.
    - `electricity+`: Filter carrier with AC bus carrier and also include water tanks and EV batteries.
    - `storage`: Filter carrier with all storage-related technologies (`storage_units`, `links`, and `stores`).
    - `storage-cap`: Filter carrier with all storage capacity-related technologies (`storage_units` and `links`).
    - `storage-energy`: Filter carrier with all storage energy-related technologies (`storage_units` and `stores`).
    - `power`: Filter carrier with power generation technologies.

  b. for energy balance:
    
    - `electricity`: Filter carrier with AC bus carrier.
    - `electricity+`: Filter carrier with AC bus carrier and also include solar rooftops, BEV chargers, and Vehicle-to-Gas.
    - `low voltage`: Filter carrier with low voltage bus carrier.
    - `storage`: Filter carrier with all storage capacity-related technologies.
    - `heat`: Filter carrier with heat bus carrier.
    - `hydrogen`: Filter carrier with hydrogen bus carrier.
  
  c. for SOC can be defined individually or in arrays of `storage_units` carriers.

4. **group_carrier**:

  - `pretty`: The first letter is capitalized, abbreviations are spelled out, and similar names are combined.
  - `sector`: Aggregate all technologies into either the power sector, heating, transport sector, primary fuel, or CCUS.
  - Leaving this option empty will lead to the carriers using `nice_names`.

5. **plot** for system cost, capacity, capacity stats, generation, and emission:

  - `detail`: A single bar plot is created with all the values annotated.
  - `overview`: Multiple bar plots are created.
  - The `plot` for energy balance and SOC is always in a time-series format.

6. **figsize**:

  - For detail plots: the default figsize is (6,8)
  - For overview, energy balance and SOC plots: the default figsize is (12,9)

7. **plot_kw**:

  - For system cost, capacity, capacity stats, generation, and emission, see `matplotlib.pyplot.bar documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.bar.html>`__.
  - For energy balance and SOC, see `pandas.DataFrame.plot documentation <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html#pandas.DataFrame.plot>`__.

