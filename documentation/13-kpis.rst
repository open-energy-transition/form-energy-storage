##########################################
KPIs Configuration
##########################################

Plot_KPIs creates both predefined and configurable figures. The predefined figures are:

- A map of curtailed energy for Germany and the wider region.
- A map of transmission line loading for Germany and the wider region.
- An overview plot of the energy trade.
- A summary CSV for storage capacities.

The configurable figures are defined based on ``config/config.kpi.yaml``.

Key components of each configuration:

1. ``extract``: Defines the network statistics for the figures (e.g., system cost, generation, storage, emissions).
2. ``include`` and/or ``exclude``: Specifies the regions or entities to include or exclude in the calculation (e.g., "DE" for Germany, "EU" for the European Union).
3. ``carrier_filter``: Specifies the energy carriers or technologies (e.g., power, electricity, storage, heat) relevant to the figure.
4. ``group_carrier``: Specifies the names used for each carrier. Same names are aggregated.
5. ``plot``: Determines the type of plot or visualization to be used (e.g., "detail" for detailed data, "overview" for broader data representation).
6. ``figsize``: Determines the figsize of the plots. If not defined, the default size are chosen.
7. ``plot_kw``: Additional keyword arguments for the plot (e.g., title, labels, and axis).

Breakdown of the components:

``extract``
------------------------

- **system cost**: Extract data from ``csvs/nodal_costs.csv``.
- **capacity**: Extract data from ``csvs/nodal_capacities.csv``.
- **capacity stats**: Extract capacity data from ``n.statistics``. See notes for more detail.
- **generation**: Extract power generation data from ``n.statistics.energy_balance``.
- **emission**: Extract emissions data based on emission links to the atmosphere.
- **energy balance**: Extract energy balance data from ``n.statistics.energy_balance``.
- **SOC**: Extract state of charge from the ``storage_units`` component.

.. note::
    **capacity stats** is a special option because it has two extra components:

    ``stats``: Three options:
    
    - **install**: Extract capacity data from ``n.statistics.installed_capacity``
    - **optimal**: Extract capacity data from ``n.statistics.optimal_capacity``
    - **expand**:  The difference between optimal minus install.

    ``storage``: If set to ``true``, only storage capacities of the component ``store`` and ``storage_units`` are taken into account.

``carrier_filter``
------------------------

a. for system cost, capacity, capacity stats, generation, and emission:

- **electricity**: Filter carrier with AC bus carrier.
- **electricity+**: Filter carrier with AC bus carrier and also include water tanks and EV batteries.
- **storage**: Filter carrier with all storage-related technologies (``storage_units``, ``links``, and ``stores``).
- **storage-cap**: Filter carrier with all storage capacity-related technologies (``storage_units`` and ``links``).
- **storage-energy**: Filter carrier with all storage energy-related technologies (``storage_units`` and ``stores``).
- **power**: Filter carrier with power generation technologies.

.. note::
    **storage**, **storage-cap**, **storage-energy**, **power** is based on a list in ``config/config.kpi.yaml``.
    You can create your own filtering scheme by adding a list in ``kpi: filter_scheme``

    .. literalinclude:: ../config/config.kpi.yaml
       :language: yaml
       :start-at: power:
       :end-at: ]

b. for energy balance:

- **electricity**: Filter carrier with AC bus carrier.
- **electricity+**: Filter carrier with AC bus carrier and also include solar rooftops, BEV chargers, and Vehicle-to-Gas.
- **low voltage**: Filter carrier with low voltage bus carrier.
- **storage**: Filter carrier with all storage capacity-related technologies.
- **heat**: Filter carrier with heat bus carrier.
- **hydrogen**: Filter carrier with hydrogen bus carrier.

c. for SOC can be defined individually or in arrays of ``storage_units`` carriers.

``group_carrier``
------------------------

- **pretty**: The first letter is capitalized, abbreviations are spelled out, and similar names are combined.
- **sector**: Aggregate all technologies into either the power sector, heating, transport sector, primary fuel, or CCUS.

Leaving this option empty will lead to the carriers using ``nice_names``.

``plot``
------------------------

for system cost, capacity, capacity stats, generation, and emission:

- **detail**: A single bar plot is created with all the values annotated.
- **overview**: Multiple bar plots are created.

The plot for energy balance and SOC is always in a time-series format.

``figsize``
------------------------

- For detail plots: the default figsize is (6,8)
- For overview, energy balance and SOC plots: the default figsize is (12,9)

``plot_kw``
------------------------

- For system cost, capacity, capacity stats, generation, and emission, see `matplotlib.pyplot.bar documentation <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.bar.html>`__.
- For energy balance and SOC, see `pandas.DataFrame.plot documentation <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.html#pandas.DataFrame.plot>`__.

