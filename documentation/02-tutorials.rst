##########################################
Tutorials
##########################################

This is a brief tutorial on how to recreate the scenarios of the study **The Role of Energy Storage in Germany**.
Unlike the original PyPSA-Eur repository, this version comes pre-configured with scenarios that you can test run.

Build the tutorial model
=========================

First, open the base configuration file for this study, called ``config.form.yaml`, using a code editor of your choice, and navigate to the ``run`` section:

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: run:
   :end-before: # docs

as you can see, ``baseline-mds`` is the name of the scenario that is chosen. 
The content of ``baseline-mds`` as well as a comprehensive list of all studied scenarios is defined in ``scenarios.form.yaml`` as described in `Scenarios Configuration <https://open-energy-transition.github.io/form-energy-storage/12-scenarios.html>`_.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: baseline-mds:
   :end-before: baseline-nomds:

To create your own scenarios, it's recommended to define the configurations in ``scenarios.form.yaml`` by either modifying an existing scenario or duplicating and renaming a scenario.

.. note::
   If you want to create a model with lower computational requirements, you can adjust the ``clustering: temporal: resolution_sector`` setting here:

   .. literalinclude:: ../config/config.form.yaml
      :language: yaml
      :start-at: clustering:
      :end-before: # docs

   The general rules of thumb are as follows:

   * ``1H``: 1-hour time resolution
   * ``4380SEG``: 3 times faster computational time compared to ``1H``
   * ``2920SEG``: 4 times faster computational time compared to ``1H``
   * ``100H``: significantly quicker to solve compared to ``1H``
   
   If you want to test run this model on your computer, replace ``4380SEG`` with ``100H``.
   Using a lower time resolution, such as that from ``SEG``, can still capture the desired variability in the model. 
   See `Segmentation Temporal Clustering <https://open-energy-transition.github.io/form-energy-storage/21-segmentation.html>`_ for more details

Run the model
=========================

Before running the model, it's recommended to test run the snakemake command:

.. code:: console

   snakemake solve_sector_networks --configfile config/config.form.yaml -n

The ``-n`` at the end indicates a dry run. If no errors occur, it will list the jobs that would be run:

.. code:: console

    Job stats:
    job                                                 count
    ------------------------------------------------  -------
    add_electricity                                         1
    add_existing_baseyear                                   1
    add_transmission_projects_and_dlr                       1
    base_network                                            1
    build_ammonia_production                                1
    build_biomass_potentials                                1
    build_central_heating_temperature_profiles              1
    build_clustered_population_layouts                      1
    build_cop_profiles                                      1
    build_daily_heat_demand                                 1
    build_direct_heat_source_utilisation_profiles           1
    build_district_heat_share                               1
    build_electricity_demand                                1
    build_electricity_demand_base                           1
    build_energy_totals                                     1
    build_existing_heating_distribution                     1
    build_gas_input_locations                               1
    build_gas_network                                       1
    build_heat_totals                                       1
    build_hourly_heat_demand                                1
    build_hydro_profile                                     1
    build_industrial_distribution_key                       1
    build_industrial_energy_demand_per_country_today        1
    build_industrial_energy_demand_per_node                 1
    build_industrial_energy_demand_per_node_today           1
    build_industrial_production_per_country                 1
    build_industrial_production_per_country_tomorrow        1
    build_industrial_production_per_node                    1
    build_industry_sector_ratios                            1
    build_industry_sector_ratios_intermediate               1
    build_population_layouts                                1
    build_population_weighted_energy_totals                 2
    build_powerplants                                       1
    build_renewable_profiles                                6
    build_salt_cavern_potentials                            1
    build_shapes                                            1
    build_ship_raster                                       1
    build_shipping_demand                                   1
    build_solar_thermal_profiles                            1
    build_temperature_profiles                              1
    build_transmission_projects                             1
    build_transport_demand                                  1
    cluster_gas_network                                     1
    cluster_network                                         1
    determine_availability_matrix                           6
    final_adjustment_myopic                                 1
    prepare_network                                         1
    prepare_sector_network                                  1
    retrieve_cost_data                                      1
    retrieve_electricity_demand                             1
    retrieve_gas_infrastructure_data                        1
    retrieve_osm_prebuilt                                   1
    retrieve_ship_raster                                    1
    simplify_network                                        1
    solve_sector_network_myopic                             1
    solve_sector_networks                                   1
    time_aggregation                                        1
    total                                                  68

.. note::
    If an error occurs due to missing data, the first solution is to check the ``enable:`` setting in ``config.form.yaml``.

    .. literalinclude:: ../config/config.form.yaml
       :language: yaml
       :start-after: #enable
       :end-before: # docs

    * ``retrieve_databundle``: set this to true if a file in ``data/`` is missing
    * ``retrieve_cost_data``: set this to true if ``resources/{run}/cost_2035.csv`` is missing
    * ``retrieve_cutout``: set this to true if ``cutouts/europe-{weather year}-sarah3-era5.nc`` is missing

    The `weather data cutouts for PyPSA-Eur <https://zenodo.org/records/14936211>`_ are available for all the weather years included in the scenarios.
    If the weather year you've selected is not on that list, you'll need to set ``build_cutout`` to true.

To create and solve the model, run the snakemake command:

.. code:: console

   snakemake solve_sector_networks --configfile config/config.form.yaml

This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs.
To create and solve the model and automatically generate all KPI plots of the study, run the snakemake command:

.. code:: console
   
   snakemake all --configfile config/config.form.yaml

Analyze the model
=========================

There are several ways to analyze the results of the Snakemake run:

* By analyzing the PyPSA network directly, see `Package_unpackage_networks.ipynb <https://github.com/open-energy-transition/form-energy-storage/tree/master/notebooks/Package_unpackage_networks.ipynb>`_ for more details.
* By using the CSV files and graphs located in ``results/{run}/csvs/`` and ``results/{run}/graphs/``
* By using a pre-configured Jupyter Notebook in ``notebooks/`` or any Jupyter Notebook of your choice.
* By utilizing the customizable KPIs in ``results/{run}/maps/``

The automatic generation of KPI plots is a feature from this project and may be upstreamed in the future. To create your own custom KPI plots, see `KPI Configuration <https://open-energy-transition.github.io/form-energy-storage/13-kpis.html>`_

To generate the KPI plots, run the snakemake command:

.. code:: console

   snakemake plot_KPIs_all --configfile config/config.form.yaml

The Jupyter Notebook ``notebooks/Scenarios_comparison.ipynb`` uses the script functions in ``plot_KPIs`` to plot and compare the results of multiple scenarios.
