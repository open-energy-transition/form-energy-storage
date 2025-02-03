##########################################
Tutorials
##########################################

This is a brief tutorial on how to use the Form Energy Storage scenarios.
Unlike the original PyPSA-Eur repository, this version comes pre-configured with scenarios that you can test run.

Test run Scenarios
=====================

First, open the ``config.form.yaml`` file using a code editor of your choice, and navigate to the ``run`` section:

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: run:
   :end-before: # docs

as you can see, ``baseline-mds`` is the name of the scenario that is chosen. 
The content of ``baseline-mds`` is defined in ``scenarios.form.yaml``.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: baseline-mds:
   :end-before: baseline-nomds:

.. note::
    If you want to create a model that minimizes computational time, you can adjust the ``clustering: temporal: resolution_sector`` setting here:

    .. literalinclude:: ../config/config.form.yaml
       :language: yaml
       :start-at: clustering:
       :end-before: # docs

    The rule of thumb is:

    The general rule of thumb are as follows:

    * ``4380SEG``: Three-hour time resolution, takes 4 hours to run on a computer cluster
    * ``2920SEG``: Four-hour time resolution, takes 2 hours to run on a computer cluster
    * ``100H``: Hundred-hour time resolution, very quick to solve.

    If you want to test run this model on your computer, replace ``4380SEG`` with ``100H``.

Before running the model, run the snakemake command:

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

To create the model, run the snakemake command:

.. code:: console

   snakemake solve_sector_networks --configfile config/config.form.yaml

This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs:

.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
        0[label = "solve_sector_networks", color = "0.50 0.6 0.85", style="rounded"];
        1[label = "solve_sector_network_myopic", color = "0.43 0.6 0.85", style="rounded"];
        2[label = "final_adjustment_myopic", color = "0.52 0.6 0.85", style="rounded"];
        3[label = "add_existing_baseyear", color = "0.63 0.6 0.85", style="rounded"];
        4[label = "prepare_sector_network", color = "0.29 0.6 0.85", style="rounded"];
        5[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        6[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.40 0.6 0.85", style="rounded"];
        7[label = "build_ship_raster\nrun: baseline-mds", color = "0.41 0.6 0.85", style="rounded"];
        8[label = "retrieve_ship_raster", color = "0.53 0.6 0.85", style="rounded"];
        9[label = "build_shapes\nrun: baseline-mds", color = "0.38 0.6 0.85", style="rounded"];
        10[label = "retrieve_naturalearth_countries", color = "0.65 0.6 0.85", style="rounded,dashed"];
        11[label = "retrieve_eez", color = "0.45 0.6 0.85", style="rounded,dashed"];
        12[label = "retrieve_nuts_shapes", color = "0.04 0.6 0.85", style="rounded,dashed"];
        13[label = "cluster_network\nclusters: 52", color = "0.05 0.6 0.85", style="rounded"];
        14[label = "simplify_network", color = "0.10 0.6 0.85", style="rounded"];
        15[label = "add_transmission_projects_and_dlr", color = "0.60 0.6 0.85", style="rounded"];
        16[label = "base_network", color = "0.34 0.6 0.85", style="rounded"];
        17[label = "retrieve_osm_prebuilt", color = "0.56 0.6 0.85", style="rounded"];
        18[label = "build_transmission_projects", color = "0.15 0.6 0.85", style="rounded"];
        19[label = "build_electricity_demand_base", color = "0.66 0.6 0.85", style="rounded"];
        20[label = "build_electricity_demand\nrun: baseline-mds", color = "0.58 0.6 0.85", style="rounded"];
        21[label = "retrieve_electricity_demand", color = "0.15 0.6 0.85", style="rounded"];
        22[label = "retrieve_synthetic_electricity_demand", color = "0.49 0.6 0.85", style="rounded,dashed"];
        23[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        24[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.40 0.6 0.85", style="rounded"];
        25[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        26[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.40 0.6 0.85", style="rounded"];
        27[label = "cluster_gas_network", color = "0.25 0.6 0.85", style="rounded"];
        28[label = "build_gas_network\nrun: baseline-mds", color = "0.64 0.6 0.85", style="rounded"];
        29[label = "retrieve_gas_infrastructure_data", color = "0.21 0.6 0.85", style="rounded"];
        30[label = "build_gas_input_locations", color = "0.06 0.6 0.85", style="rounded"];
        31[label = "retrieve_gem_europe_gas_tracker", color = "0.35 0.6 0.85", style="rounded,dashed"];
        32[label = "time_aggregation\nsector_opts: ", color = "0.12 0.6 0.85", style="rounded"];
        33[label = "prepare_network\nll: v1.11\nopts: ", color = "0.26 0.6 0.85", style="rounded"];
        34[label = "add_electricity", color = "0.62 0.6 0.85", style="rounded"];
        35[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        36[label = "determine_availability_matrix\ntechnology: solar", color = "0.40 0.6 0.85", style="rounded"];
        37[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        38[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.40 0.6 0.85", style="rounded"];
        39[label = "build_renewable_profiles", color = "0.13 0.6 0.85", style="rounded"];
        40[label = "determine_availability_matrix\ntechnology: onwind", color = "0.40 0.6 0.85", style="rounded"];
        41[label = "build_hydro_profile", color = "0.23 0.6 0.85", style="rounded"];
        42[label = "retrieve_cost_data\nrun: baseline-mds\nyear: 2035", color = "0.18 0.6 0.85", style="rounded"];
        43[label = "build_powerplants", color = "0.31 0.6 0.85", style="rounded"];
        44[label = "build_hourly_heat_demand", color = "0.30 0.6 0.85", style="rounded"];
        45[label = "build_daily_heat_demand", color = "0.54 0.6 0.85", style="rounded"];
        46[label = "build_population_layouts", color = "0.18 0.6 0.85", style="rounded"];
        47[label = "retrieve_worldbank_urban_population", color = "0.30 0.6 0.85", style="rounded,dashed"];
        48[label = "build_solar_thermal_profiles", color = "0.51 0.6 0.85", style="rounded"];
        49[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.39 0.6 0.85", style="rounded"];
        50[label = "build_energy_totals", color = "0.55 0.6 0.85", style="rounded"];
        51[label = "build_clustered_population_layouts", color = "0.04 0.6 0.85", style="rounded"];
        52[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.39 0.6 0.85", style="rounded"];
        53[label = "build_heat_totals", color = "0.39 0.6 0.85", style="rounded"];
        54[label = "build_shipping_demand", color = "0.20 0.6 0.85", style="rounded"];
        55[label = "build_transport_demand", color = "0.01 0.6 0.85", style="rounded"];
        56[label = "build_temperature_profiles", color = "0.25 0.6 0.85", style="rounded"];
        57[label = "build_biomass_potentials\nplanning_horizons: 2035", color = "0.02 0.6 0.85", style="rounded"];
        58[label = "retrieve_jrc_enspreso_biomass", color = "0.02 0.6 0.85", style="rounded,dashed"];
        59[label = "build_salt_cavern_potentials", color = "0.17 0.6 0.85", style="rounded"];
        60[label = "build_industrial_energy_demand_per_node", color = "0.47 0.6 0.85", style="rounded"];
        61[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2035", color = "0.11 0.6 0.85", style="rounded"];
        62[label = "build_industry_sector_ratios", color = "0.49 0.6 0.85", style="rounded"];
        63[label = "build_ammonia_production\nrun: baseline-mds", color = "0.31 0.6 0.85", style="rounded"];
        64[label = "retrieve_usgs_ammonia_production", color = "0.24 0.6 0.85", style="rounded,dashed"];
        65[label = "build_industrial_energy_demand_per_country_today", color = "0.41 0.6 0.85", style="rounded"];
        66[label = "build_industrial_production_per_country", color = "0.32 0.6 0.85", style="rounded"];
        67[label = "build_industrial_production_per_node", color = "0.48 0.6 0.85", style="rounded"];
        68[label = "build_industrial_distribution_key", color = "0.13 0.6 0.85", style="rounded"];
        69[label = "retrieve_hotmaps_industrial_sites", color = "0.03 0.6 0.85", style="rounded,dashed"];
        70[label = "retrieve_gem_steel_plant_tracker", color = "0.36 0.6 0.85", style="rounded,dashed"];
        71[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2035", color = "0.06 0.6 0.85", style="rounded"];
        72[label = "build_industrial_energy_demand_per_node_today", color = "0.03 0.6 0.85", style="rounded"];
        73[label = "build_district_heat_share\nplanning_horizons: 2035", color = "0.52 0.6 0.85", style="rounded"];
        74[label = "build_cop_profiles", color = "0.61 0.6 0.85", style="rounded"];
        75[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2035", color = "0.16 0.6 0.85", style="rounded"];
        76[label = "build_direct_heat_source_utilisation_profiles", color = "0.59 0.6 0.85", style="rounded"];
        77[label = "build_existing_heating_distribution", color = "0.08 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        42 -> 1
        3 -> 2
        4 -> 3
        43 -> 3
        14 -> 3
        13 -> 3
        51 -> 3
        42 -> 3
        74 -> 3
        77 -> 3
        50 -> 3
        5 -> 4
        23 -> 4
        25 -> 4
        27 -> 4
        30 -> 4
        32 -> 4
        33 -> 4
        49 -> 4
        52 -> 4
        54 -> 4
        55 -> 4
        44 -> 4
        50 -> 4
        57 -> 4
        42 -> 4
        59 -> 4
        14 -> 4
        13 -> 4
        51 -> 4
        60 -> 4
        67 -> 4
        73 -> 4
        56 -> 4
        74 -> 4
        48 -> 4
        76 -> 4
        6 -> 5
        9 -> 5
        13 -> 5
        7 -> 6
        9 -> 6
        13 -> 6
        8 -> 7
        10 -> 9
        11 -> 9
        12 -> 9
        14 -> 13
        19 -> 13
        15 -> 14
        16 -> 14
        16 -> 15
        18 -> 15
        17 -> 16
        9 -> 16
        16 -> 18
        9 -> 18
        14 -> 19
        9 -> 19
        20 -> 19
        21 -> 20
        22 -> 20
        24 -> 23
        9 -> 23
        13 -> 23
        7 -> 24
        9 -> 24
        13 -> 24
        26 -> 25
        9 -> 25
        13 -> 25
        7 -> 26
        9 -> 26
        13 -> 26
        28 -> 27
        13 -> 27
        29 -> 28
        31 -> 30
        29 -> 30
        13 -> 30
        33 -> 32
        44 -> 32
        48 -> 32
        34 -> 33
        42 -> 33
        35 -> 34
        37 -> 34
        39 -> 34
        5 -> 34
        23 -> 34
        25 -> 34
        41 -> 34
        13 -> 34
        42 -> 34
        43 -> 34
        19 -> 34
        36 -> 35
        9 -> 35
        13 -> 35
        9 -> 36
        13 -> 36
        38 -> 37
        9 -> 37
        13 -> 37
        9 -> 38
        13 -> 38
        40 -> 39
        9 -> 39
        13 -> 39
        9 -> 40
        13 -> 40
        9 -> 41
        13 -> 43
        45 -> 44
        46 -> 45
        13 -> 45
        9 -> 46
        47 -> 46
        46 -> 48
        13 -> 48
        50 -> 49
        51 -> 49
        9 -> 50
        46 -> 51
        13 -> 51
        53 -> 52
        51 -> 52
        50 -> 53
        9 -> 54
        13 -> 54
        50 -> 54
        51 -> 55
        49 -> 55
        50 -> 55
        56 -> 55
        46 -> 56
        13 -> 56
        58 -> 57
        12 -> 57
        13 -> 57
        9 -> 57
        13 -> 59
        61 -> 60
        67 -> 60
        72 -> 60
        62 -> 61
        65 -> 61
        66 -> 61
        63 -> 62
        64 -> 63
        50 -> 65
        66 -> 65
        63 -> 66
        68 -> 67
        71 -> 67
        13 -> 68
        51 -> 68
        69 -> 68
        70 -> 68
        66 -> 71
        68 -> 72
        65 -> 72
        50 -> 73
        51 -> 73
        75 -> 74
        56 -> 74
        13 -> 74
        56 -> 75
        13 -> 75
        75 -> 76
        51 -> 77
        49 -> 77
        73 -> 77
    }
|

To generate the KPI plots, run the snakemake command:

.. code:: console

   snakemake plot_KPIs_all --configfile config/config.form.yaml


.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
        0[label = "plot_KPIs_all", color = "0.55 0.6 0.85", style="rounded"];
        1[label = "plot_KPIs", color = "0.21 0.6 0.85", style="rounded"];
        2[label = "solve_sector_network_myopic", color = "0.18 0.6 0.85", style="rounded"];
        3[label = "final_adjustment_myopic", color = "0.55 0.6 0.85", style="rounded"];
        4[label = "add_existing_baseyear", color = "0.56 0.6 0.85", style="rounded"];
        5[label = "prepare_sector_network", color = "0.16 0.6 0.85", style="rounded"];
        6[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        7[label = "determine_availability_matrix\ntechnology: offwind-ac", color = "0.44 0.6 0.85", style="rounded"];
        8[label = "build_ship_raster\nrun: baseline-mds", color = "0.63 0.6 0.85", style="rounded"];
        9[label = "retrieve_ship_raster", color = "0.64 0.6 0.85", style="rounded"];
        10[label = "build_shapes\nrun: baseline-mds", color = "0.34 0.6 0.85", style="rounded"];
        11[label = "retrieve_naturalearth_countries", color = "0.29 0.6 0.85", style="rounded,dashed"];
        12[label = "retrieve_eez", color = "0.16 0.6 0.85", style="rounded,dashed"];
        13[label = "retrieve_nuts_shapes", color = "0.03 0.6 0.85", style="rounded,dashed"];
        14[label = "cluster_network\nclusters: 52", color = "0.49 0.6 0.85", style="rounded"];
        15[label = "simplify_network", color = "0.04 0.6 0.85", style="rounded"];
        16[label = "add_transmission_projects_and_dlr", color = "0.61 0.6 0.85", style="rounded"];
        17[label = "base_network", color = "0.10 0.6 0.85", style="rounded"];
        18[label = "retrieve_osm_prebuilt", color = "0.06 0.6 0.85", style="rounded"];
        19[label = "build_transmission_projects", color = "0.57 0.6 0.85", style="rounded"];
        20[label = "build_electricity_demand_base", color = "0.26 0.6 0.85", style="rounded"];
        21[label = "build_electricity_demand\nrun: baseline-mds", color = "0.43 0.6 0.85", style="rounded"];
        22[label = "retrieve_electricity_demand", color = "0.51 0.6 0.85", style="rounded"];
        23[label = "retrieve_synthetic_electricity_demand", color = "0.11 0.6 0.85", style="rounded,dashed"];
        24[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        25[label = "determine_availability_matrix\ntechnology: offwind-dc", color = "0.44 0.6 0.85", style="rounded"];
        26[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        27[label = "determine_availability_matrix\ntechnology: offwind-float", color = "0.44 0.6 0.85", style="rounded"];
        28[label = "cluster_gas_network", color = "0.14 0.6 0.85", style="rounded"];
        29[label = "build_gas_network\nrun: baseline-mds", color = "0.07 0.6 0.85", style="rounded"];
        30[label = "retrieve_gas_infrastructure_data", color = "0.22 0.6 0.85", style="rounded"];
        31[label = "build_gas_input_locations", color = "0.04 0.6 0.85", style="rounded"];
        32[label = "retrieve_gem_europe_gas_tracker", color = "0.17 0.6 0.85", style="rounded,dashed"];
        33[label = "time_aggregation\nsector_opts: ", color = "0.65 0.6 0.85", style="rounded"];
        34[label = "prepare_network\nll: v1.11\nopts: ", color = "0.12 0.6 0.85", style="rounded"];
        35[label = "add_electricity", color = "0.39 0.6 0.85", style="rounded"];
        36[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        37[label = "determine_availability_matrix\ntechnology: solar", color = "0.44 0.6 0.85", style="rounded"];
        38[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        39[label = "determine_availability_matrix\ntechnology: solar-hsat", color = "0.44 0.6 0.85", style="rounded"];
        40[label = "build_renewable_profiles", color = "0.58 0.6 0.85", style="rounded"];
        41[label = "determine_availability_matrix\ntechnology: onwind", color = "0.44 0.6 0.85", style="rounded"];
        42[label = "build_hydro_profile", color = "0.44 0.6 0.85", style="rounded"];
        43[label = "retrieve_cost_data\nrun: baseline-mds\nyear: 2035", color = "0.57 0.6 0.85", style="rounded"];
        44[label = "build_powerplants", color = "0.31 0.6 0.85", style="rounded"];
        45[label = "build_hourly_heat_demand", color = "0.60 0.6 0.85", style="rounded"];
        46[label = "build_daily_heat_demand", color = "0.25 0.6 0.85", style="rounded"];
        47[label = "build_population_layouts", color = "0.08 0.6 0.85", style="rounded"];
        48[label = "retrieve_worldbank_urban_population", color = "0.66 0.6 0.85", style="rounded,dashed"];
        49[label = "build_solar_thermal_profiles", color = "0.22 0.6 0.85", style="rounded"];
        50[label = "build_population_weighted_energy_totals\nkind: energy", color = "0.27 0.6 0.85", style="rounded"];
        51[label = "build_energy_totals", color = "0.41 0.6 0.85", style="rounded"];
        52[label = "build_clustered_population_layouts", color = "0.26 0.6 0.85", style="rounded"];
        53[label = "build_population_weighted_energy_totals\nkind: heat", color = "0.27 0.6 0.85", style="rounded"];
        54[label = "build_heat_totals", color = "0.52 0.6 0.85", style="rounded"];
        55[label = "build_shipping_demand", color = "0.45 0.6 0.85", style="rounded"];
        56[label = "build_transport_demand", color = "0.00 0.6 0.85", style="rounded"];
        57[label = "build_temperature_profiles", color = "0.02 0.6 0.85", style="rounded"];
        58[label = "build_biomass_potentials\nplanning_horizons: 2035", color = "0.12 0.6 0.85", style="rounded"];
        59[label = "retrieve_jrc_enspreso_biomass", color = "0.18 0.6 0.85", style="rounded,dashed"];
        60[label = "build_salt_cavern_potentials", color = "0.48 0.6 0.85", style="rounded"];
        61[label = "build_industrial_energy_demand_per_node", color = "0.13 0.6 0.85", style="rounded"];
        62[label = "build_industry_sector_ratios_intermediate\nplanning_horizons: 2035", color = "0.58 0.6 0.85", style="rounded"];
        63[label = "build_industry_sector_ratios", color = "0.54 0.6 0.85", style="rounded"];
        64[label = "build_ammonia_production\nrun: baseline-mds", color = "0.01 0.6 0.85", style="rounded"];
        65[label = "retrieve_usgs_ammonia_production", color = "0.11 0.6 0.85", style="rounded,dashed"];
        66[label = "build_industrial_energy_demand_per_country_today", color = "0.64 0.6 0.85", style="rounded"];
        67[label = "build_industrial_production_per_country", color = "0.23 0.6 0.85", style="rounded"];
        68[label = "build_industrial_production_per_node", color = "0.59 0.6 0.85", style="rounded"];
        69[label = "build_industrial_distribution_key", color = "0.46 0.6 0.85", style="rounded"];
        70[label = "retrieve_hotmaps_industrial_sites", color = "0.50 0.6 0.85", style="rounded,dashed"];
        71[label = "retrieve_gem_steel_plant_tracker", color = "0.24 0.6 0.85", style="rounded,dashed"];
        72[label = "build_industrial_production_per_country_tomorrow\nplanning_horizons: 2035", color = "0.39 0.6 0.85", style="rounded"];
        73[label = "build_industrial_energy_demand_per_node_today", color = "0.15 0.6 0.85", style="rounded"];
        74[label = "build_district_heat_share\nplanning_horizons: 2035", color = "0.45 0.6 0.85", style="rounded"];
        75[label = "build_cop_profiles", color = "0.49 0.6 0.85", style="rounded"];
        76[label = "build_central_heating_temperature_profiles\nplanning_horizons: 2035", color = "0.21 0.6 0.85", style="rounded"];
        77[label = "build_direct_heat_source_utilisation_profiles", color = "0.59 0.6 0.85", style="rounded"];
        78[label = "build_existing_heating_distribution", color = "0.33 0.6 0.85", style="rounded"];
        79[label = "make_summary", color = "0.62 0.6 0.85", style="rounded"];
        80[label = "plot_power_network_clustered", color = "0.53 0.6 0.85", style="rounded"];
        81[label = "plot_power_network", color = "0.25 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        14 -> 1
        79 -> 1
        3 -> 2
        43 -> 2
        4 -> 3
        5 -> 4
        44 -> 4
        15 -> 4
        14 -> 4
        52 -> 4
        43 -> 4
        75 -> 4
        78 -> 4
        51 -> 4
        6 -> 5
        24 -> 5
        26 -> 5
        28 -> 5
        31 -> 5
        33 -> 5
        34 -> 5
        50 -> 5
        53 -> 5
        55 -> 5
        56 -> 5
        45 -> 5
        51 -> 5
        58 -> 5
        43 -> 5
        60 -> 5
        15 -> 5
        14 -> 5
        52 -> 5
        61 -> 5
        68 -> 5
        74 -> 5
        57 -> 5
        75 -> 5
        49 -> 5
        77 -> 5
        7 -> 6
        10 -> 6
        14 -> 6
        8 -> 7
        10 -> 7
        14 -> 7
        9 -> 8
        11 -> 10
        12 -> 10
        13 -> 10
        15 -> 14
        20 -> 14
        16 -> 15
        17 -> 15
        17 -> 16
        19 -> 16
        18 -> 17
        10 -> 17
        17 -> 19
        10 -> 19
        15 -> 20
        10 -> 20
        21 -> 20
        22 -> 21
        23 -> 21
        25 -> 24
        10 -> 24
        14 -> 24
        8 -> 25
        10 -> 25
        14 -> 25
        27 -> 26
        10 -> 26
        14 -> 26
        8 -> 27
        10 -> 27
        14 -> 27
        29 -> 28
        14 -> 28
        30 -> 29
        32 -> 31
        30 -> 31
        14 -> 31
        34 -> 33
        45 -> 33
        49 -> 33
        35 -> 34
        43 -> 34
        36 -> 35
        38 -> 35
        40 -> 35
        6 -> 35
        24 -> 35
        26 -> 35
        42 -> 35
        14 -> 35
        43 -> 35
        44 -> 35
        20 -> 35
        37 -> 36
        10 -> 36
        14 -> 36
        10 -> 37
        14 -> 37
        39 -> 38
        10 -> 38
        14 -> 38
        10 -> 39
        14 -> 39
        41 -> 40
        10 -> 40
        14 -> 40
        10 -> 41
        14 -> 41
        10 -> 42
        14 -> 44
        46 -> 45
        47 -> 46
        14 -> 46
        10 -> 47
        48 -> 47
        47 -> 49
        14 -> 49
        51 -> 50
        52 -> 50
        10 -> 51
        47 -> 52
        14 -> 52
        54 -> 53
        52 -> 53
        51 -> 54
        10 -> 55
        14 -> 55
        51 -> 55
        52 -> 56
        50 -> 56
        51 -> 56
        57 -> 56
        47 -> 57
        14 -> 57
        59 -> 58
        13 -> 58
        14 -> 58
        10 -> 58
        14 -> 60
        62 -> 61
        68 -> 61
        73 -> 61
        63 -> 62
        66 -> 62
        67 -> 62
        64 -> 63
        65 -> 64
        51 -> 66
        67 -> 66
        64 -> 67
        69 -> 68
        72 -> 68
        14 -> 69
        52 -> 69
        70 -> 69
        71 -> 69
        67 -> 72
        69 -> 73
        66 -> 73
        51 -> 74
        52 -> 74
        76 -> 75
        57 -> 75
        14 -> 75
        57 -> 76
        14 -> 76
        76 -> 77
        52 -> 78
        50 -> 78
        74 -> 78
        2 -> 79
        43 -> 79
        80 -> 79
        81 -> 79
        14 -> 80
        2 -> 81
        14 -> 81
    }       
|

Repository structure
=====================

* ``benchmarks``: will store ``snakemake`` benchmarks (does not exist initially)
* ``config``: configurations used in the study
* ``cutouts``: will store raw weather data cutouts from ``atlite`` (does not exist initially)
* ``data``: includes input data that is not produced by any ``snakemake`` rule
* ``doc``: includes all files necessary to build the ``readthedocs`` documentation of PyPSA-Eur
* ``envs``: includes all the ``mamba`` environment specifications to run the workflow
* ``logs``: will store log files (does not exist initially)
* ``notebooks``: includes all the ``notebooks`` used for ad-hoc analysis
* ``report``: contains all files necessary to build the report; plots and result files are generated automatically
* ``rules``: includes all the ``snakemake``rules loaded in the ``Snakefile``
* ``resources``: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* ``results``: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* ``scripts``: includes all the Python scripts executed by the ``snakemake`` rules to build the model