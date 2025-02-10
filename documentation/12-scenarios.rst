##########################################
Scenarios Configuration
##########################################

This is the scenario outline used in the analysis **The Role of Energy Storage in Germany**.
For each scenario involving medium-duration storages, 
there is always a corresponding counterfactual scenario that excludes medium-duration storages.

``Baseline``
=============

The Baseline scenario is based on the climate year 2013, with the cost of iron-air set at 2,350 EUR/kWh.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: baseline-mds:
   :end-before: mid-capex-mds:

``Mid-Capex``
=============

The Mid-Capex scenario sets the cost of iron-air at 2,000 EUR/kWh, 
slightly lower than the Baseline scenario.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: mid-capex-mds:
   :end-before: low-capex-mds:

``Low-Capex``
=============

The Low-Capex scenario sets the cost of iron-air at 1,525 EUR/kWh, 
which is half the cost of the Baseline scenario.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: low-capex-mds:
   :end-before: cy2010-mds:

``Climate Year 2010``
==========================

This scenario is based on the climate year 2010, which is characterised by the lowest annual wind resources 
and the second highest annual heating demand from 1960 - 2021 (`Gøtske et al., 2024 <https://www.nature.com/articles/s41467-024-54853-3>`_)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy2010-mds:
   :end-before: cy2012-mds:

``Climate Year 2012``
==========================

This scenario is based on the climate year 2012, which is the most stressful climate year in terms of 
2-week Dunkelflaute situation at European aggregated level (`2024 TYNDP Scenarios Methodology Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2024/07/TYNDP_2024_Scenarios_Methodology_Report_240708.pdf>`_, page 30)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy2012-mds:
   :end-before: cy1996-mds:

``Climate Year 1996``
==========================

This scenario is based on the climate year 1996, which is characterised by low yield of PV and wind, 
as well as a pronounced “dark doldrums” for about 10 days in January (`TransnetBW Energy System 2050 Report <https://www.energysystem2050.net/content/TransnetBW-Study_EnergySystem2050.pdf?v2>`_, page 19)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy1996-mds: