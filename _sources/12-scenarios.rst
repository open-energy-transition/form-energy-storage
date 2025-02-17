##########################################
Scenarios Configuration
##########################################

This is the scenario outline used in the analysis **The Role of Energy Storage in Germany**.
For each scenario involving multi-day storages, 
there is always a corresponding counterfactual scenario that excludes multi-day storages.

.. note::
   The PyPSA-Eur configuration files follow a pyramid-like structure, where the parameters in the highest configuration file add to and override those in the configuration file below it. 
   The order is as follows:

   1. scenarios.form.yaml (**this section**)
   2. `config.form.yaml <https://open-energy-transition.github.io/form-energy-storage/11-baseline.html>`_
   3. `config.default.yaml <https://pypsa-eur.readthedocs.io/en/latest/configuration.html>`_
   4. `config.kpi.yaml <https://open-energy-transition.github.io/form-energy-storage/13-kpis.html>`_

Thus, for example changes specified in `scenario.form.yaml` will add to and override configurations in `config.form.yaml` and so on.

Thus, for example changes specified in `scenario.form.yaml` will add to and override configurations in `config.form.yaml` and so on.

``Baseline``
============

The Baseline scenario is based on the climate year 2013, with the cost of iron-air set at 2,350 EUR/kWh (2024).

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: baseline-mds:
   :end-before: mid-capex-mds:

``Mid-Capex``
=============

The Mid-Capex scenario sets the cost of iron-air at 2,000 EUR/kWh (2024), 
slightly lower than the Baseline scenario.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: mid-capex-mds:
   :end-before: low-capex-mds:

``Low-Capex``
=============

The Low-Capex scenario sets the cost of iron-air at 1,525 EUR/kWh (2024), 
representing the lowest cost for iron-air when compared to the Baseline scenario.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: low-capex-mds:
   :end-before: cy2010-mds:

``Climate Year 2010``
=====================

This scenario is based on the climate year 2010, which is characterised by the lowest annual wind resources 
and the second highest annual heating demand from 1960 - 2021 (`Gøtske et al., 2024 <https://www.nature.com/articles/s41467-024-54853-3>`_)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy2010-mds:
   :end-before: cy2012-mds:

``Climate Year 2012``
=====================

This scenario is based on the climate year 2012, which is the most stressful climate year in terms of 
2-week Dunkelflaute situation at European aggregated level (`2024 TYNDP Scenarios Methodology Report <https://2024.entsos-tyndp-scenarios.eu/wp-content/uploads/2024/07/TYNDP_2024_Scenarios_Methodology_Report_240708.pdf>`_, page 30)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy2012-mds:
   :end-before: cy1996-mds:

``Climate Year 1996``
=====================

This scenario is based on the climate year 1996, which is characterised by low yield of PV and wind, 
as well as a pronounced “dark doldrums” for about 10 days in January (`TransnetBW Energy System 2050 Report <https://www.energysystem2050.net/content/TransnetBW-Study_EnergySystem2050.pdf?v2>`_, page 19)

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: cy1996-mds:
   :end-before: nogasppde-mid-capex-mds:

``Germany Gas Power Plant Phase-Out``
==========================

These scenarios assume that Germany, driven by policy decision, successfully phases out gas-fired power plants and gas CHPs in 2035. 
The impact of the cost difference between iron-air batteries at 2,000 EUR/kWh (2024) and 1,525 EUR/kWh (2024) 
is also examined.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: nogasppde-mid-capex-mds:
   :end-before: delayed-nep-mid-capex-mds:

``Transmission Sensitivity``
============================

These scenarios assume that Germany's 'NEP' transmission projects face delays, with only confirmed AC lines and DC links available that are to be built before and including 2032. 
The impact of the cost difference between iron-air batteries at 2,000 EUR/kWh (2024) and 1,525 EUR/kWh (2024) 
is also examined.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: delayed-nep-mid-capex-mds:
   :end-before: gasprice-peak-mid-capex-mds:

``High Gas Price``
==================

These scenarios assume that the EU's gas market faces a 50% price increase compared to the baseline scenario.
The impact of the cost difference between iron-air batteries at 2,000 EUR/kWh (2024) and 1,525 EUR/kWh (2024) 
is also examined.

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: gasprice-peak-mid-capex-mds:
   :end-before: gasprice-peak-nomds: