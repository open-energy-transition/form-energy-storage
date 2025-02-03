##########################################
Base Configuration
##########################################

This is the base configuration outline used in the analysis "The Role of Energy Storage in Germany".
The complete explaination of each configuration can be found in `PyPSA-Eur: Configuration <https://pypsa-eur.readthedocs.io/en/latest/configuration.html>`_.


``run``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#run>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: run:
   :end-before: # docs

This is the only few section in the file that needs to be changed in order to run the scenario.

* prefix: (optional) directory for output results
* name: Scenario name from ``config/scenarios.form.yaml``

``foresight``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#foresight>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: foresight:
   :end-at: foresight:

The scope of this work is based on myopic foresight.

``scenario``
============

* `Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#scenario>`_
* `Wildcard Documentation <https://pypsa-eur.readthedocs.io/en/latest/wildcards.html>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: scenario:
   :end-before: # docs

* **ll**: v1.11 indicates that the existing transmission network capacity can be increased by 11%.
* **clusters**: The model results are grouped into 52 nodes; see the clustering section for the disaggregation list.
* **planning_horizons**: The model is simulating the year 2035.

``countries``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#countries>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: countries:
   :end-before: # docs

The analysis includes Germany and its neighboring countries, along with Italy. The complete list:

* AT: Austria
* BE: Belgium
* CH: Switzerland
* CZ: Czech Republic
* DE: Germany
* DK: Denmark
* FR: France
* IT: Italy
* LU: Luxembourg
* NL: Netherlands
* PL: Poland
* SE: Sweden

``snapshots``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#snapshots>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: snapshots:
   :end-before: # docs

The baseline scenario is based on the climate year 2013, with the snapshots adjusted accordingly.

``enable``
==========

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#enable>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-after: #enable
   :end-before: # docs

For the first run, it is recommended to set ``retrieve_databundle``, ``retrieve_cost_data`` and ``retrieve_cutout`` as true.

``co2 budget``
==============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#co2-budget>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: co2_budget:
   :end-before: # docs

The CO2 budget is derived from individual countries carbon emission targets for included sectors excluding domestic transport emissions from scope
This means that the total budget is: 

.. math::
   1.487e+9 * 0.107 = 

The complete calculation can be seen in the Assumption processing section.

``electricity``
===============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#electricity>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: electricity:
   :end-before: # docs

Configuration changes made:

* Lithium-ion are disaggregated into 1h, 2h, 4h and 8h for a more detail analysis.
* Iron-air batteries have 100h max hours
* Nuclear, Hard Coal, and Lignite power plants are excluded from the powerplantmatching database.
* Instead, we use a validated dataset of power plants from ``data/custom_powerplants.csv``

``atlite``
==========

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#atlite>`_

Define and specify the ``atlite.Cutout`` used for calculating renewable potentials and time-series. All options except for ``features`` are directly used as `cutout parameters <https://atlite.readthedocs.io/en/latest/ref_api.html#cutout>`__.

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: atlite:
   :end-before: # docs

The ``default_cutout`` is derived from SARAH-3 weather data, with missing values filled in using ERA-5 weather data.

``renewable``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#renewable>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: renewable:
   :end-before: # docs

The modification in max hours of pumped hydro storages are based on Form Energy. 

``conventional``
================

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#conventional>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: conventional:
   :end-before: # docs

The minimum energy dispatch for nuclear power plants is set to 50%, instead of the default 20%.

``lines``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#lines>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: lines:
   :end-before: # docs

The maximum line capacity is set to 50% of its total capacity, rather than 70%, 
to better approximate security and reserve capacity for reactive power flows.

``transmission projects``
=========================

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transmission_projects>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: transmission_projects:
   :end-before: pypsa_eur:

TYNDP 2020 dataset is not intergrated in this analysis.

``existing_capacities``
=======================

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#existing_capacities>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: existing_capacities:
   :end-before: # docs

The content is identical to the default file.

``sector``
=======================

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#sector>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: sector:
   :end-before: # docs

Configuration changes made:

* This analysis includes the transport and heating sectors, while excluding biomass, agriculture, and industry.
* The availability of battery electric vehicles (BEV) for demand-side management is reduced from 50% to 20%.
* The charging rate for BEVs is reduced from 1.1% to 0.7%.
* The predicted share of land electric vehicles is reduced from 45% to 30% by 2035.
* Land transport emissions are excluded from the model, as they are determined exogenously.
* The reduction in space heating demand is increased from 11% to 21% by 2035.
* Li-ion batteries, vanadium, liquid air, CAES, and Iron-air batteries are defined as storage unit components rather than stores.
* The hydrogen network is excluded from the model.

``load``
==============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#load>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: load:
   :end-before: # docs

The model's electricity load is based on data from 2013. 
By 2035, electricity demand is projected to increase by 26% compared to 2013. `Source from Agora Energiewende <https://www.agora-energiewende.de/fileadmin/Projekte/2021/2021_11_DE_KNStrom2035/2022-06-23_Praesentation_Klimaneutrales_Stromsystem_2035.pdf>`_

``costs``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#costs>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: costs:
   :end-before: # docs

Configuration changes made:

* The cost data from PyPSA technology cost database is based on the year 2035, instead of 2030.
* Gas price is set at 38.84 EUR/MWh. `Source from the EU Comission (Annex) <https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:52022SC0230>`_
* Coal price is set at 24.57 EUR/MWh `Source from Business Analytiq (2024 price) <https://businessanalytiq.com/procurementanalytics/index/subbituminous-coal-price-index/>`_
* Lignite price is set at 22.11 EUR/MWh `Source from Business Analytiq (2024 price) <https://businessanalytiq.com/procurementanalytics/index/lignite-coal-price-index/>`_
* Uranium price is set at 1.75 EUR/MWh `Source from Business Insider (2024 price) <https://markets.businessinsider.com/commodities/uranium-price>`_
* Iron-air battery price is set at 23,500 EUR/MWh.
* Storage energy cost for CAES is 30,000 EUR/MWh (2022).
* Storage power cost for CAES is 1,725,000 EUR/MW (2022).
* The lifetime of CAES is reduced to 40 years.

``clustering``
==============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#clustering>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: clustering:
   :end-before: # docs

Germany has a focus weight of 0.6, which means that out of 52 nodes, the distribution is as follows:

* Germany: 31 nodes
* France: 6 nodes
* Italy: 5 nodes
* Denmark: 2 nodes

The following countries each have 1 node, contributing to a total of 9 nodes:

* Austria
* Belgium
* Switzerland
* Czech Republic
* Luxembourg
* Netherlands
* Poland
* Sweden

The temporal resolution is reduced to 3-hour time steps.

``plotting``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#plotting>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: plotting:
   :end-before: # docs

Not much change here other than the coloring scheme

``solving``
=============

`Documentation <https://pypsa-eur.readthedocs.io/en/latest/configuration.html#solving>`_

.. literalinclude:: ../config/config.form.yaml
   :language: yaml
   :start-at: solving: