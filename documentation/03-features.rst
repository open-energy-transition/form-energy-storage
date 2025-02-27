##########################################
Features
##########################################

Here is a list of changes made to this repository specifically for this work, with the potential to be upstreamed to the main PyPSA-Eur repository.

**Model Features for Energy Storage Technologies**

* Add structure in ``prepare_sector_network`` to specify a modelling approach of storages as either ``Stores`` or ``StorageUnits`` (https://github.com/open-energy-transition/pypsa-eur/pull/15)

* Add iron-air battery as storage technology with implementation as ``StorageUnit`` or ``Store`` (https://github.com/open-energy-transition/pypsa-eur/pull/20)

* Add Lithium Iron Phosphate, Vanadium, Liquid-air, Compressed-air energy storage technologies as ``StorageUnit`` (https://github.com/open-energy-transition/pypsa-eur/pull/21)

**Model Adjustment Features**

* Add an option to overwrite cost table attributes with the config in ``prepare_costs`` and ``load_costs`` (https://github.com/open-energy-transition/pypsa-eur/pull/23)

* Add ramping limit options for conventional powerplants in ``prepare_sector_network``. (https://github.com/open-energy-transition/form-energy-storage/pull/28)

* Adjust emission to only include the sector scope (https://github.com/open-energy-transition/form-energy-storage/pull/29)

* Add option for ``final_adjustment`` before solving the network. This script is used to limit the electricity grid interconnecting capacity for each country (https://github.com/open-energy-transition/form-energy-storage/pull/33)

* Add a constraint to limit the use of Direct Air Capture (DAC) (https://github.com/open-energy-transition/form-energy-storage/pull/42)

* Include heat DSM implementation and configuration (https://github.com/open-energy-transition/form-energy-storage/pull/47)

* Add option to disable CHP capacity extendability for Germany in ``final_adjustment`` (https://github.com/open-energy-transition/form-energy-storage/pull/62)

* Add a ``filter_year`` configuration option that is used to filter out projects that have a build_year until and including the given filter_year (https://github.com/open-energy-transition/form-energy-storage/pull/65)

**Model Visualization Features**

* Add KPI visualization script and functions. (https://github.com/open-energy-transition/form-energy-storage/pull/30)

* Add new battery techs to nice names and plotting for cost map (https://github.com/open-energy-transition/form-energy-storage/pull/31)

* Add plotting script for FE KPIs (https://github.com/open-energy-transition/form-energy-storage/pull/37)

* Include optional use of Latex for plotting (https://github.com/open-energy-transition/form-energy-storage/pull/59)

* Create notebook with scenarios comparison. All plots and PyPSA networks can be exported to csvs (https://github.com/open-energy-transition/form-energy-storage/pull/61)

* Add csvs to all map plots and fixes in line loading plots (https://github.com/open-energy-transition/form-energy-storage/pull/68)

**Model Documentation Features**

* Provide documentation using github pages (https://github.com/open-energy-transition/form-energy-storage/pull/57)

**Model Calibration Runs**

* Calibration run 2023/2024 (https://github.com/open-energy-transition/form-energy-storage/pull/45)

**Bugfixes**

* Fix deprecation warning inside ``prepare_cost`` and ``load_cost`` functions (https://github.com/open-energy-transition/pypsa-eur/pull/24)

* Fix bugs of adding nuclear capacities twice (https://github.com/open-energy-transition/form-energy-storage/pull/32)

* Fix for ``heat_dsm_profile`` for leap year (https://github.com/open-energy-transition/form-energy-storage/pull/60)
