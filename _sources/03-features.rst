##########################################
Features
##########################################

Here is a list of changes made to this repository specifically for this work, with the potential to be upstreamed to the main PyPSA-Eur repository.

**Model Features for Energy Storage Technologies**

* Add structure in ``prepare_sector_network`` to specify a modelling approaches of storage as either ``Stores`` or ``StorageUnits`` (https://github.com/open-energy-transition/pypsa-eur/pull/15)

* Add iron-air battery as storage technology with implementation as ``StorageUnit`` or ``Store`` (https://github.com/open-energy-transition/pypsa-eur/pull/20)

* Add Lithium Iron Phosphate, Vanadium, Liquid-air, Compressed-air energy storage technologies as ``StorageUnit`` (https://github.com/open-energy-transition/pypsa-eur/pull/21)

**Model Adjustment Features**

* Add an option to overwrite cost table attributes with the config in ``prepare_costs`` and ``load_costs`` (https://github.com/open-energy-transition/pypsa-eur/pull/23)

* Add ramping limit options for conventional powerplants in ``prepare_sector_network``. (https://github.com/open-energy-transition/form-energy-storage/pull/28)

* Adjust emission to only include the sector scope (https://github.com/open-energy-transition/form-energy-storage/pull/29)

* Add option for ``final_adjustment`` before solving the network. This script is used to limit the grid capacity for each country (https://github.com/open-energy-transition/form-energy-storage/pull/33)

* Add a constraint to limit the use of Direct Air Capture (DAC) (https://github.com/open-energy-transition/form-energy-storage/pull/42)

* Merge heat DSM implementation and configuration (https://github.com/open-energy-transition/form-energy-storage/pull/47)

* Add solid biomass and biogas sources for bioenergy power plants (https://github.com/open-energy-transition/form-energy-storage/pull/54)


**Model Visualization Features**

* Add KPI visualization script and functions. (https://github.com/open-energy-transition/form-energy-storage/pull/30)

* Add new battery techs to nice names and plotting for cost map (https://github.com/open-energy-transition/form-energy-storage/pull/31)

* Add plotting script for FE KPIs (https://github.com/open-energy-transition/form-energy-storage/pull/37)

* Make the use of Latex for plotting optional (https://github.com/open-energy-transition/form-energy-storage/pull/59)

**Model Calibration Runs**

* Calibration run 2023/2024 (https://github.com/open-energy-transition/form-energy-storage/pull/45)

* Update calibrated run branch with the recent addition in ``form_energy_storage_dev`` (https://github.com/open-energy-transition/form-energy-storage/pull/48)

**Bugfixes**

* Fix deprecation warning inside ``prepare_cost`` and ``load_cost`` functions (https://github.com/open-energy-transition/pypsa-eur/pull/24)

* Fix bugs of adding nuclear capacities twice (https://github.com/open-energy-transition/form-energy-storage/pull/32)

* Fix for ``heat_dsm_profile`` for leap year (https://github.com/open-energy-transition/form-energy-storage/pull/60)
