##########################################
Scenarios Configuration
##########################################

``Baseline``
=============

.. literalinclude:: ../config/scenarios.form.yaml
   :language: yaml
   :start-at: baseline-mds-1H:
   :end-at: baseline-nomds-1H:

The differences between this scenario and the following are:

* **baseline-nomds-1H**: It excludes both iron-air batteries and H2 storage.
* **baseline-mds-100H**: The temporal resolution is reduced to 100 hours per model.
* **baseline-mds-4380SEG**: The temporal resolution is reduced to 2-hour time steps.