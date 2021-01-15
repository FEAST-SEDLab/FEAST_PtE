Welcome to FEAST's documentation!
========================================

 .. toctree::
    :maxdepth: 5
    :caption: Contents:
 
FEAST modules
========================================

DetectionModules
----------------------------------------
abstract_detection_method
****************************************
.. automodule:: feast.DetectionModules.abstract_detection_method

DetectionMethod
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.abstract_detection_method.DetectionMethod
    :members:

comp_survey
****************************************
.. automodule:: feast.DetectionModules.comp_survey

CompSurvey
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.comp_survey.CompSurvey
    :members:

ldar_program
****************************************
.. automodule:: feast.DetectionModules.ldar_program

LDARProgram
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.ldar_program.LDARProgram
    :members:

repair
****************************************
.. automodule:: feast.DetectionModules.repair

Repair
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.repair.Repair
    :members:

site_monitor
****************************************
.. automodule:: feast.DetectionModules.site_monitor

SiteMonitor
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.site_monitor.SiteMonitor
    :members:

site_survey
****************************************
.. automodule:: feast.DetectionModules.site_survey

SiteSurvey
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.site_survey.SiteSurvey
    :members:

EmissionSimModules
----------------------------------------
emission_class_functions
****************************************
.. automodule:: feast.EmissionSimModules.emission_class_functions

Emission
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.emission_class_functions.Emission
    :members:

bootstrap_emission_maker
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.EmissionSimModules.emission_class_functions.bootstrap_emission_maker

comp_indexes_fcn
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.EmissionSimModules.emission_class_functions.comp_indexes_fcn

emission_objects_generator
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.EmissionSimModules.emission_class_functions.emission_objects_generator

permitted_emission
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.EmissionSimModules.emission_class_functions.permitted_emission

infrastructure_classes
****************************************
.. automodule:: feast.EmissionSimModules.infrastructure_classes

Component
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.infrastructure_classes.Component
    :members:

GasField
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.infrastructure_classes.GasField
    :members:

Site
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.infrastructure_classes.Site
    :members:

result_classes
****************************************
.. automodule:: feast.EmissionSimModules.result_classes

ResultAggregate
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultAggregate
    :members:

ResultContinuous
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultContinuous
    :members:

ResultDiscrete
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultDiscrete
    :members:

simulation_classes
****************************************
.. automodule:: feast.EmissionSimModules.simulation_classes

Scenario
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.simulation_classes.Scenario
    :members:

Time
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.simulation_classes.Time
    :members:

GeneralClassesFunctions
----------------------------------------
InputData
----------------------------------------
RawDataProcessingScripts
****************************************
input_data_classes
----------------------------------------
.. automodule:: feast.input_data_classes

DataFile
****************************************
.. autoclass:: feast.input_data_classes.DataFile
    :members:

LeakData
****************************************
.. autoclass:: feast.input_data_classes.LeakData
    :members:

ProductionData
****************************************
.. autoclass:: feast.input_data_classes.ProductionData
    :members:

RepairData
****************************************
.. autoclass:: feast.input_data_classes.RepairData
    :members:

ResultsProcessing
----------------------------------------
plotting_functions
****************************************
.. automodule:: feast.ResultsProcessing.plotting_functions

abatement_cost_plotter
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.plotting_functions.abatement_cost_plotter

plot_fixer
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.plotting_functions.plot_fixer

time_series
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.plotting_functions.time_series

results_analysis_functions
****************************************
.. automodule:: feast.ResultsProcessing.results_analysis_functions

npv_calculator
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.results_analysis_functions.npv_calculator

results_analysis
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.results_analysis_functions.results_analysis

