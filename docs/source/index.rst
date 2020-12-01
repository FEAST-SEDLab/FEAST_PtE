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

comp_survey
****************************************
.. automodule:: feast.DetectionModules.comp_survey

CompSurvey
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.comp_survey.CompSurvey

ldar_program
****************************************
.. automodule:: feast.DetectionModules.ldar_program

LDARProgram
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.ldar_program.LDARProgram

repair
****************************************
.. automodule:: feast.DetectionModules.repair

Repair
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.repair.Repair

site_monitor
****************************************
.. automodule:: feast.DetectionModules.site_monitor

SiteMonitor
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.site_monitor.SiteMonitor

site_survey
****************************************
.. automodule:: feast.DetectionModules.site_survey

SiteSurvey
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.DetectionModules.site_survey.SiteSurvey

EmissionSimModules
----------------------------------------
emission_class_functions
****************************************
.. automodule:: feast.EmissionSimModules.emission_class_functions

Emission
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.emission_class_functions.Emission

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

GasField
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.infrastructure_classes.GasField

Site
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.infrastructure_classes.Site

result_classes
****************************************
.. automodule:: feast.EmissionSimModules.result_classes

ResultAggregate
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultAggregate

ResultContinuous
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultContinuous

ResultDiscrete
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.result_classes.ResultDiscrete

simulation_classes
****************************************
.. automodule:: feast.EmissionSimModules.simulation_classes

Scenario
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.simulation_classes.Scenario

Time
++++++++++++++++++++++++++++++++++++++++
.. autoclass:: feast.EmissionSimModules.simulation_classes.Time

input_data_classes
----------------------------------------
.. automodule:: feast.input_data_classes

DataFile
****************************************
.. autoclass:: feast.input_data_classes.DataFile

LeakData
****************************************
.. autoclass:: feast.input_data_classes.LeakData

ProductionData
****************************************
.. autoclass:: feast.input_data_classes.ProductionData

RepairData
****************************************
.. autoclass:: feast.input_data_classes.RepairData

MEET_1_importer
----------------------------------------
.. automodule:: feast.MEET_1_importer

gas_comp_to_dict
****************************************
.. autofunction:: feast.MEET_1_importer.gas_comp_to_dict

gascomp_reader
****************************************
.. autofunction:: feast.MEET_1_importer.gascomp_reader

gc_dat_to_gas_field
****************************************
.. autofunction:: feast.MEET_1_importer.gc_dat_to_gas_field

load_gas_comp_file
****************************************
.. autofunction:: feast.MEET_1_importer.load_gas_comp_file

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

summary_plotter
++++++++++++++++++++++++++++++++++++++++
.. autofunction:: feast.ResultsProcessing.plotting_functions.summary_plotter

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

