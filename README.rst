=====
FEAST v3.0
=====

Introduction:
-------------
FEAST is a Fugitive Emissions Abatement Simulation Toolkit for evaluating natural gas leak detection and repair (LDAR) programs. FEAST can be used to estimate natural gas savings and net present value of LDAR programs. FEAST 1.0 was developed by the Environmental Assessment and Optimization group at Stanfor University and released in 2016. It was initially written in the Matlab programming language. Version 2.0 was written in Python and released in 2017. Version 3.0 incorporated additional emission data sets, eliminated explicit plume models in favor probability of detection curves to simulate detection technologies and incorporated permitted emissions that are not repaired by leak detection and repair programs.

Tutorial:
---------
Clicking the link below will direct you to a jupyter notebook server that demonstrates a few basic simulations in FEAST. The server will not support advanced simulations and your results will not be maintained indefinitely, so it should be used as an educational tool rather than an analysis platform. Enter any username and password to access the server. By recording your username and password you will be able to reaccess your work after the server times out.

http://64.225.7.26/

File Structure:
---------------
Python FEAST consists of a directory containing over 30 python module and object files. The file map below illustrates where each file is stored. The map is followed by a short description of the files.

::
	
	FEAST_PtE
	|----field_simulation.py
	|----Glossary.txt
	|----README.rst
	|----requirements.txt
	|----setup.cfg
	|----setup.py
	|----FEAST
		|----DetectionModules
			|----__init__.py
			|----abstract_detection_method.py
			|----null.py
			|----tech_detect.py
			|----tiered_detect.py
		|----GeneralClassesFunctions
			|----__init__.py
			|----leak_class_functions.py
			|----plotting_functions.py
			|----results_analysis_functions.py
			|----simulation_classes.py
			|----simulation_functions.py
			|----site_emission_methods.py
		|----InputData
			|----input_data_classes.py
			|----DataObjectInstances
				|----allen_leaks.p
				|----COGCC_site_prod_2019.p
				|----fernandez_leak_reapair_costs_2006.p
				|----fort_worth_leaks.p
				|----fort_worth_notank.p
				|----fort_worth_tank.p
				|----production-emissions.p
			|----RawData
				|----Allen_leakdata_2013.csv
				|----COGCC-2019-Production.xlsx
				|----COGCC-2019-well-locations.dbf
				|----FernandezRepairCost.csv
				|----FortWorth.csv
				|----ProductionSite-ComponentEmissions.xlsx
			|----RawDataProcessingScripts
				|----allend_data_prep.py
				|----COGCC_2019_prodctiondata.py
				|----fernandez_repair_cost_reader.py
				|----fort_worth_data_prep.py
				|----fort_worth_tank_notank.py
				|----production_emission_data.py
				|----README.txt
				|----repair_cost_data_reader.py


File descriptions
-----------------
field_simulation.py 
	contains one function of the same name (field_simulation). One call to field_simulation() creates one realization of a FEAST scenario. field_simulation() accepts several optional input arguments to change parameters from their default settings.

DetectionModules:
=================
DetectionModules is the directory containing all of the LDAR program files:

abstract_detection_method.py 
	defines a parent class with the attributes and methods that all LDAR programs have. 

null.py 
	defines the null detection method. 
	
tech_detect.py
	defines a class of LDAR programs that perform periodic surveys and identify emission sources at the component level.

tiered_detect.py
	defines a class of LDAR programs that perform periodic surveys that identify emission sources at the site level. Secondary surveys are used to identify the emitting components at high emitting sites.

GeneralClassesFunctions:
------------------------
GeneralClassesFunctions contains files that define classes and functions that store simulation settings and gas field states:

leak_class_functions.py
	defines the Leak class used to store all the data required to define a set of leaks. The module also contains function definitions used to create and manipulate leak objects.

plotting_functions.py 
	defines functions for plotting simulation results.

results_analysis_functions.py 
	defines functions that compile results from numerous realizations of a scenario to calculate mean net present value, detected	leak size distributions and other statistics. plotting_functions.py calls results_analysis_functions.py to produce plots.

simulation_classes.py 
	defines classes that are necessary for a simulation. These classes are Component, GasField, FinanceSettings, Results, Site and Time.

simulation_functions.py 
	defines functions that are necessary for a simulation but are neither part of a LDAR program nor methods of a class. The functions are listed below:
	
	-save_results	        Generates a Results object at the end of a simulation and saves it.	
	-set_kwargs_attrs	Allows any attribute to be set using key word arguments.

InputData:
----------
InputData is a directory containing raw data files, scripts for processing those raw data files and python object files created from the raw data. FEAST v3.0 only uses the python object files, but the raw files and processing files are included for transparency and to allow for alternative processing files to be added in the future. The following list describes the subdirectories and class file in InputData.

input_data_classes.py    
	Defines all of the input data classes used by FEAST.
	
DataObjectInstances    
	Contains python data object files used by FEAST
	
RawData    
	Contains raw csv files for leak data sets and other inputs to PyFEAST.
	
RawDataProcessingScripts    
	Contains the scripts used to produce the objects in DataObjectInstaces from the files in RawData.

Author:
-------
Chandler Kemp https://github.com/ChandlerKemp

Acknowledgments:
----------------
JP Addison reviewed all code developed for the Python implementation of FEAST.
