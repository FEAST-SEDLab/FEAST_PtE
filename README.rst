Fugitive Emissions Abatement Simulation Tool (FEAST) v3.0
==========================================================

Creators: Chandler Kemp, Arvind Ravikumar
_________________________________________

Introduction: What is FEAST?
----------------------------
The Fugitive Emissions Abatement Simulation Toolkit or FEAST is a model to evaluate the effectiveness of methane leak detection and repair (LDAR) programs at oil and gas facilities. Recent advances in the development of new fixed and mobile (truck-, drone-, plane-, and satellite-based) methane leak detection technologies have led to growing interest in alternative LDAR programs. Thus, FEAST can also be used to compare the relative effectiveness of new technologies and methods as part of LDAR programs. FEAST uses publicly available data-sets on methane emissions and recent data from field tests of new methane leak detection technologies to simulate the performance of LDAR programs. 

A Brief History of FEAST
------------------------
FEAST was initially developed at the Environmental Assessment and Optimization group at Stanford University and released in 2016 (C.Kemp et al. Environ. Sci. Tech. 50 4546 http://dx.doi.org/10.1021/acs.est.5b06068). The latest version of FEAST, version 3.0, incorporates additional emission data sets, improves parametric representations of new technologies and methods, and explicitly accounts for leaks (unintentional) and vents (intentional) within the simulated oil and gas facility. 

FEAST-related Resources
------------------------
More details on FEAST: [Will Be Available Soon]
Webinar link: https://www.youtube.com/watch?v=RWi8CmmGOPA
FEAST publications:  [Will be Available Soon] 

Contact Information
-------------------
For any FEAST-related technical help, please send your questions to: feast.help [at] gmail [dot] com

Tutorial:
---------
Clicking the link below will direct you to a jupyter notebook server that demonstrates a few basic simulations in FEAST. The server will not support advanced simulations and your results will not be maintained indefinitely, so it should be used as an educational tool rather than an analysis platform. Enter any username and password to access the server. By recording your username and password you will be able to reaccess your work after the server times out.

http://feast.harrisburgu.cloud

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
	defines functions that are necessary for a simulation but are neither part of a LDAR program nor methods of a class. The functions are listed below::

	-save_results	        Generates and saves a Results object.
	-set_kwargs_attrs	Allows any attribute to be set using key word arguments.

site_emission_methdods.py
	allows initial gas field emisions to be controlled at the site level.

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

Results Files:
---------------
Results files can be loaded with the "pickle" package in Python's standard library as shown below:
    re = pickle.load(open('results_realization0.p', 'rb'))

The resulting object (re) contain the emissions timeseries from every LDAR program in the simulation, as well as all of
the settings used to run the simulation. The paths to most high level attributes are shown below:
    re.econ_settings----Defines economic variables, including the discount rate assumed in net present value analyses
    re.tech_dict----A dictionary of LDAR programs and their associated properties
        re.tech_dict['tech name'].emissions----list of total emissions at each time step (g/s)
        re.tech_dict['tech name'].find_cost----list of leak detection costs associated with each time step
        re.tech_dict['tech name'].repair_cost----list of repair costs associated with each time step
    re.time----Defines time variables, including simulation duration (endtime) and time step size (delta_t)
    re.gas_field----Stores all gas field settings and properties
        re.gas_field.input_leaks----List of emissions added at each time step.
