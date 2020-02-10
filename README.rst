=====
FEAST v3.0
=====

Introduction:
-------------
FEAST is a Fugitive Emissions Abatement Simulation Toolkit for evaluating natural gas leak detection and repair (LDAR) programs. FEAST can be used to estimate natural gas savings and net present value of LDAR programs. FEAST 1.0 was developed by the Environmental Assessment and Optimization group at Stanfor University and released in 2016. It was initially written in the Matlab programming language. Version 2.0 was written in Python and released in 2017. Version 3.0 incorporated additional emission data sets, eliminated explicit plume models in favor probability of detection curves to simulate detection technologies and incorporated permitted emissions that are not repaired by leak detection and repair programs.

Tutorial:
---------
Clicking the "launch/binder" link below will open the FEAST directory on a server that is able to run Python. Once the directory opens, select the file "Tutorial.ipynb" to open an interactive Jupyter notebook that illustrates what FEAST can do.

.. image:: http://mybinder.org/badge.svg 
	:target: http://mybinder.org:/repo/eaogroup/feast

Note: this tutorial is supported by mybinder.org. Occassionally, the mybinder service is unavailable. In that case the notebook will hang on the mybinder loading screen. You can check the status of mybinder here: http://mybinder.org/status If the service is down, you will have to either check back later or install python and clone the repository to your machine in order to run the tutorial. 

File Structure:
---------------
Python FEAST consists of a directory containing over 30 python module and object files. The file map below illustrates where each file is stored. The map is followed by a short description of the files.

::

	FEAST
	|----field_simulation.py
	|----Glossary.txt
	|----README.rst
    |----requirements.txt
    |----setup.cfg
    |----setup.py
	|----DetectionModules
		|----__init__.py
		|----abstract_detection_method.py
		|----helper_functions
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
			|----arpae_wind.p
            |----COGCC_site_prod_2019.p
			|----fernandez_leak_reapair_costs_2006.p
			|----fort_worth_leaks.p
            |----fort_worth_tank.p
            |----fort_worth_notank.p
			|----fort_worth_wind.p
			|----hitran_database.p
			|----pnnl_methane.p
            |----production-emissions.p
		|----RawData
            |----Allen_leakdata_2013.csv
            |----COGCC-2019-Production.xlsx
            |----COGCC-2019-well-locations.dbf
			|----ARPAEWind.csv
			|----FernandezRepairCost.csv
			|----FortWorth.csv
			|----FortWorthwindData.csv
			|----HITRAN_database.csv
            |----octave-leak-data.csv
			|----pnnl_methane.csv
            |----ProductionSite-ComponentEmissions.xlsx
		|----RawDataProcessingScripts
            |----allend_data_prep.py
			|----arpae_wind_reader.py
			|----COGCC_2019_prodctiondata.py
			|----fernandez_repair_cost_reader.py
			|----fort_worth_data_prep.py
			|----fort_worth_tank_notank.py
			|----fort_worth_wind_reader.py
			|----HITRAN_reader.py
			|----pnnl_reader.py
			|----production_emission_data.py
			|----README.txt
			|----repair_cost_data_reader.py
			|----wind_data_reader.py

File descriptions
-----------------
field_simulation.py 
	contains one function of the same name (field_simulation). One call to field_simulation() creates one realization of a FEAST 		scenario. field_simulation() accepts several optional input arguments to change parameters from their default settings.

DetectionModules:
-----------------
DetectionModules is the directory containing all of the LDAR program files:

abstract_detection_method.py 
	defines a parent class with the attributes and methods that all LDAR programs have. 

helper_functions.py 
	contains a few short functions that are used by multiple LDAR programs. 

null.py 
	defines the null detection method. 

dd.py 
	defines a "Distributed Detector" LDAR program

ir.py
	defines LDAR programs that use an infrared camera. An a manual subclass and an airborne, automated subclass are included.

fid.py
	defines an LDAR program based on a flame ionization detector

GeneralClassesFunctions:
------------------------
GeneralClassesFunctions contains files that define classes and functions that are not directly specified by LDAR programs or input data. Each module in the directory is described below:

leak_class_functions.py
	defines the Leak class used to store all the data required to define a set of leaks. The module also contains function
	definitions used to create and manipulate leak objects.

plotting_functions.py 
	defines functions for plotting simulation results.

results_analysis_functions.py 
	defines functions that compile results from numerous realizations of a scenario to calculate mean net present value, detected
	leak size distributions and other statistics. plotting_functions.py calls results_analysis_functions.py to produce plots.

simulation_classes.py 
	defines classes that are necessary for a simulation. These classes are GasField, FinanceSettings, Atmosphere, Time and Results.

simulation_functions.py 
	defines functions that are necessary for a simulation but are neither part of a LDAR program nor methods of a class. The
	functions are listed below:
	
	-new_leak_count      Calculates the number of new leaks to generate at a time step
	-save_results        Generates a Results object at the end of a simulation and saves it.
	-set_kwargs_attrs    Allows any attribute specified in a class to be set using key word arguments
	-gauss_leak_model    Calculates the concentration of gas due to a leak at specified location and conditions.


InputData:
----------
InputData is a directory containing raw data files, scripts for processing those raw data files and python object files created from the raw data. PyFEAST only uses the python object files, but the raw files and processing files are included for transparency and to allow for alternative processing files to be added in the future. The following list describes the subdirectories and class file in InputData.

input_data_classes.py    
	Defines all of the input data classes used by PyFEAST.
	
DataObjectInstances    
	Contains python data object files used by PyFEAST
	
RawData    
	Contains raw csv files for wind speed, leak data sets and other inputs to PyFEAST.
	
RawDataProcessingScripts    
	Contains the scripts used to produce the objects in DataObjectInstaces from the csv files in RawData.

Author:
-------
Chandler Kemp https://github.com/ChandlerKemp

Acknowledgments:
----------------
JP Addison reviewed all code developed for the Python implementation of FEAST.
