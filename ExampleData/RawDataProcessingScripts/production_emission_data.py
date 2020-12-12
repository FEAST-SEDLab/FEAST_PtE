"""
    production_emission_data loads emissions data from numerous component level studies that are compiled in the
    ProductionSite-ComponentEmissions.xlsx spreadsheet.
"""

# -------------- reading the spreadsheet --------------
import pandas as pd
from feast.input_data_classes import LeakData
import pickle
import numpy as np
from os.path import dirname, abspath
import os

rsc_path = dirname(abspath(__file__))
rsc_path, _ = os.path.split(rsc_path)
file_in = os.path.join(rsc_path, 'RawData', 'ProductionSite-ComponentEmissions.xlsx')
file_out = os.path.join(rsc_path, 'DataObjectInstances', 'production_emissions.p')

dat = pd.read_excel(file_in,
                    sheet_name="AllSources-kgday-methane",
                    header=3)
em_array = []
for study in dat:
    if 'Ravikumar-Measured' in study:
        continue
    em_array.extend(dat[study])
em_array = np.array(em_array)
em_array = em_array[np.invert(np.isnan(em_array))]
em_array = em_array * 1000 / 24 / 3600  # convert from kg/day to g/s
em_array = em_array[em_array > 0]

notes = \
    """
    Data extracted from the compilation spreadsheet ProductionSite-ComponentEmissions.xlsx
    The number of components surveyed at each well generally were not recorded. Therefore, the number of components is 
    estimated by assuming 650 components per well.    
    """

emissions = LeakData(notes=notes, raw_file_name=file_in.split('/')[-1], data_prep_file='production_emission_data.py')

# The dict structure allows for multiple types of detection methods used in the study
leak_data = {'All': em_array}
well_counts = {'All': 2612}  # Number of wells in the study
comp_counts = {'All': 650}  # Assumed components per well

emissions.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)

pickle.dump(emissions, open(file_out, 'wb'))

print('Successfully completed production-emission-data processing.')
