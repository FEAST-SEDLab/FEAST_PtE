"""
    This script can be used to generate a data object that for FEAST to use as an emission distribution dataset.
"""

# -------------- reading the spreadsheet --------------
import pandas as pd
from feast.input_data_classes import LeakData
import pickle
import os
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]
print(file_in)
# def gen_dat_file(file_in, file_out):
"""
Loads a raw data file following the raw_data_basic_template.xlsx format and creates a LeakData object instance.
:param file_in: Path to a raw data excel file
:param file_out: Path at which to save the output file
:return: None
"""
dat = pd.read_excel(
    file_in,
    header=1
)

notes = dat['Data'][0]

emissions = LeakData(notes=notes, raw_file_name=os.path.split(file_in)[-1], data_prep_file='template_processor.py')

# The dict structure allows for multiple types of detection methods used in the study
leak_data = {'All': dat['Data'][4:]}
well_counts = {'All': dat['Data'][1]}  # Number of wells in the study
comp_counts = {'All': dat['Data'][2]}  # Assumed components per well

emissions.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)

pickle.dump(emissions, open(file_out, 'wb'))

print('Successfully created emission size distribution object.')
