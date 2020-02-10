"""
octave_matching_input loads data employed in Ravikumar's matlab code to a DataFile object so that the new program can be
compared
"""

# -------------- reading the csv file --------------
import csv
from InputData.input_data_classes import LeakData
import pickle
# -------------- Hard coded values --------------
# In some cases, a leak would be detected with an FID or IR camera, but no flux could be measured with the HI-FLOW
# sampler. In these cases, the study team assigned a flux of 0.001 cfm to the leak. These data are omitted from FEAST.
cfm_unmeasured_value = 0.001  # cfm
# Number of wells surveyed (set to equal 8 leaks per well)
n_wells = 344

flux_array = []  # g/s
counter = 0
flux = 0
with open('InputData/RawData/octave-leak-data.csv') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    source_name = ''
    for row in data:
        if row[0][0].isdigit():
            flux = float(row[0])
            if flux > 0:
                flux_array.append(flux)

notes = \
    """These data are for program testing purposes only. They are a conglomeration of several datasets"""

octave_leaks = LeakData(notes=notes, raw_file_name='octave-leak-data.csv', data_prep_file='octave-matching-input.py')

leak_data = {'all': flux_array}
well_counts = {'all': n_wells}
octave_leaks.define_data(leak_data=leak_data, well_counts=well_counts)

pickle.dump(octave_leaks, open('InputData/DataObjectInstances/octave_leaks.p', 'wb'))
