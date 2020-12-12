""" allen_data_prep loads data from the Allen emissions data spreadsheet and processes it before saving the
    results to a DataFile object. The Allen emissions data are from the following resource:
    D. Allen, V. M. Torres, et al. Measurements of methane emissions at natural
    gas production sites in the United States David". In: Proceedings of the National
    Academy of Sciences 110.44 (Sept. 2013), pp. 18025{18030. issn: 0027-8424.
    doi: 10.1073/pnas.1315099110. url: http://www.pnas.org/cgi/doi/10.
    1073/pnas.1315099110
"""

# -------------- reading the csv file --------------
import csv
from feast.input_data_classes import LeakData
import pickle
from os.path import dirname, abspath
import os

# -------------- Hard coded values --------------
# In some cases, a leak would be detected with an FID or IR camera, but no flux could be measured with the HI-FLOW
# sampler. In these cases, the study team assigned a flux of 0.001 cfm to the leak. These data are omitted from FEAST.
cfm_unmeasured_value = 0.001  # cfm
# Number of wells surveyed with an IR camera and FID in the Fort Worth study
n_wells_IR = 292
# Unit conversion from cfm to g/s (assuming standard conditions and pure methane)
cfm_to_gps = 0.0283/60*1e5/8.314/293*16

flux_IR = []  # g/s
counter = 0
flux = 0

rsc_path, _ = os.path.split(dirname(abspath(__file__)))
file_in = os.path.join(rsc_path, 'RawData', 'Allen_leakdata_2013.csv')
file_out = os.path.join(rsc_path, 'DataObjectInstances', 'allen_leaks.p')

with open(file_in) as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    source_name = ''
    for row in data:
        if row[8][0].isdigit():
            flux = float(row[8])*cfm_to_gps
            if flux > 0:
                flux_IR.append(float(row[8])*cfm_to_gps)

notes = \
    """Data extracted from D. Allen, V. M. Torres, et al. Measurements of methane emissions at natural
    gas production sites in the United States David". In: Proceedings of the National
    Academy of Sciences 110.44 (Sept. 2013), pp. 18025{18030. issn: 0027-8424.
    doi: 10.1073/pnas.1315099110. Flux data are recorded in grams/second."""

allen_leaks = LeakData(notes=notes, raw_file_name='Allen_leakdata_2013.csv', data_prep_file='allen_data_prep.py')

leak_data = {'IR': flux_IR}
well_counts = {'IR': n_wells_IR}
comp_counts = {'IR': None}

allen_leaks.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)

pickle.dump(allen_leaks, open(file_out, 'wb'))
