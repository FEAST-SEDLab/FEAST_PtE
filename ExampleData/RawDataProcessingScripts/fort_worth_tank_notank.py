""" fortWorthDataPrep loads data from the Fort Worth emissions data spreadsheet and processes it before saving the
    results to a DataFile object. The Fort Worth emissions data are available from the following resource:
    City of Fort Worth Natural Gas Air Quality Study. Tech. rep. Fort Worth, TX:
    Eastern Research Group and Sage Environmental Consulting LP. for the City
    of Fort Worth, 2011. <http://fortworthtexas.gov/gaswells/air-quality-study/final> (Emissions Calculations Workbook)
    For convenience, the sheet titled '2.2 Emission Data' was exported from the original document, commas were stripped
    and the file is saved as a csv in the DataFiles directory of PyFEAST.

    Column index 28 (AC in Excel) is used to define flux. This column is labeled "%CFM." That label can be confusing, so
    more detail is provided here. At every emission point in the data set, a high flow sampler was used to measure the
    emission rate. The high flow sampler records the total flow through the sampler (CFM), and the "the leak rate as %
    of total sample flow (% CFM)" (page 3-17 in the Fort Worth Natural Gas Air Quality Study Final Report). The
    definition of % CFM provided in the report seems to say that the percentage of the total gas flow that is methane is
    recorded as % CFM. However, page 3-18 implies that % CFM is equivalent to the actual emission rate: "The % CFM
    reading obtained with the Hi Flow Sampler from that point exceeded the daily rolling average % CFM for all Hi Flow
    Sampler tests conducted thus far. In other words, the emission rate had to equal or exceed the average emission
    rate." Furthermore subsequent correlation equations use %CFM as a the HiFlow Sampler's measure of emissions on page
    3-26. Finally, within the spreadsheet itself, %CFM is calculated as the product of the total flow through the
    sampler and the methane concentration. Therefore, %CFM is a measure of the total leak rate by the high flow sampler.
"""

# -------------- reading the csv file --------------
import csv
from feast.input_data_classes import LeakData
import pickle
from os.path import dirname, abspath
import os

rsc_path, _ = os.path.split(dirname(abspath(__file__)))
file_in = os.path.join(rsc_path, 'RawData', 'FortWorth.csv')
file_out_tank = os.path.join(rsc_path, 'DataObjectInstances', 'fort_worth_tank.p')
file_out_notank = os.path.join(rsc_path, 'DataObjectInstances', 'fort_worth_notank.p')


# -------------- Hard coded values --------------
# In some cases, a leak would be detected with an FID or IR camera, but no flux could be measured with the HI-FLOW
# sampler. In these cases, the study team assigned a flux of 0.001 cfm to the leak. These data are omitted from FEAST.
cfm_unmeasured_value = 0.001  # cfm
# Number of wells surveyed with an IR camera and FID in the Fort Worth study
n_wells_IR, n_wells_FID = 1138, 114
# Number of components per well surveyed with an IR camera and FID in the Fort Worth Study (See EmissionsData.xslx,
# sheet "Executive PS Site Summary")
n_comp_IR, n_comp_FID = round(736659 / 1138), round(736659 / 1138)

# Unit conversion from cfm to g/s (assuming standard conditions and pure methane)
cfm_to_gps = 0.0283/60*1e5/8.314/293*16

flux_IR, flux_FID = [], []  # g/s
# There are some characters in the csv file that raise errors. These characters are in text and are ignored by the
# errors='ignore' flag
tank = {
    'IR': [],
    'FID': []
}

notank = {
    'IR': [],
    'FID': []
}

with open(file_in, encoding='utf-8', errors='ignore') as csvfile:
    data = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in data:
        if row[0].isdigit():
            flux = float(row[28])
            if flux != cfm_unmeasured_value:
                flux *= cfm_to_gps
                if row[7] == 'TK':
                    if row[29] == '--':
                        tank['FID'].append(flux)
                    else:
                        tank['IR'].append(flux)
                else:
                    if row[29] == '--':
                        notank['FID'].append(flux)
                    else:
                        notank['IR'].append(flux)

note_base = \
    """Data extracted from City of Fort Worth Natural Gas Air Quality Study. Tech. rep. Fort Worth, TX:
    Eastern Research Group and Sage Environmental Consulting LP. for the City
    of Fort Worth, 2011. Flux data are recorded in grams/second."""

note_tank = note_base + " Only data labeled TK are included."
note_notank = note_base + " Only data labeled NTK are included."
fname = 'fort_worth_tank_notank.py'
fort_worth_tank = LeakData(notes=note_tank, raw_file_name=file_in.split('/')[-1], data_prep_file=fname)
fort_worth_notank = LeakData(notes=note_notank, raw_file_name=file_in.split('/')[-1], data_prep_file=fname)

well_counts = {'IR': n_wells_IR, 'FID': n_wells_FID}
comp_counts = {'IR': n_comp_IR, 'FID': n_comp_FID}

fort_worth_tank.define_data(leak_data=tank, well_counts=well_counts, comp_counts=comp_counts)
fort_worth_notank.define_data(leak_data=notank, well_counts=well_counts, comp_counts=comp_counts)

pickle.dump(fort_worth_tank, open(file_out_tank, 'wb'))
pickle.dump(fort_worth_notank, open(file_out_notank, 'wb'))
