import csv
import numpy as np
from ..input_data_classes import HITRAN
import pickle
from os.path import dirname, abspath

rsc_path = '/'.join(dirname(abspath(__file__)).split('/')[:-1])
file_in = '/'.join([rsc_path, 'RawData', 'HITRAN_database.csv'])
file_out = '/'.join([rsc_path, 'DataObjectInstances/hitran_database.p'])

def hitran_reader(file_in=None, notes=None, file_out=None):
    """
    reads a csv file containing pnnl methane data.
    file_in:      Path to a raw data file
    notes:        Comments on the data
    file_out:     Path to store the data at
    """

    file = open(file_in)
    data = csv.reader(file, delimiter=',')
    n_lines = sum(1 for row in data)
    file.seek(0)
    index = 0
    nu, s, gamma, temp = np.zeros(n_lines), np.zeros(n_lines), np.zeros(n_lines), np.zeros(n_lines)
    delta, lower_e = np.zeros(n_lines), np.zeros(n_lines)
    for row in data:
        if row[0][0].isdigit():
            nu[index] = float(row[0])
            s[index] = float(row[1])
            gamma[index] = float(row[2])
            temp[index] = float(row[3])
            delta[index] = float(row[4])
            lower_e[index] = float(row[5])
            index += 1
    nu += delta
    if notes is None:
        notes = "Data from the HITRAN 2012 molecular spectroscopic database"
    hitran = HITRAN(notes, 'pnnl_methane.csv')
    hitran.define_data(nu, s, gamma, temp, lower_e)
    pickle.dump(hitran, open(file_out, 'wb'))


hitran_reader(file_in=file_in, notes=None, file_out=file_out)