from .repair_cost_data_reader import repair_cost_data_reader
from os.path import dirname, abspath
import os

rsc_path = dirname(abspath(__file__))
rsc_path, _ = os.path.split(rsc_path)
file_in = os.path.join(rsc_path, 'RawData', 'FernandezRepairCost.csv')
file_out = os.path.join(rsc_path, 'DataObjectInstances', 'fernandez_leak_repair_costs_2006.p')

Notes = 'Data are from Fernandez, R. "Cost Effective Directed Inspection and Maintenance..." National Gas Machinery' + \
        'Laboratory. 2006'

repair_cost_data_reader(file_in, Notes, file_out)
