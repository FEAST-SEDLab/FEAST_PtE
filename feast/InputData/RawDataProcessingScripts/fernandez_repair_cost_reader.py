from .repair_cost_data_reader import repair_cost_data_reader
from os.path import dirname, abspath

rsc_path = '/'.join(dirname(abspath(__file__)).split('/')[:-1])
file_in = '/'.join([rsc_path, 'RawData', 'FernandezRepairCost.csv'])
file_out = '/'.join([rsc_path, 'DataObjectInstances/fernandez_leak_repair_costs_2006.p'])

Notes = 'Data are from Fernandez, R. "Cost Effective Directed Inspection and Maintenance..." National Gas Machinery' + \
        'Laboratory. 2006'

repair_cost_data_reader(file_in, Notes, file_out)
