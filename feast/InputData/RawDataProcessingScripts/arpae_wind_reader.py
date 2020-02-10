from .wind_data_reader import wind_data_reader
from os.path import dirname, abspath

rsc_path = '/'.join(dirname(abspath(__file__)).split('/')[:-1])
file_in = '/'.join([rsc_path, 'RawData', 'ARPAEWind.csv'])
file_out = '/'.join([rsc_path, 'DataObjectInstances/arpae_wind.p'])
Notes = 'Data are the US Department of Energy Methane Observation Networks with Innovative Technology to Obtain' +\
        'Reductions – MONITOR wind speed data, 2014.'
speed_col = 5
time_col = 2

wind_data_reader(file_in, speed_col=5, direction_col=None, time_col=2, notes=Notes, file_out=file_out)
