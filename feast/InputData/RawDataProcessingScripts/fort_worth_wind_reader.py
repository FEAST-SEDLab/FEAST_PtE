from .wind_data_reader import wind_data_reader
from os.path import dirname, abspath
import os

rsc_path, _ = os.path.split(dirname(abspath(__file__)))
file_in = os.path.join(rsc_path, 'RawData', 'FortWorthWindData.csv')
file_out = os.path.join(rsc_path, 'DataObjectInstances', 'fort_worth_wind.p')

Notes = 'Data are from NOAA National Weather Service Forecast Office Dallas/FortWorth, TX Climate Data. National' + \
            'Climatic Data Center. 2014'
speed_col = 4
direction_col = 3
time_col = 6

wind_data_reader(file_in, speed_col, direction_col, time_col, Notes, file_out=file_out)
