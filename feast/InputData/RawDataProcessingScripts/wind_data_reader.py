def wind_data_reader(file_in=None, speed_col=None, direction_col=None, time_col=None, notes=None,
                     file_out=None):
    """
    windDataReader reads a csv file containing wind data and returns a WindData object. It also saves the WindData
    object to file_out.
    Inputs:
        file_path_in     Path to a data file
        wind_speed       column number for wind speed data (should be provided in m/s)
        wind_direction   column number for wind direction data (should be provided in degrees East of North)
        time_col         column number of time data (should be provided in units of days)
        return           wind data object
    """
    import csv
    from ..input_data_classes import WindData
    import pickle

    wind_dict = dict()
    wind_dict[speed_col] = [] if speed_col else None
    wind_dict[time_col] = [] if time_col else None
    wind_dict[direction_col] = [] if direction_col else None

    with open(file_in) as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            if row[0][0].isdigit():
                for key in wind_dict.keys():
                    if key:
                        wind_dict[key].append(float(row[key-1]))

    wind_data_out = WindData(raw_file_name=file_in, notes=notes)
    wind_data_out.define_data(wind_speed=wind_dict[speed_col], wind_direction=wind_dict[direction_col],
                              time=wind_dict[time_col])

    if file_out is not None:
        pickle.dump(wind_data_out, open(file_out, 'wb'))
    return wind_data_out
