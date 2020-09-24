"""
This module defines all classes used to store input data
"""
from feast.EmissionSimModules.simulation_functions import set_kwargs_attrs


class DataFile:
    """
    DataFile is an abstract super class used to store all data files that may be called in FEAST.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None, **kwargs):
        """
        notes           A string containing notes on the object created
        raw_file_name   path (or list of paths) to a raw input file(s)
        kwargs          optional input dictionary that will override default parameters or add new parameters
        """
        self.notes = notes
        self.raw_file_name = raw_file_name
        set_kwargs_attrs(self, kwargs, only_existing=False)


class LeakData(DataFile):
    """
       LeakData is designed to store all leak size data from a reference. It accommodates multiple detection methods
       within a single instance
    """
    def __init__(self, notes='No notes provided', raw_file_name=None, data_prep_file=None, leak_sizes=None):
        """
        notes           A string containing notes on the object created
        raw_file_name   path to a raw input file
        leak_sizes      list of leak sizes. If leaks were detected using multiple methods, leak_sizes must be a dict
                        with one key for each detection method
        """
        DataFile.__init__(self, notes, raw_file_name, leak_sizes=leak_sizes, data_prep_file=data_prep_file)
        self.leak_sizes = dict()
        self.well_counts = dict()
        self.comp_counts = dict()

    def define_data(self, leak_data=None, well_counts=None, comp_counts=None, keys=None):
        """
        Check data formatting and set keys
        leak_data       leak_data must be a dict if there are multiple detection methods. If there is exactly one
                        detection method, leak_data may be a list.
        well_counts     well_counts lists the number of wells monitored by each detection method in keys
        keys            each detection method should have a unique key associated with it
        """
        if keys is None:
            # leak_data should be a dict if there are multiple detection methods. Each key will define a detection
            # method and the associated dictionary entry will be the list of fluxes from leaks found using that
            # detection method
            if type(leak_data) is dict:
                keys = leak_data.keys()
            else:
                keys = ['all_leaks']
        if keys != leak_data.keys() or leak_data.keys() != well_counts.keys():
            error_str = "The 'keys' argument passed to LeakData.defineData(), leak_data.keys() and " \
                        "well_counts.keys() must all be equivalent"
            raise ValueError(error_str)

        for key in keys:
            self.leak_sizes[key] = leak_data[key]
            self.well_counts[key] = well_counts[key]
            self.comp_counts[key] = comp_counts[key]


class WindData(DataFile):
    """
    WindData is designed to store all wind data from a reference. It accommodates wind speed, direction and time.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None):
        """
        notes           A string containing notes on the object created
        raw_file_name   path to a raw input file
        """
        DataFile.__init__(self, notes, raw_file_name)
        self.wind_speed, self.time, self.wind_direction = [], [], []  # m/s, hour, degrees

    def define_data(self, wind_speed=None, time=None, wind_direction=None):
        """
        Inputs:
            wind_speed          list of wind speeds
            time                list of times at which the wind speeds were recorded
            wind_direction      list of wind directions
        """
        self.wind_speed = wind_speed
        self.time = time
        self.wind_direction = wind_direction


class RepairData(DataFile):
    """
    RepairData is designed to store the costs of repairing leaks from a particular reference an associated notes.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None):
        """
        Inputs:
            notes           A string containing notes on the object created
            raw_file_name   path to a raw input file
        """
        DataFile.__init__(self, notes, raw_file_name)
        self.repair_costs = []

    def define_data(self, repair_costs=None):
        """
        Inputs:
            repair_costs    list of costs to repair leaks
        """
        if repair_costs is None:
            repair_costs = []
        self.repair_costs = repair_costs


class PNNLData(DataFile):
    """
    Class for storing PNNL spectra data
    """
    def __init__(self, notes, raw_file_name):
        DataFile.__init__(self, notes, raw_file_name)
        self.km = None
        self.nu = None

    def define_data(self, km, nu):
        self.km = km
        self.nu = nu


class HITRAN(DataFile):
    """
    Class for storing PNNL spectra data
    """
    def __init__(self, notes, raw_file_name):
        DataFile.__init__(self, notes, raw_file_name)
        self.nu = None
        self.s = None
        self.gamma = None
        self.temp = None
        self.lower_e = None

    def define_data(self, nu, s, gamma, temp, lower_e):
        self.nu = nu
        self.s = s
        self.gamma = gamma
        self.temp = temp
        self.lower_e = lower_e
