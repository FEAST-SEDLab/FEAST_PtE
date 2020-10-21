"""
This module defines all classes used to store input data.
"""


class DataFile:
    """
    DataFile is an abstract super class that all data file types inherit from FEAST.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None, data_prep_file=None):
        """
        :param notes: A string containing notes on the object created
        :param raw_file_name: path (or list of paths) to a raw input file(s)
        :param data_prep_file: path to the file used to process input data and create the DataFile object
        """
        self.notes = notes
        self.raw_file_name = raw_file_name
        self.data_prep_file = data_prep_file


class LeakData(DataFile):
    """
   LeakData is designed to store all leak size data from a reference. It accommodates multiple detection methods
   within a single instance.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None, data_prep_file=None, leak_sizes=None):
        """
        Creates a LeakData object

        :param notes: a string containing notes on the object created
        :param raw_file_name: path to a raw input file
        :param data_prep_file: path to a script used to build the object from the raw data file
        :param leak_sizes: list of leak sizes. If leaks were detected using multiple methods, leak_sizes must be a dict
                        with one key for each detection method
        """
        DataFile.__init__(self, notes, raw_file_name, data_prep_file=data_prep_file)
        self.leak_sizes = leak_sizes or dict()
        self.well_counts = dict()
        self.comp_counts = dict()

    def define_data(self, leak_data=None, well_counts=None, comp_counts=None, detect_methods=None):
        """
        Check data formatting and set keys...there is exactly one detection method, leak_data may be a list.

        :param leak_data: leak_data must be a dict of emission rates if there are multiple detection methods. If
            there is exactly one detection method, leak_data may be a list.
        :param well_counts: lists the number of wells inspected by each detection method in keys
        :param comp_counts: lists the number of components inspected by each detection method in keys
        :param detect_methods: each detection method should have a unique key associated with it
        :return: None
        """

        if detect_methods is None:
            # leak_data should be a dict if there are multiple detection methods. Each key will define a detection
            # method and the associated dictionary entry will be the list of fluxes from leaks found using that
            # detection method
            if type(leak_data) is dict:
                detect_methods = leak_data.keys()
            else:
                detect_methods = ['all_leaks']
        if detect_methods != leak_data.keys() or leak_data.keys() != well_counts.keys():
            error_str = "The 'keys' argument passed to LeakData.defineData(), leak_data.keys() and " \
                        "well_counts.keys() must all be equivalent"
            raise ValueError(error_str)

        for key in detect_methods:
            self.leak_sizes[key] = leak_data[key]
            self.well_counts[key] = well_counts[key]
            self.comp_counts[key] = comp_counts[key]


class RepairData(DataFile):
    """
    RepairData is designed to store the costs of repairing leaks from a particular reference andd associated notes.
    """
    def __init__(self, notes='No notes provided', raw_file_name=None):
        """
        :param notes: A string containing notes on the object created
        :param raw_file_name: path to a raw input file
        """
        DataFile.__init__(self, notes, raw_file_name)
        self.repair_costs = []

    def define_data(self, repair_costs=None):
        """
        :param repair_costs: list of costs to repair leaks
        """
        if repair_costs is None:
            repair_costs = []
        self.repair_costs = repair_costs


class ProductionData(DataFile):
    """
    Stores an array of production rates that may be associated with gas production sites
    """
    def __init__(self, site_prod=None, **kwargs):
        """
        :param prod_dat: an array of production rates
        :param kwargs: pass through
        """
        DataFile.__init__(self, kwargs)
        self.site_prod = site_prod
