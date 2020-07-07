import copy
import numpy as np
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class DetectionMethod:
    """
    DetectionMethod is an abstract super class that defines the form required for all detection methods
    """
    def __init__(self, time, **kwargs):
        """
        Inputs:
            time         a time object (Defined in simulation_classes)
            kwargs       optional input dict that will override default parameters
        """
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        # Set all attributes defined in kwargs, regardless of whether they already exist
        set_kwargs_attrs(self, kwargs, only_existing=True)

    def check_time(self, time):
        """
        Determines whether or not the detection method is active during the present time step
        :param time:
        :param gas_field:
        :return:
        """
        oktime = self.ophrs['begin'] <= np.mod(time.current_time, 1) * 24 < self.ophrs['end']
        # accounts for a delta_t that is greater than the daily working hours
        timesize = time.delta_t * 24 > self.ophrs['end'] - self.ophrs['begin']
        return oktime or timesize
