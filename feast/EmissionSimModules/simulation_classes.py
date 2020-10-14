"""
    simulation_classes stores the classes used to represent time, results and financial settings in simulations.
"""
import os
import pickle


class Time:
    """
    Instances of the time class store all time related information during a simulation
    """
    def __init__(self, delta_t=1, end_time=365, current_time=0):
        """
        :param delta_t: length of one timestep (days)
        :param end_time: length of the simulation (days)
        :param current_time: current time in a simulation (days)
        """
        self.n_timesteps = int(end_time/delta_t+1)
        self.end_time = end_time
        self.delta_t = delta_t
        self.current_time = current_time
        self.time_index = 0


# FinanceSettings stores all parameters relating to economic calculations
class FinanceSettings:
    """
    Finance settings stores the natural gas price and discount rate to be used in cost modeling.
    """
    def __init__(self, gas_price=0, discount_rate=0):
        """
        :param gas_price: The value of recovered gas ($/gram). 2e-4 $/g=$5/mcf methane at STP.
        :param discount_rate: The rate at which future cash flows are discounted
        """
        self.gas_price = gas_price
        self.discount_rate = discount_rate


class Results:
    """
    Class in which to save results
    """
    def __init__(self, time, gas_field, ldar_program_dict, econ_set):
        """
        :param time: Time object
        :param gas_field: GasField object
        :param ldar_program_dict: dict of detection methods and associated data
        :param econ_set: Economic settings defined for the simulation
        """
        self.time = time
        self.gas_field = gas_field
        self.ldar_program_dict = ldar_program_dict
        self.econ_settings = econ_set

    def save(self, dir_out):
        """
        Save results to a file

        :param dir_out: Name of directory in which to save output file.
        """

        if not os.path.exists(dir_out):
            os.makedirs(dir_out)
        n_realization = len(os.listdir(dir_out))
        file_out = dir_out + '/realization' + str(n_realization) + '.p'
        pickle.dump(self, open(file_out, 'wb'))


