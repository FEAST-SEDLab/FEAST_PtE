"""
    simulation_classes stores the classes used to represent time, results and financial settings in simulations.
"""


class Time:
    """
    Instances of the time class store all time related information during a simulation
    """
    def __init__(self, delta_t=3650/4000, end_time=10 * 365, current_time=0):
        """
        Inputs:
        delta_t                 length of one timestep (days)
        end_time                length of the simulation (days)
        current_time            current time in a simulation (days)
        """
        self.n_timesteps = int(end_time/delta_t+1)
        self.end_time = end_time
        self.delta_t = delta_t
        self.current_time = current_time
        self.time_index = 0


# FinanceSettings stores all parameters relating to economic calculations
class FinanceSettings:
    def __init__(self, gas_price=2E-4, discount_rate=0.08):
        self.gas_price = gas_price  # dollars/gram (2e-4 $/g=$5/mcf methane at STP)
        self.discount_rate = discount_rate


class Results:
    """
    Class in which to save results
    """
    def __init__(self, time, gas_field, ldar_program_dict, econ_set):
        """
        Inputs:
        time                    Time object
        gas_field               GasField object
        ldar_program_dict       dict of detection methods and associated data
        econ_set                Economic settings defined for the simulation
        """
        self.time = time
        self.gas_field = gas_field
        self.ldar_program_dict = ldar_program_dict
        self.econ_settings = econ_set


