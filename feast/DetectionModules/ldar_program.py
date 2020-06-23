"""
This module defines the LDARProgram class.
"""

import numpy as np
import copy
from .comp_detect import CompDetect
from .repair import Repair


class LDARProgram:
    """
    An LDAR program contains any detection methods that will be applied as well as the repair method
    """
    def __init__(self, time, gas_field, tech_list=[], repair=Repair()):
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        self.emissions = copy.deepcopy(gas_field.initial_emissions)
        self.emissions_timeseries = []
        self.vents_timeseries = []
        self.tech_list = tech_list
        self.repair = repair
        if not self.tech_list:
            self.tech_list.append(CompDetect(time, gas_field))

    def end_emissions(self, time):
        """
        Sets the emissions to zero for all emissions that have exceeded their duration
        :param time: a Time object defined in simulation_classes
        :return repair_costs:
        """
        nem = self.emissions.n_leaks
        null_repair_cond = (self.emissions.flux[:nem] > 0) & \
            self.emissions.reparable[:nem] & \
            (self.emissions.endtime[:nem] <= time.current_time)
        self.repair_cost[time.time_index] += np.sum(self.emissions.repair_cost[null_repair_cond])

        cond = np.where(self.emissions.endtime[:self.emissions.n_leaks] <= time.current_time)[0]

        self.emissions.delete_leaks(cond)

    def action(self, time, gas_field):
        """
        Runs the detect method for every tech in tech_list and runs the repair method
        :param time: the simulation time object
        :param gas_field: the simulation gas_field object
        :return:
        """
        for tech in self.tech_list:
            tech.detect(time, gas_field, self.emissions, self.find_cost)
        self.repair.repair(time, self.emissions)
        self.end_emissions(time)
