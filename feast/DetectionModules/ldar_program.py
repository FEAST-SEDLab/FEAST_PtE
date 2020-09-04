"""
This module defines the LDARProgram class.
"""

import numpy as np
import copy
from .comp_survey import CompSurvey
from .repair import Repair


class LDARProgram:
    """
    An LDAR program contains any detection methods that will be applied as well as the repair method
    """
    def __init__(self, time, gas_field, tech_dict):
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        self.emissions = copy.deepcopy(gas_field.initial_emissions)
        self.emissions_timeseries = []
        self.vents_timeseries = []
        self.tech_dict = tech_dict
        self.repair = []
        for tech in tech_dict.values():
            if type(tech.dispatch_object) is Repair:
                self.repair.append(tech.dispatch_object)

    def end_emissions(self, time):
        """
        Sets the emissions to zero for all emissions that have exceeded their duration
        :param time: a Time object defined in simulation_classes
        :return repair_costs:
        """
        nem = self.emissions.n_leaks
        null_repair_cond = np.where((self.emissions.flux[:nem] > 0) &
                                    self.emissions.reparable[:nem] &
                                    (self.emissions.endtime[:nem] <= time.current_time))[0]
        self.repair_cost[time.time_index] += np.sum(self.emissions.repair_cost[null_repair_cond])

        cond = np.where(self.emissions.endtime[:self.emissions.n_leaks] <= time.current_time)[0]

        self.emissions.delete_leaks(cond)

    def action(self, time, gas_field):
        """
        Runs the detect method for every tech in tech_dict and runs the repair method
        :param time: the simulation time object
        :param gas_field: the simulation gas_field object
        :return:
        """
        for tech in self.tech_dict.values():
            if tech.survey_interval and np.mod(time.current_time, tech.survey_interval) < time.delta_t:
                tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
            tech.detect(time, gas_field, self.emissions, self.find_cost)
        for rep in self.repair:
            rep.repair(time, self.emissions)
        self.end_emissions(time)
