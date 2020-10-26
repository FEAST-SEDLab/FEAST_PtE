"""
This module defines the LDARProgram class.
"""

import numpy as np
import copy
from .repair import Repair


class LDARProgram:
    """
    An LDAR program contains one or more detection methods and one or more repair methods. Each LDAR program records
    the find and repair costs associated with all detection and repair methods in the program. The LDAR program
    deploys runs the action methods of each detection and repair method contained in the program. The detection and
    repair methods determine their own behavior at each time step.
    """
    def __init__(self, time, gas_field, tech_dict):
        """
        :param time: a Time object
        :param gas_field: a GasField object
        :param tech_dict: a dict containing all of the detection methods to be employed by the LDAR program. The dict
            must have the form {"name": DetectionMethod}. All of the relationships between detection methods and between
            detection methods and repair methods must be defined by the dispatch_objects specified for each method.
        """
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        self.emissions = copy.deepcopy(gas_field.emissions)
        self.emissions_timeseries = []
        self.vents_timeseries = []
        self.tech_dict = tech_dict
        self.repair = []
        for tech in tech_dict.values():
            if type(tech.dispatch_object) is Repair:
                self.repair.append(tech.dispatch_object)

    # def end_emissions(self, time):
    #     """
    #     Sets the emissions to zero for all emissions that have exceeded their duration
    #     :param time: a Time object defined in simulation_classes
    #     :return repair_costs:
    #     """
    #     nem = self.emissions.n_em
    #     null_repair_cond = np.where((self.emissions.flux[:nem] > 0) &
    #                                 self.emissions.reparable[:nem] &
    #                                 (self.emissions.endtime[:nem] <= time.current_time))[0]
    #     self.repair_cost[time.time_index] += np.sum(self.emissions.repair_cost[null_repair_cond])
    #
    #     cond = np.where(self.emissions.endtime[:self.emissions.n_em] <= time.current_time)[0]
    #
    #     self.emissions.delete_leaks(cond)

    def action(self, time, gas_field):
        """
        Runs the detect method for every tech in tech_dict and runs the repair method
        :param time: the simulation time object
        :param gas_field: the simulation gas_field object
        :return:
        """
        for tech in self.tech_dict.values():
            if hasattr(tech, 'survey_interval') and tech.survey_interval \
                    and np.mod(time.current_time, tech.survey_interval) < time.delta_t:
                tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
            tech.detect(time, gas_field, self.emissions.get_current_emissions(time), self.find_cost)
        for rep in self.repair:
            rep.repair(time, self.emissions)
