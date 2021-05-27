"""
This module defines the LDARProgram class.
"""

import numpy as np
import copy
from .repair import Repair
from ..EmissionSimModules.result_classes import ResultDiscrete


class LDARProgram:
    """
    An LDAR program contains one or more detection methods and one or more repair methods. Each LDAR program records
    the find and repair costs associated with all detection and repair methods in the program. The LDAR program
    deploys runs the action methods of each detection and repair method contained in the program. The detection and
    repair methods determine their own behavior at each time step.
    """
    def __init__(self, gas_field, tech_dict):
        """
        :param gas_field: a GasField object
        :param tech_dict: a dict containing all of the detection methods to be employed by the LDAR program. The dict
            must have the form {"name": DetectionMethod}. All of the relationships between detection methods and between
            detection methods and repair methods must be defined by the dispatch_objects specified for each method.
        """
        self.emissions = copy.deepcopy(gas_field.emissions)
        self.emissions_timeseries = []
        self.vents_timeseries = []
        self.tech_dict = tech_dict
        self.repair = {}
        self.repair_cost = ResultDiscrete(units='USD')
        for tech_name, tech in tech_dict.items():
            if type(tech.dispatch_object) is Repair:
                self.repair[tech_name + ' ' + tech.dispatch_object.name] = tech.dispatch_object

    def action(self, time, gas_field):
        """
        Runs the detect method for every tech in tech_dict and runs the repair method
        :param time: the simulation time object
        :param gas_field: the simulation gas_field object
        :return:
        """
        for i, tech in enumerate(self.tech_dict.values()):
            # Check for tiered detection
            if (len(self.tech_dict.values()) > 1) and (i >= 1):
                det_ct = self.tech_dict[list(self.tech_dict.keys())[i - 1]].detection_count.time_value
                if len(det_ct) > 0:
                    det_ct_time = det_ct[-1][0]
                else:
                    det_ct_time = None
                if (det_ct_time == time.current_time) or (det_ct_time is not None):
                    if hasattr(tech, 'survey_interval') and tech.survey_interval \
                            and np.mod(time.current_time, tech.survey_interval) < time.delta_t:
                        tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
                    tech.detect(time, gas_field, self.emissions.get_current_emissions(time))
            else:
                if hasattr(tech, 'survey_interval') and tech.survey_interval \
                        and np.mod(time.current_time, tech.survey_interval) < time.delta_t:
                    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
                tech.detect(time, gas_field, self.emissions.get_current_emissions(time))
        for rep in self.repair.values():
            rep.repair(time, self.emissions)

    def calc_rep_costs(self, time):
        """
        Calculates the total repair costs up to time.current_time, assuming that all reparable emissions that have a
        max end_time less than time.current_time have been repaired.

        :param time: a FEAST time object
        :return: None
        """
        for em in self.emissions.emissions.index.unique():
            empdf_temp = self.emissions.emissions.loc[[em]]
            max_row = empdf_temp[empdf_temp.end_time == empdf_temp.end_time.max()].iloc[0]
            if max_row.reparable & (max_row.end_time < time.current_time):
                self.repair_cost.append_entry([max_row.end_time, max_row.repair_cost])
