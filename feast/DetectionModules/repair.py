"""
This module defines the Repair class. Repair may be called by detection objects as follow up actions.
"""

import numpy as np
from ..EmissionSimModules.result_classes import ResultDiscrete


class Repair:
    """
    Defines a repair process. A repair process determines when emissions are ended by an LDAR program and the
    associated costs.
    """
    def __init__(self, repair_delay=0, name=None):
        """
        :param repair_delay: The time between when an emission is passed to Repair and when it is removed from the
            simulation (float--days)
        """
        if name is None:
            name = 'repair'
        self.name = name
        self.repair_delay = repair_delay
        self.to_repair = []
        self.repair_count = ResultDiscrete(units='Count')
        self.repair_cost = ResultDiscrete(units='USD')

    def repair(self, time, emissions):
        """
        Adjusts the emission end time based on the current time and the repair delay time
        If the null emission end time comes before the repair time, the end time is not changed

        :param time: a Time object
        :param emissions: an Emission object
        :return: None
        """
        # todo: Check Null scenario repair costs
        if len(self.to_repair) > 0:
            rep_cond = emissions.emissions.reparable & emissions.emissions.index.isin(self.to_repair) & \
                       (emissions.emissions.end_time > time.current_time + self.repair_delay)
            emissions.emissions.loc[rep_cond, 'end_time'] = time.current_time + self.repair_delay
            time_cond = emissions.emissions['start_time'] >= emissions.emissions['end_time']
            emissions.emissions.loc[time_cond, 'start_time'] = emissions.emissions.loc[time_cond, 'end_time']
            self.repair_count.append_entry([time.current_time + self.repair_delay,
                                            len(emissions.emissions.loc[rep_cond].index.unique())])
            self.repair_cost.append_entry([time.current_time + self.repair_delay,
                                           np.sum(emissions.emissions.loc[rep_cond, 'repair_cost'])])
            self.to_repair = []

    def action(self, site_inds=None, emit_inds=None):
        """
        adds emissions to the to_repair queue.

        :param site_inds: not used
        :param emit_inds: A list of emission indexes to repair
        :return: None
        """
        self.to_repair.extend(emit_inds)
