"""
This module defines the Repair class. Repair may be called by detection objects as follow up actions.
"""

import numpy as np


class Repair:
    """
    Defines a repair process. A repair process determines when emissions are ended by an LDAR program and the
    associated costs.
    """
    def __init__(self, repair_delay=0):
        """
        :param repair_delay: The time between when an emission is passed to Repair and when it is removed from the
            simulation (float--days)
        """
        self.repair_delay = repair_delay
        self.to_repair = []

    def repair(self, time, emissions):
        """
        Adjusts the emission end time based on the current time and the repair delay time
        If the null emission end time comes before the repair time, the end time is not changed

        :param time: a Time object
        :param emissions: an Emission object
        :return: None
        """
        if len(self.to_repair) > 0:
            rep_cond = np.array(self.to_repair)[emissions.reparable[self.to_repair]]
            rep_cond = rep_cond[emissions.endtime[rep_cond] > time.current_time + self.repair_delay]
            emissions.endtime[rep_cond] = time.current_time + self.repair_delay
            self.to_repair = []

    def action(self, site_inds=None, emit_inds=None):
        """
        adds emissions to the to_repair queue.

        :param site_inds: not used
        :param emit_inds: A list of emission indexes to repair
        :return: None
        """
        self.to_repair.extend(emit_inds)
