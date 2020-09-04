"""
This module defines Repair classes.
These classes may be called by detection objects as follow up actions.
"""

import numpy as np


class Repair:
    """
    Defines a repair process
    Includes two methods:
    __init__
    repair
    """
    def __init__(self, repair_delay=0):
        self.repair_delay = repair_delay
        self.to_repair = []

    def repair(self, time, emissions):
        """
        Adjusts the emission end time based on the current time and the repair delay time
        If the null emission end time comes before the repair time, the end time is not changed

        :param time: a Time object as defined in simulation classes
        :param emissions: an Emission object, as defined in simulation classes
        :return:
        """
        if len(self.to_repair) > 0:
            rep_cond = np.array(self.to_repair)[emissions.reparable[self.to_repair]]
            rep_cond = rep_cond[emissions.endtime[rep_cond] > time.current_time + self.repair_delay]
            emissions.endtime[rep_cond] = time.current_time + self.repair_delay
            self.to_repair = []

    def action(self, site_inds=None, emit_inds=None):
        self.to_repair.extend(emit_inds)
