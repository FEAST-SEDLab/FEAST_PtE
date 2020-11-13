"""
The site_monitor module defines the site monitor detection class, SiteMonitor.
"""

import numpy as np
from .abstract_detection_method import DetectionMethod

class SiteMonitor(DetectionMethod):
    """
    This class specifies a site level continuous monitoring method. A site monitor continuously observes
    emissions from an entire site and determines when an action should be dispatched at the site.
    The method has three essential characteristics:

    1. A list of the sites where the method applies
    2. A time-to-detect surface specified as a list of conditions and associated mean detection times
    3. The ability to dispatch a follow up action
    """
    def __init__(self, time, dispatch_object, time_to_detect_points, time_to_detect_days,
                 ophrs=None, capital=0, site_queue=None, **kwargs):
        """
        :param time: a Time object
        :param dispatch_object: the object to dispatch for follow up actions
        :param time_to_detect_points: The conditions at which the time to detection was measured. (NxM array,
            where N is the number of distinct conditions and M is the number of variables (up to two)).
        :param time_to_detect_days: The list of probabilities of detection associated with every point in
            detection_probability_points (array of shape N, where N is the number of conditions with an associated
            probability of detection).
        :param ophrs: The times of day when the SiteMonitor is operational. Should be a dict:\n
            {'begin': hour integer, 'end': hour integer}
        :param capital: The total cost of installing the site monitor system in the simulation (float--$)
        :param site_queue: A list of sites where the site monitor system is installed
        """

        DetectionMethod.__init__(self, time, **kwargs)
        if self.ophrs is None:
            self.ophrs = {'begin': 0, 'end': 24}
        self.dispatch_object = dispatch_object

        # --------------- Process Variables -------------------
        self.ophrs = ophrs
        self.deployment_cost.append_entry([time.current_time, capital])
        self.site_queue = site_queue or []  # list of sites where this method applies

        # --------------- Detection Variables -----------------
        self.time_to_detect_points = time_to_detect_points  # conditions for which the mean time to detect is specified
        self.time_to_detect_days = time_to_detect_days  # Mean detection times given the conditions listed in _points

        # -------------- Internal variables -----------------
        if type(self.site_queue) is not list:
            self.site_queue = list(self.site_queue)
        self.time_to_detect_points = np.array(self.time_to_detect_points)
        self.time_to_detect_days = np.array(self.time_to_detect_days)

    @staticmethod
    def prob_detection(time, ttd):
        """
        Calculates the probability of detection during a timestep of length time.delta_t and a mean time to detection
        ttd

        :param time: Simulation time object
        :param ttd: mean time to detection (float--days)
        :return: the probability of detection during this timestep
        """
        if ttd == 0:
            return 1
        else:
            return 1 - np.exp(-time.delta_t / ttd)

    def detect_prob_curve(self, time, gas_field, site_inds, emissions):
        """
        Determines which sites are passed to the dispatch method.
        In this case, the sites to pass are determined by calculating a probability of detection based on the
        simulation time resolution (time.delta_t) and the mean time to detection

        :param time: simulation Time object
        :param gas_field: simulation GasField object
        :param site_inds: the set of sites to be considered
        :param emissions: an object storing all emissions in the simulation
        :return detect: the indexes of detected leaks
        """
        n_scores = len(site_inds)
        if n_scores == 0:
            return site_inds
        probs = np.zeros(n_scores)
        counter = 0
        for site_ind in site_inds:
            vals = np.zeros([1, len(self.detection_variables)])
            ind = 0
            for v, im in self.detection_variables.items():
                if v in gas_field.met:
                    vals[0, ind] = gas_field.get_met(time, v, interp_modes=im, ophrs=self.ophrs)[v]
                else:
                    # sum all emission variables needed for detection
                    vals[0, ind] = np.sum(emissions[v][emissions.site_index == site_ind])
                ind += 1
            ttd = self.empirical_interpolator(self.time_to_detect_points, self.time_to_detect_days, vals)
            probs[counter] = self.prob_detection(time, ttd)
            counter += 1
        scores = np.random.uniform(0, 1, n_scores)
        detect = np.array(site_inds)[scores <= probs]
        return detect

    def detect(self, time, gas_field, emissions):
        """
        The detection method implements a continuous monitor detection method model

        :param time: a Time object
        :param gas_field: a GasField object
        :param emissions: a DataFrame containing emission data to evaluate
        :return: None
        """
        # enforces the operating hours
        if self.check_time(time):
            # choose sites accounts for the operating envelope
            site_inds = self.choose_sites(gas_field, time, len(self.site_queue), clear_sites=False)
            if len(site_inds) > 0:
                detect = self.detect_prob_curve(time, gas_field, site_inds, emissions)
                if len(detect) > 0:
                    self.detection_count.append_entry([time.current_time, len(detect)])
                # Deploy follow up action
                self.dispatch_object.action(detect, None)

    def action(self, site_inds=None, emit_inds=None):
        """
        Action to add sites to queue. Expected to be called by another detection method or by an LDAR program

        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return: None
        """
        self.extend_site_queue(site_inds)
