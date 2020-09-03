import numpy as np
from .abstract_detection_method import DetectionMethod
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class SiteMonitor(DetectionMethod):
    """
    This class specifies a site level continuous monitoring method
    The method has three essential characteristics:
    1.) A list of the sites where the method applies
    2.) A time-to-detect surface specified as a list of conditions and associated mean detection times.
    3.) A dict of the variables used to define the time-to-detect surface (the attribute is specified in the
    DetectionMethod class inherited here)
    """
    def __init__(self, time, **kwargs):

        DetectionMethod.__init__(self, time)
        self.dispatch_object = None

        # --------------- Process Variables -------------------
        self.ophrs = {'begin': 0, 'end': 24}
        self.capital = 0
        self.site_queue = []  # list of sites where this method applies

        # --------------- Detection Variables -----------------
        self.time_to_detect_points = None  # conditions for which the mean time to detect is specified
        self.time_to_detect_days = None  # Mean detection times given the conditions listed in _points

        # -------------- Internal variables -----------------

        # Set all attributes defined in kwargs
        set_kwargs_attrs(self, kwargs, only_existing=True)
        self.time_to_detect_points = np.array(self.time_to_detect_points)
        self.time_to_detect_days = np.array(self.time_to_detect_days)
        # -------------- Set calculated parameters --------------

    @staticmethod
    def prob_detection(time, ttd):
        """
        Calculates the probability of detection during a timestep of length time.delta_t and a mean time to detection
        ttd
        :param time: Simulation time object
        :param ttd: mean time to detection (days)
        :return: the probability of detection during this timestep
        """
        if ttd == 0:
            return 1
        else:
            return 1 - np.exp(-time.delta_t / ttd)

    def detect_prob_curve(self, time, gas_field, site_inds, emissions):
        """
        Determines which sites are passed to the dispatch method
        In this case, the sites to pass are determined by calculating a probability of detection based on the
        simulation time resolution delta_t and the mean time to detection
        :param time: Simulation time object
        :param gas_field: Simulation gas_field object
        :param site_inds: The set of sites to be considered
        :param emissions: an object storing all emissions in the simulation
        :return detect: the indexes of detected leaks
        """
        n_scores = len(site_inds)
        if n_scores == 0:
            return site_inds
        probs = np.zeros(n_scores)
        counter = 0
        for site_ind in site_inds:
            vals = np.zeros(len(self.detection_variables))
            ind = 0
            cond = np.where(emissions.site_index[:emissions.n_leaks] == site_ind)[0]
            for v, im in self.detection_variables.items():
                if v in gas_field.met:
                    vals[ind] = gas_field.get_met(time, v, interp_modes=im, ophrs=self.ophrs)[v]
                else:
                    # sum all emission variables needed for detection
                    vals[ind] = np.sum(emissions.__getattribute__(v)[cond])
                ind += 1
            ttd = self.empirical_interpolator(self.time_to_detect_points, self.time_to_detect_days, vals)
            probs[counter] = self.prob_detection(time, ttd)
            counter += 1
        scores = np.random.uniform(0, 1, n_scores)
        detect = np.array(site_inds)[scores <= probs]
        return detect

    def detect(self, time, gas_field, emissions, find_cost):
        """
        The detection method implements a continuous monitor detection method model
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        # enforces the operating hours
        if self.check_time(time):
            # choose sites accounts for the operating envelope
            site_inds = self.choose_sites(gas_field, time, len(self.site_queue))
            if len(site_inds) > 0:
                detect = self.detect_prob_curve(time, gas_field, site_inds, emissions)
                # Deploy follow up action
                self.dispatch_object.action(detect, None)

    def action(self, site_inds=[], emit_inds=[]):
        """
        Action to add sites to queue. Expected to be called by another detection method or by an LDAR program
        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return:
        """
        self.site_queue.extend(site_inds)
