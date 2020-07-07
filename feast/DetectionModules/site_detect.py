"""
This module defines the component level survey based detection class, CompDetect.
"""
import numpy as np
from .abstract_detection_method import DetectionMethod
from .comp_detect import CompDetect
from .repair import Repair
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class SiteDetect(DetectionMethod):
    """
    This class specifies a site level, survey based detection method.
    The class has three essential attributes:
    1.) An operating envelope function to determine if the detection method can be applied
    2.) A probability of detection surface function to determine which emissions are detected
    3.) The ability to dispatch a follow up action
    """
    def __init__(self, time, **kwargs):
        self.dispatch_object = CompDetect(time)

        # --------------- Process Variables -------------------
        self.ophrs = {'begin': 8, 'end': 17}
        self.survey_interval = None
        self.sites_per_day = 200  # sites_per_day
        self.site_cost = 100  # $/site

        # --------------- Detection Variables -----------------
        self.mu = 0.474  # g/s (emission size with 50% probability of detection)
        self.sigma = 1.36  # ln(g/s) (standard deviation of emission detection probability curve in log space)
        self.sites_to_survey = []  # queue of sites to survey
        self.site_survey_index = None

        # Set all attributes defined in kwargs
        set_kwargs_attrs(self, kwargs, only_existing=True)

        # -------------- Set calculated parameters --------------
        work_time = (self.ophrs['end'] - self.ophrs['begin']) / 24
        self.sites_per_timestep = int(self.sites_per_day * time.delta_t * np.min([1, time.delta_t / work_time]))
        if self.sites_per_timestep < 1 and self.sites_per_Day > 0:
            print("WARNING: expecting less than 1 site surveyed per timestep. May lead to unexpected behavior.")
        self.logmu = np.log(self.mu)

    def detect_prob_curve(self, site_inds, emissions):
        """
        This function determines which leaks are found given an array of indexes defined by "cond"
        In this case, the detect leaks are determined using a probability of detection curve
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
            cond = np.where(emissions.site_index[:emissions.n_leaks] == site_ind)[0]
            site_flux = np.sum(emissions.flux[cond])
            if site_flux > 0:
                probs[counter] = 0.5 + 0.5 * np.math.erf((np.log(site_flux) - self.logmu) / (self.sigma * np.sqrt(2)))
            counter += 1
        scores = np.random.uniform(0, 1, n_scores)
        detect = np.array(site_inds)[scores <= probs]
        return detect

    def sites_surveyed(self, time, find_cost):
        """
        Determines which sites are surveyed during the current time step.
        Accounts for the number of sites surveyed per timestep
        :param time:
        :param find_cost: the find_cost array associated with the ldar program
        :return site_inds: indexes of sites scanned at this timestep
        """
        n_sites = np.min([self.sites_per_timestep, len(self.sites_to_survey)])
        find_cost[time.time_index] += n_sites * self.site_cost
        site_inds = []
        for ind in range(n_sites):
            site_inds.append(self.sites_to_survey.pop(0))
        return site_inds

    def detect(self, time, gas_field, emissions, find_cost):
        """
        The detection method implements a survey-based detection method model
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        # enforces the operating hours
        if self.check_time(time):
            site_inds = self.sites_surveyed(time, find_cost)
            if len(site_inds) > 0:
                detect = self.detect_prob_curve(site_inds, emissions)
                # Deploy follow up action
                self.dispatch_object.action(detect, None)

    def action(self, site_inds=[], emit_inds=[]):
        """
        Action to add sites to queue. Expected to be called by another detection method or by an LDAR program
        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return:
        """
        self.sites_to_survey.extend(site_inds)
