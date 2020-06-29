"""
This module defines the component level survey based detection class, CompDetect.
"""
import numpy as np
import copy
from .repair import Repair
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class CompDetect:
    """
    This class specifies a component level, survey based detection method.
    The class has three essential attributes:
    1.) An operating envelope function to determine if the detection method can be applied
    2.) A probability of detection surface function to determine which emissions are detected
    3.) The ability to call a follow up action
    """
    def __init__(self, time, **kwargs):
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        self.dispatch_object = Repair()

        # --------------- Process Variables -------------------
        self.ophrs = {'begin': 800, 'end': 1700}
        self.survey_interval = None
        self.survey_speed = 150  # components/hr
        self.labor = 100  # $/hr

        # --------------- Detection Variables -----------------
        self.mu = 0.02  # g/s (median detection threshold)
        self.sigma = 0.80  # ln(g/s) (standard deviation of emission detection probability curve in log space)
        self.sites_to_survey = []  # queue of sites to survey
        self.comp_survey_index = 0
        self.site_survey_index = None

        # Set all attributes defined in kwargs
        set_kwargs_attrs(self, kwargs, only_existing=True)

        # -------------- Set calculated parameters --------------
        work_time = (self.ophrs['end'] - self.ophrs['begin']) / 2400
        self.comps_per_timestep = self.survey_speed * 24 * time.delta_t * np.min([1, work_time / time.delta_t])
        self.logmu = np.log(self.mu)

    def detect_prob_curve(self, cond, emissions):
        """
        This function determines which leaks are found given an array of indexes defined by "cond"
        In this case, the detect leaks are determined using a probability of detection curve
        :param cond: The set of indexes to be considered
        :param emissions: an object storing all emissions in the simulation
        :return detect: the indexes of detected leaks
        """
        n_scores = len(cond)
        if n_scores == 0:
            return cond
        cond = cond[emissions.flux[cond] > 0]
        n_scores = len(cond)
        scores = np.random.uniform(0, 1, n_scores)
        probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - self.logmu) / (self.sigma * np.sqrt(2))) for f
                                      in emissions.flux[cond]])
        detect = cond[scores <= probs]
        return detect

    @staticmethod
    def find_site_name(gas_field, site_index):
        for sitename, site in gas_field.sites.items():
            if site['parameters'].site_inds[0] <= site_index < site['parameters'].site_inds[1]:
                return sitename
        return -1

    @staticmethod
    def find_comp_name(gas_field, sitename, comp_index):
        for compname, comp in gas_field.sites[sitename]['parameters'].comp_dict.items():
            if comp.comp_inds[0] <= comp_index < comp.comp_inds[1]:
                return compname
        return -1

    def check_time(self, time):
        """
        Determines whether or not the detection method is active during the present time step
        :param time:
        :param gas_field:
        :return:
        """
        oktime = self.ophrs['begin'] <= np.mod(time.current_time, 1) * 2400 < self.ophrs['end']
        # accounts for a delta_t that is greater than the daily working hours
        timesize = time.delta_t * 2400 > self.ophrs['end'] - self.ophrs['begin']
        return oktime or timesize

    def emitters_surveyed(self, time, gas_field, emissions, find_cost):
        """
        Determines which emitters are surveyed during the current time step.
        Accounts for the number of components surveyed per timestep, the number of components at each site, and the
        component and site at which the survey left off in the previous time step
        :param time:
        :param gas_field:
        :param emissions: emissions object
        :param find_cost: the find_cost array associated with the ldar program
        :return emitter_inds: indexes of emissions to re-evaluate
        """
        remaining_comps = self.comps_per_timestep
        find_cost[time.time_index] += self.comps_per_timestep / self.survey_speed * self.labor
        # only considering nonzero leaks is important because cleaning leaks that have been set to 0 from the leak
        # set happens periodically in the simulation and would otherwise cause indexing errors
        emitter_inds = []
        while remaining_comps > 0:
            if self.site_survey_index is None:
                self.site_survey_index = self.sites_to_survey.pop(0)
            # site_name is used to flag all sites with the same properties
            # One site name can refer to multiple sites
            site_name = self.find_site_name(gas_field, self.site_survey_index)
            emitter_inds.extend(np.where((emissions.site_index == self.site_survey_index) &
                                (emissions.comp_index >= self.comp_survey_index) &
                                (emissions.comp_index < self.comp_survey_index + remaining_comps))[0])
            if remaining_comps + self.comp_survey_index > gas_field.sites[site_name]['parameters'].max_comp_ind:
                remaining_comps -= (gas_field.sites[site_name]['parameters'].max_comp_ind - self.comp_survey_index)
                self.comp_survey_index = 0
                if len(self.sites_to_survey) > 0:
                    self.site_survey_index = self.sites_to_survey.pop(0)
                else:
                    self.site_survey_index = None
                    self.comp_survey_index = 0
            else:
                self.comp_survey_index += remaining_comps
                remaining_comps = 0
        return emitter_inds

    def detect(self, time, gas_field, emissions, find_cost):
        """
        The detection method implements a survey-based detection method model
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        # enforces the operating hours
        if self.check_time(time):
            emitter_inds = self.emitters_surveyed(time, gas_field, emissions, find_cost)
            if len(emitter_inds) > 0:
                detect = self.detect_prob_curve(np.array(emitter_inds), emissions)
                # Deploy follow up action
                self.dispatch_object.action(None, detect)

    def action(self, site_inds=[], emit_inds=[]):
        """
        Action to add sites to queue. Expected to be called by another detection method or by an LDAR program
        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return:
        """
        self.sites_to_survey.extend(site_inds)
