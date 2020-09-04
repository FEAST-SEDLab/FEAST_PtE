"""
This module defines the component level survey based detection class, CompDetect.
"""
import numpy as np
from .abstract_detection_method import DetectionMethod
from ..EmissionSimModules.simulation_functions import set_kwargs_attrs


class CompSurvey(DetectionMethod):
    """
    This class specifies a component level, survey based detection method.
    The class has three essential attributes:
    1.) An operating envelope function to determine if the detection method can be applied
    2.) A probability of detection surface function to determine which emissions are detected
    3.) The ability to call a follow up action
    """
    def __init__(self, time, **kwargs):
        DetectionMethod.__init__(self, time)
        self.dispatch_object = None

        # --------------- Process Variables -------------------
        self.ophrs = {'begin': 8, 'end': 17}
        self.survey_interval = None
        self.survey_speed = 150  # components/hr
        self.labor = 100  # $/hr
        self.op_env_wait_time = 7  # days (amount of time to wait for operating envelope conditions at a partly surveyed
                                   # site)

        # --------------- Detection Variables -----------------
        self.site_queue = []  # queue of sites to survey
        self.comp_survey_index = 0
        self.site_survey_index = None
        self.detection_probability_points = None
        self.detection_probabilities = None

        # -------------- Internal variables -----------------
        self.mid_site_fail_time = np.infty

        # Set all attributes defined in kwargs
        set_kwargs_attrs(self, kwargs, only_existing=True)

        # -------------- Set calculated parameters --------------
        work_time = (self.ophrs['end'] - self.ophrs['begin']) / 24
        self.comps_per_timestep = self.survey_speed * 24 * time.delta_t * np.min([1, work_time / time.delta_t])

    def detect_prob_curve(self, time, gas_field, cond, emissions):
        """
        This function determines which leaks are found given an array of indexes defined by "cond"
        In this case, the detected leaks are determined using a probability of detection curve
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
        vals = self.get_current_conditions(time, gas_field, emissions, cond)
        probs = self.empirical_interpolator(self.detection_probabilities, self.detection_probability_points, vals)
        detect = cond[scores <= probs]
        return detect

    def emitters_surveyed(self, time, gas_field, emissions, find_cost):
        """
        Determines which emitters are surveyed during the current time step.
        Accounts for the number of components surveyed per timestep, the number of components at each site, and the
        component and site at which the survey left off in the previous time step
        :param time:
        :param gas_field:
        :param emissions: emissions object
        :param find_cost: the find_cost array associated with the ldar program
        :return emitter_inds: indexes of emissions to evaluate at this timestep
        """
        remaining_comps = self.comps_per_timestep
        find_cost[time.time_index] += self.comps_per_timestep / self.survey_speed * self.labor
        # only considering nonzero leaks is important because cleaning leaks that have been set to 0 from the leak
        # set happens periodically in the simulation and would otherwise cause indexing errors
        emitter_inds = []
        while remaining_comps > 0:
            if (time.current_time - self.mid_site_fail_time) > self.op_env_wait_time:
                # if the survey has been stuck part way through a site for seven days due to operating envelope
                # conditions, move on to the next valid site.
                self.comp_survey_index = 0
            if self.comp_survey_index == 0:
                site_inds = self.choose_sites(gas_field, time, 1)
                if len(site_inds) > 0:
                    self.site_survey_index = site_inds[0]
                else:
                    self.site_survey_index = None
                    break
            else:
                op_env = self.check_op_envelope(gas_field, time, self.site_survey_index)
                if 'fail' in op_env:
                    if self.mid_site_fail_time > time.current_time:
                        self.mid_site_fail_time = time.current_time
                    break
            # site_name is used to flag all sites with the same properties
            # One site name can refer to multiple sites
            site_name = self.find_site_name(gas_field, self.site_survey_index)
            site_cond = emissions.site_index == self.site_survey_index
            comp_cond = (emissions.comp_index >= self.comp_survey_index) & \
                        (emissions.comp_index < self.comp_survey_index + remaining_comps)
            emitter_inds.extend(np.where(site_cond & comp_cond)[0])
            if remaining_comps + self.comp_survey_index > gas_field.sites[site_name]['parameters'].max_comp_ind:
                remaining_comps -= (gas_field.sites[site_name]['parameters'].max_comp_ind - self.comp_survey_index)
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
                detect = self.detect_prob_curve(time, gas_field, np.array(emitter_inds), emissions)
                # Deploy follow up action
                self.dispatch_object.action(None, detect)

    def action(self, site_inds=[], emit_inds=[]):
        """
        Action to add sites to queue. Expected to be called by another detection method or by an LDAR program
        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return:
        """
        self.site_queue.extend(site_inds)
