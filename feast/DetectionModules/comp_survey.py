"""
This module defines the component level survey based detection class, CompSurvey.
"""
import numpy as np
from .abstract_detection_method import DetectionMethod


class CompSurvey(DetectionMethod):
    """
    This class specifies a component level, survey based detection method.
    A component level method identifies the specific component that is the source of the
    emission at the time of detection. Examples of components include connectors, valves, open ended lines, etc.
    A survey based method inspects emissions at specific moments in time (as opposed to a monitor method that
    continuously monitors for emissions).
    The class has three essential attributes:

    1. An operating envelope function to determine if conditions satisfy requirements for the method to be deployed
    2. A probability of detection surface function to determine which emissions are detected
    3. The ability to call a follow up action
    """
    def __init__(self, time, dispatch_object, survey_interval, survey_speed, labor, site_queue,
                 detection_probability_points, detection_probabilities,
                 comp_survey_index=0, site_survey_index=0,
                 op_env_wait_time=7, ophrs={'begin': 8, 'end': 17}, **kwargs):
        """
        :param time: a Time object
        :param dispatch_object: an object to dispatch for follow-up actions (typically a Repair method)
        :param survey_interval: Time between surveys (float--days)
        :param survey_speed: Speed of surveys (float--components/hr)
        :param labor: Cost of surveys (float--$/hr)
        :param site_queue: Sites to survey (list of ints)
        :param detection_probability_points: Set of conditions at which the probability was measured (array)
        :param detection_probabilities: Set of probabilities that were measured (array)
        :param comp_survey_index: Index of the component to be surveyed next (int)
        :param site_survey_index: Index of the site to be surveyed next (or currently under survey) (int)
        :param op_env_wait_time: Time to wait if operating envelope conditions fail part way through a site before
            moving to the next site (float--days)
        :param ophrs: range of times of day when the method performs survey (dict--hours). eg: {'begin': 8, 'end': 17}
        """
        DetectionMethod.__init__(self, time, **kwargs)
        self.dispatch_object = dispatch_object

        # --------------- Process Variables -------------------
        self.ophrs = ophrs
        self.survey_interval = survey_interval
        self.survey_speed = survey_speed  # components/hr
        self.labor = labor  # $/hr
        # days (amount of time to wait for operating envelope conditions at a partly surveyed site)
        self.op_env_wait_time = op_env_wait_time

        # --------------- Detection Variables -----------------
        self.site_queue = site_queue  # queue of sites to survey
        self.comp_survey_index = comp_survey_index
        self.site_survey_index = site_survey_index
        self.detection_probability_points = np.array(detection_probability_points)
        self.detection_probabilities = np.array(detection_probabilities)

        # -------------- Internal variables -----------------
        self.mid_site_fail_time = np.infty

        # -------------- Set calculated parameters --------------
        work_time = (self.ophrs['end'] - self.ophrs['begin']) / 24
        self.comps_per_timestep = self.survey_speed * 24 * time.delta_t * np.min([1, work_time / time.delta_t])

    def detect_prob_curve(self, time, gas_field, cond, emissions):
        """
        This function determines which leaks are found given an array of indexes defined by "cond."
        The method uses attributes of the DetectionMethod and interpolation.

        :param time: a Time object
        :param gas_field: a GasField object
        :param cond: The array of indexes to be considered (array of ints)
        :param emissions: an Emission object storing all emissions in the simulation
        :return detect: the indexes of detected leaks (array of ints)
        """
        n_scores = len(cond)
        if n_scores == 0:
            return cond
        cond = cond[emissions.flux[cond] > 0]
        n_scores = len(cond)
        scores = np.random.uniform(0, 1, n_scores)
        vals = self.get_current_conditions(time, gas_field, emissions, cond)
        probs = self.empirical_interpolator(self.detection_probability_points, self.detection_probabilities, vals)
        detect = cond[scores <= probs]
        return detect

    def emitters_surveyed(self, time, gas_field, emissions, find_cost):
        """
        Determines which emitters are surveyed during the current time step.
        Accounts for the number of components surveyed per timestep, the number of components at each site, and the
        component and site at which the survey left off in the previous time step

        :param time: a Time object
        :param gas_field: a GasField object
        :param emissions: an Emissions object
        :param find_cost: the find_cost array associated with the ldar program
        :return emitter_inds: indexes of emissions to evaluate at this timestep (list of ints)
        """
        remaining_comps = self.comps_per_timestep
        find_cost[time.time_index] += self.comps_per_timestep / self.survey_speed * self.labor
        # only considering nonzero leaks is important because cleaning leaks that have been set to 0 from the leak
        # set happens periodically in the simulation and would otherwise cause indexing errors
        emitter_inds = []
        while remaining_comps > 0:
            if (time.current_time - self.mid_site_fail_time) > self.op_env_wait_time:
                # if the survey has been stuck part way through a site for op_env_wait_time due to operating envelope
                # conditions, move on to the next valid site.
                self.comp_survey_index = 0
                self.mid_site_fail_time = np.infty
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
                    # TODO: append time.current_time, n site failures
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
        The detect method checks that the current time is within operating hours, selects emitters to inspect,
        determines which emissions are detected and dispatches the follow up action for detected emissions.

        :param time: a Time object
        :param gas_field: a GasField object
        """
        # enforces the operating hours
        if self.check_time(time):
            emitter_inds = self.emitters_surveyed(time, gas_field, emissions, find_cost)
            # TODO: append time.current_time, number of components surveyed to self.comp_survey_count
            # TODO append time.current_time, cost at this time step
            # TODO: append time.current_time, number of sites surveyed to self.site_survey_count
            if len(emitter_inds) > 0:
                detect = self.detect_prob_curve(time, gas_field, np.array(emitter_inds), emissions)
                # TODO: append time.current_time, len(detect) to self.detection_count
                # Deploy follow up action
                self.dispatch_object.action(None, detect)

    def action(self, site_inds=None, emit_inds=None):
        """
        Adds sites to the queue for future inspections. This method is expected to be called by another detection
        method or by an LDAR program.

        :param site_inds: List of sites to add to the queue
        :param emit_inds: Not used.
        :return: None
        """
        self.extend_site_queue(site_inds)
