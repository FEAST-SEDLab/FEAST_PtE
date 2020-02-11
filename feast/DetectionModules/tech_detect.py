from .abstract_detection_method import DetectionMethod
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
import numpy as np
from . import helper_functions


class TechDetect(DetectionMethod):
    """
    This class specifies a general detection method.
    The detection method relies on a "probability of detection curve" that is independent of the a specific plume model.
    It includes a detection method and several parameters.
    """
    def __init__(self, time, gas_field, **kwargs):
        """
        Inputs:
              gas_field    a gas_field object (Defined in feast_classes)
              time         a time object (Defined in feast_classes)
              kwargs       optional input dictionary that will override default parameters
        """
        DetectionMethod.__init__(self, time, gas_field)
        # -------------- Hardware variables --------------
        self.lifetime = 10 * 365  # days

        # -------------- Process Variables --------------
        self.survey_interval = 365*0.5  # days
        self.survey_speed = 150  # components/hour
        self.drive_speed = 15  # m/s
        self.setup_time = 0.1  # hours
        self.labor = 100  # dollars/hour

        # -------------- Detection Variables -------------
        self.mu = 0.0185  # g/s
        self.lam = 0.9  # g/s
        self.ophrs = {'begin': 800, "end": 1700}
        self.insurvey = False
        self.site_survey_index = 0
        self.comp_survey_index = 0
        set_kwargs_attrs(self, kwargs)
        # -------------- Set calculated parameters --------------
        self.survey_time = (gas_field.n_comps / self.survey_speed +
                            gas_field.site_spacing / self.drive_speed / 3600 * gas_field.n_sites + self.setup_time *
                            gas_field.n_sites)  # Time per survey
        work_time = (self.ophrs['end'] - self.ophrs['begin']) / 2400
        self.comps_per_timestep = self.survey_speed * 24 * time.delta_t * np.min([1, work_time / time.delta_t])
        # hours per site
        # time_factor accounts for the finite simulation size. The effective capital cost is
        # reduced in the simulation based on the ratio of the wells in the
        # simulation to the number of wells that a single FID could survey.
        self.time_factor = self.survey_time / (self.survey_interval * work_time)
        # leaks_per_timestep is calculated based on the survey speed and number of leaks at the beginning of each survey
        self.leaks_per_timestep = 0
        self.loglam = np.log(self.lam)
        self.logmu = np.log(self.mu)
        # -------------- Financial Properties --------------
        self.capital_0 = 0 * self.time_factor  # dollars (defaults to zero)
        self.maintenance_0 = self.capital_0 * 0.1  # dollars/year

        self.capital = np.zeros(time.n_timesteps)
        helper_functions.replacement_cap(time, self)

        # maintenance costs are estimated as 10% of capital per year
        self.maintenance = [self.maintenance_0 * time.delta_t / 365, ] * time.n_timesteps  # $

        # survey_cost is the cost to survey all wells in the natural gas field
        self.survey_cost = self.labor * self.survey_time

        # find_cost is the cost of searching for leaks
        for ind in range(0, time.n_timesteps):
            curr_time = ind * time.delta_t
            if curr_time % self.survey_interval < time.delta_t:
                self.find_cost[ind] = self.survey_cost

    def detect_prob_curve(self, cond):
        """
        This function determines which leaks are found given an array of indexes defined by "cond"
        In this case, the detect leaks are determined using a probability of detection curve
        :param cond: The set of indexes to be considered
        :return detect: the indexes of detected leaks
        """
        n_scores = len(cond)
        if n_scores == 0:
            return cond
        cond = cond[self.leaks.flux[cond] > 0]
        n_scores = len(cond)
        scores = np.random.uniform(0, 1, n_scores)
        probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - self.logmu) / (self.loglam * np.sqrt(2))) for f
                                      in self.leaks.flux[cond]])
        detect = cond[scores <= probs]
        return detect

    def detection(self, time, gas_field):
        """
        The detection method implements a survey-based detection method model
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        self.null_detection(time, gas_field)
        self.survey_detection(time, gas_field, self.detect_prob_curve)
