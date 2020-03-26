from .abstract_detection_method import DetectionMethod
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs
import numpy as np
import scipy.special


class TieredDetect(DetectionMethod):
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
        self.sites_per_day = 4  # defines first tier speed...can be a dict for multiple site types
        self.labor = 100  # dollars/hour

        # -------------- Detection Variables -------------
        self.mu = 0.474
        self.lam = 3.88
        self.mu2 = 0.00185
        self.lam2 = 2.23
        self.ophrs = {'begin': 800, "end": 1700}
        self.insurvey = False
        self.surveyed_index = 0
        self.prelim_survey_time = 0
        self.secondary_survey_time = 0
        self.secondary_comps_hr = 150
        set_kwargs_attrs(self, kwargs)
        # -------------- Set calculated parameters --------------
        self.work_time = (self.ophrs['end'] - self.ophrs['begin']) / 100  # hours/day on average
        self.survey_time = gas_field.n_sites / self.sites_per_day * self.work_time  # hours
        # time_factor accounts for the finite simulation size. The effective capital cost is
        # reduced in the simulation based on the ratio of the sites in the
        # simulation to the number of sites that could be surveyed if there were enough wells to use the tech every day.
        self.time_factor = self.survey_time / (self.survey_interval * self.work_time)
        self.secondary_survey_cost = np.zeros(time.n_timesteps)
        # leaks_per_timestep is calculated based on the survey speed and number of leaks at the beginning of each survey
        self.leaks_per_timestep = 0
        self.loglam = np.log(self.lam)
        self.logmu = np.log(self.mu)
        self.loglam2 = np.log(self.lam2)
        self.logmu2 = np.log(self.mu2)
        self.site_survey_index = 0
        self.comp_survey_index = 0
        self.repair_delay = 0  # days
        # -------------- Financial Properties --------------
        self.capital_0 = 0 * self.time_factor  # dollars (defaults to zero)
        self.maintenance_0 = self.capital_0 * 0.1  # dollars/year

        self.capital = np.zeros(time.n_timesteps)
        self.replacement_cap(time)

        # maintenance costs are estimated as 10% of capital per year
        self.maintenance = [self.maintenance_0 * time.delta_t / 365, ] * time.n_timesteps  # $

        # survey_cost is the cost to survey all wells in the natural gas field
        self.survey_cost = self.labor * self.survey_time

        # find_cost is the cost of searching for leaks
        for ind in range(0, time.n_timesteps):
            curr_time = ind * time.delta_t
            if curr_time % self.survey_interval < time.delta_t:
                self.find_cost[ind] = self.survey_cost

    @staticmethod
    def time_in_ophrs(ct, et, op_begin, op_end):
        return max(min(et, op_end) - max(ct, op_begin), 0) / 100

    def sites_surveyed(self, time):
        """
        Computes the number of sites surveyed in a time step, allowing for time steps that are any fraction of a day
        or longer than a day and enforcing the operating hours defined for the detection technology
        :param time:
        :return nsites: the number of sites expected to be surveyed this time step
        """
        ct = 0
        nsites = 0
        while ct + 1 <= time.delta_t:
            ct += 1
            nsites += self.sites_per_day
        et = (time.delta_t - ct) * 2400
        ct = np.mod(time.current_time, 1) * 2400
        et += ct
        surv_hrs = self.time_in_ophrs(ct, et, self.ophrs['begin'], self.ophrs['end'])
        if et > 2400:
            ct = 0
            et -= 2400
            surv_hrs += self.time_in_ophrs(ct, et, self.ophrs['begin'], self.ophrs['end'])
        nsites += surv_hrs / self.work_time * self.sites_per_day
        return int(nsites) + np.random.binomial(1, np.mod(nsites, 1))

    def detection(self, time, gas_field):
        """
        The detection method applies a probability of detection curve
        Inputs:
            time        an object of type Time (defined in feast_classes)
            gas_field   an object of type GasField (defined in feast_classes)
        """
        if time.current_time % self.survey_interval < time.delta_t:
            self.insurvey = True
        if self.insurvey:
            end_survey = False
            n_sites_surveyed = self.sites_surveyed(time)
            if n_sites_surveyed > 0:
                if n_sites_surveyed + self.site_survey_index > gas_field.n_sites:
                    n_sites_surveyed = gas_field.n_sites - self.site_survey_index
                    end_survey = True
                site_flux = np.zeros(n_sites_surveyed)
                for site_ind in range(self.site_survey_index, self.site_survey_index + n_sites_surveyed):
                    site_flux[site_ind - self.site_survey_index] = \
                        np.sum(self.leaks.flux[self.leaks.site_index == site_ind])
                scores = np.random.uniform(0, 1, n_sites_surveyed)
                cond = site_flux > 0
                probs = np.zeros(n_sites_surveyed)
                probs[cond] = 0.5 + 0.5 * scipy.special.erf((np.log(site_flux[cond]) - self.logmu) /
                                                            (self.loglam * np.sqrt(2)))
                sites_flagged = np.where(scores < probs)[0] + self.site_survey_index
                self.prelim_survey_time += n_sites_surveyed / self.sites_per_day * self.work_time
                if end_survey:
                    self.insurvey = False
                    self.site_survey_index = 0
                self.site_survey_index += n_sites_surveyed
                secondary_survey_time = 0
                for site_ind in sites_flagged:
                    site_name = self.find_site_name(gas_field, site_ind)
                    n_comps = gas_field.sites[site_name]['parameters'].max_comp_ind
                    secondary_survey_time += n_comps / self.secondary_comps_hr
                    cond = np.where((self.leaks.site_index == site_ind) &
                                    (self.leaks.flux > 0))[0]
                    scores = np.random.uniform(0, 1, len(cond))
                    probs = 0.5 + 0.5 * scipy.special.erf((np.log(self.leaks.flux[cond]) - self.logmu2) /
                                                          (self.loglam2 * np.sqrt(2)))
                    detect = cond[scores < probs]
                    # self.repair_cost[time.time_index] += \
                    #     np.sum(np.random.choice(gas_field.repair_cost_dist.repair_costs,
                    #                             np.sum(self.leaks.reparable[detect])))
                    self.leaks.endtime[detect[self.leaks.reparable[detect]]] = time.current_time + self.repair_delay
                self.find_cost[time.time_index] += self.labor * secondary_survey_time
        self.null_detection(time, gas_field)
