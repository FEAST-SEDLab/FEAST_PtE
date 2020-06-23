import copy
import numpy as np
from ..GeneralClassesFunctions.simulation_functions import set_kwargs_attrs


class DetectionMethod:
    """
    DetectionMethod is an abstract super class that defines the form required for all detection methods
    """
    def __init__(self, time, gas_field, notes=None, **kwargs):
        """
        Inputs:
            time         a time object (Defined in simulation_classes)
            gas_field    a gas_field object (Defined in simulation_classes)
            notes        a description of the object created
            kwargs       optional input dict that will override default parameters
        """
        self.notes = notes
        self.null_repaired = []
        self.repair_delay = 0  # days
        self.ophrs = {'begin': 0, 'end': 2400}
        self.survey_interval = 365  # days
        self.comps_per_timestep = 0
        self.repair_cost = np.zeros(time.n_timesteps)
        self.capital = np.zeros(time.n_timesteps)
        self.find_cost = np.zeros(time.n_timesteps)
        # leaks will be updated throughout simulations. initial_leaks should remain constant, so a copy is needed.
        self.leaks = copy.deepcopy(gas_field.initial_emissions)
        self.insurvey = False
        self.site_survey_index = 0
        self.comp_survey_index = 0
        self.count_found = []
        self.emissions = []
        self.vents = []
        # Set all attributes defined in kwargs, regardless of whether they already exist
        set_kwargs_attrs(self, kwargs, only_existing=True)
        
    def null_detection(self, time, gas_field):
        """
        Every detection method shares the same null_detection method defined here
        Inputs:
            time         a time object (Defined in simulation_classes)
        """
        rep_cond = (self.leaks.flux > 0) & self.leaks.reparable & (self.leaks.endtime <= time.current_time)
        self.repair_cost[time.time_index] += np.sum(self.leaks.repair_cost[rep_cond])
        self.end_intermittents(time)

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

    def survey_detection(self, time, gas_field, fcn_detection):
        oktime = self.ophrs['begin'] <= np.mod(time.current_time, 1) * 2400 < self.ophrs['end']
        timesize = time.delta_t * 2400 > self.ophrs['end'] - self.ophrs['begin']
        if np.mod(time.current_time, self.survey_interval) < time.delta_t:
            self.insurvey = True
        if self.insurvey and (oktime or timesize):
            remaining_comps = self.comps_per_timestep
            # only considering nonzero leaks is important because cleaning leaks that have been set to 0 from the leak
            # set happens periodically in the simulation and would otherwise cause indexing errors
            end_survey = False
            cond = []
            while remaining_comps > 0:
                site_name = self.find_site_name(gas_field, self.site_survey_index)
                if site_name == -1:
                    end_survey = True
                    break
                cond.extend(np.where((self.leaks.site_index == self.site_survey_index) &
                                     (self.leaks.comp_index >= self.comp_survey_index) &
                                     (self.leaks.comp_index < self.comp_survey_index + remaining_comps))[0])
                if remaining_comps + self.comp_survey_index > gas_field.sites[site_name]['parameters'].max_comp_ind:
                    remaining_comps -= (gas_field.sites[site_name]['parameters'].max_comp_ind - self.comp_survey_index)
                    self.comp_survey_index = 0
                    self.site_survey_index += 1
                else:
                    self.comp_survey_index += remaining_comps
                    remaining_comps = 0
            if len(cond) > 0:
                detect = fcn_detection(np.array(cond))
                # self.repair_cost[time.time_index] += \
                #     np.sum(np.random.choice(gas_field.repair_cost_dist.repair_costs,
                #                             np.sum(self.leaks.reparable[detect])))
                # Delete found leaks that are reparable
                self.leaks.endtime[detect[self.leaks.reparable[detect]]] = time.current_time + self.repair_delay
            if end_survey:
                self.insurvey = False
                self.site_survey_index = 0
                self.comp_survey_index = 0

    def end_intermittents(self, time):
        """
        Sets the emissions to zero for all leaks that have exceeded their duration
        :param time: a time object
        :return:
        """
        self.leaks.delete_leaks(np.where(self.leaks.endtime <= time.current_time)[0])

    def replacement_cap(self, time):
        """
        Enters the replacement value of the technology into tech.capital at the end of its lifetime.
        Inputs:
            time        Object defining time parameters for the simulation
            tech        LDAR technology object
        """
        lifetimes = 0
        while lifetimes < time.end_time:
            self.capital[max(1, round(lifetimes / time.delta_t))] = self.capital_0
            lifetimes += self.lifetime
