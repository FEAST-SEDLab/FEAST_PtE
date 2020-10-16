"""
This module stores component, gasfield and site classes to represent infrastructure in a simulation
"""

import pickle
import numpy as np
import pandas as pd

from feast.EmissionSimModules import emission_class_functions as ecf
from feast.EmissionSimModules.emission_class_functions import emission_objects_generator as leak_obj_gen
from feast.EmissionSimModules.simulation_classes import Time


class Component:
    """
    A class to store parameters defining a component (for example, name, leak production rate, leak size
    distribution, etc)
    """
    def __init__(self, repair_cost_path=None, emission_data_path=None, base_reparable=None, custom_emission_maker=None,
                 emission_production_rate=0, emission_per_comp=None, episodic_emission_sizes=[0],
                 episodic_emission_per_day=0, episodic_emission_duration=0, vent_sizes=[0],
                 vent_period=np.infty, vent_starts=np.array([]), vent_duration=0, name='default',
                 null_repair_rate=None, dist_type='bootstrap'):
        """
        :param repair_cost_path: path to a repair cost data file
        :param emission_data_path: path to an emission data file
        :param base_reparable: Defines whether emissions generated are reparable with a boolean true/false
        :param custom_emission_maker: Optional custom defined function for creating new emissions
        :param emission_production_rate: The rate at which new emissions are created (emissions per day per component)
        :param emission_per_comp: The number of emissions expected per component (must be less than 1)
            If emission_per_comp is left as None, then emission_per_comp is set equal to the emissions per component
            recorded in the file at emission_data_path.
        :param episodic_emission_sizes: A list of emission sizes to draw from for episodic emissions (g/s)
        :param episodic_emission_per_day: The average frequency at which episodic emissions occur (1/days)
        :param episodic_emission_duration: The duration of episodic emissions (days)
        :param vent_sizes: A list of emission sizes for periodic emissions (g/s)
        :param vent_period: The time between emissions (days)
        :param vent_duration: the time that a periodic vent persits (days)
        :param vent_starts: the time at which the first periodic vent occurs at each component in the simulation
        :param name: A name for the instance of Component
        :param null_repair_rate: the rate at which fugitive emissions are repaired. If None, a steady state
            assumption is enforced based on emission_production_rate and emission_per_comp.
        :param dist_type: The type of distribution to be used in determining emission rates for new emissions
        """
        self.name = name
        self.repair_cost_path = repair_cost_path
        self.dist_type = dist_type
        self.emission_data_path = emission_data_path
        self.emission_production_rate = emission_production_rate
        self.emission_per_comp = emission_per_comp
        self.custom_emission_maker = custom_emission_maker
        self.base_reparable = base_reparable
        self.episodic_emission_sizes = episodic_emission_sizes
        self.episodic_emission_per_day = episodic_emission_per_day
        self.episodic_emission_duration = episodic_emission_duration
        self.vent_sizes = vent_sizes
        self.vent_period = vent_period
        self.vent_duration = vent_duration
        self.vent_starts = vent_starts
        # ---- Distribution of leak repair costs
        if self.repair_cost_path:
            rsc_path = self.repair_cost_path
            with open(rsc_path, 'rb') as f:
                self.repair_cost_dist = pickle.load(f)
        if self.emission_data_path:
            self.emission_size_maker, self.emission_params, self.emission_per_well, emission_per_comp \
                = leak_obj_gen(self.dist_type, self.emission_data_path, self.custom_emission_maker)
            if self.emission_per_comp is None:
                self.emission_per_comp = emission_per_comp
        if null_repair_rate is None:
            if self.emission_per_comp == 0 or self.emission_production_rate==0:
                self.null_repair_rate = 0
            else:
                self.null_repair_rate = self.emission_production_rate / self.emission_per_comp
        else:
            self.null_repair_rate = null_repair_rate
        self.intermittent_emission_maker = ecf.permitted_emission


class GasField:
    """
    GasField accommodates all data that defines a gas field at the beginning of a simulation.
    """
    def __init__(self, time=None, sites=None, initial_emissions=None, new_emissions=None,
                 start_times=None, met_data_path=None):
        """
        :param time: A FEAST time object
        :param sites: a dict of sites like this: {'name': {'number': n_sites, 'parameters': site_object}}
        :param initial_emissions: A FEAST emission object defining emissions at the beginning of the simulation
        :param new_emissions: A list of FEAST emission objects with length time.n_timesteps
        :param start_times: The times at which new emissions begin
        :param met_data_path: A path to a met data file
        """
        self.sites = sites
        self.initial_emissions = initial_emissions
        self.new_emissions = new_emissions
        self.start_times = start_times
        self.met_data_path = met_data_path

        # -------------- Calculated parameters --------------
        self.n_comps, self.n_sites = 0, 0
        self.comp_dict = {}
        # dict to store met data
        self.met = {}

        # Define indexes and initialize emissions
        self.set_indexes()
        if self.met_data_path:
            self.met_data_maker()
        if self.initial_emissions is None:
            self.initialize_emissions(time)
        if self.new_emissions is None:
            self.emerging_emissions(time)
        self.input_emissions = []
        for ind in range(time.n_timesteps):
            cond = np.where(self.start_times == ind)[0]
            self.input_emissions.append(ecf.Emission(flux=self.new_emissions.flux[cond],
                                                     reparable=self.new_emissions.reparable[cond],
                                                     site_index=self.new_emissions.site_index[cond],
                                                     comp_index=self.new_emissions.comp_index[cond],
                                                     endtime=self.new_emissions.endtime[cond],
                                                     repair_cost=self.new_emissions.repair_cost[cond]))

    def initialize_emissions(self, time):
        """
        Create emissions that exist at the beginning of the simulation

        :param time:
        :return:
        """
        cap_est = 0
        for site_dict in self.sites.values():
            cap_est += site_dict['number'] * 10
        self.initial_emissions = ecf.Emission(capacity=cap_est)
        # This generates new leaks for each component type in each site type
        for sitedict in self.sites.values():
            site = sitedict['parameters']
            for comp_name in site.comp_dict:
                compobj = site.comp_dict[comp_name]['parameters']
                n_comp = sitedict['number'] * site.comp_dict[comp_name]['number']
                if compobj.emission_production_rate > 0:
                    n_leaks = np.random.binomial(n_comp, compobj.emission_per_comp)
                else:
                    n_leaks = 0
                self.emission_maker(n_leaks, self.initial_emissions, comp_name, n_comp, time, site)

    def set_indexes(self):
        """
        Counts components for each site and assigns appropriate indexes
        """
        site_ind = 0
        for site_dict in self.sites.values():
            site = site_dict['parameters']
            self.n_sites += site_dict['number']
            # This ensures that site indexes do not overlap between site types.
            site.site_inds = [site_ind, site_ind + site_dict['number']]
            if site.prod_dat is not None:
                site.production = np.random.choice(site.prod_dat, site_dict['number'])
            for compname, comp_d in site.comp_dict.items():
                comp = comp_d['parameters']
                self.n_comps += comp_d['number'] * site_dict['number']
                if compname in self.comp_dict:
                    raise ValueError("All component names must be unique.")
                self.comp_dict[compname] = comp
            site_ind += site_dict['number']
        self.n_sites = site_ind

    def emerging_emissions(self, time):
        """
        Defines emissions that emerge during a simulation
        :param time:
        :return:
        """
        self.new_emissions = ecf.Emission()
        for site_dict in self.sites.values():
            site = site_dict['parameters']
            for compname, comp in site.comp_dict.items():
                if comp['parameters'].vent_duration > 0:
                    comp['parameters'].vent_starts = np.random.uniform(0, comp['parameters'].vent_period,
                                                                       comp['number'])
                n_comp = site_dict['number'] * comp['number']
                n_leaks = np.random.poisson(n_comp * comp['parameters'].emission_production_rate * time.end_time)
                n_episodic = np.random.poisson(n_comp * comp['parameters'].episodic_emission_per_day * time.end_time)
                self.emission_maker(n_leaks, self.new_emissions, compname, n_comp, time, site, n_episodic=n_episodic)
        self.start_times = np.random.randint(0, time.n_timesteps, self.new_emissions.n_em, dtype=int)
        self.new_emissions.endtime += self.start_times * time.delta_t

    def emission_size_maker(self, time):
        """
        Creates a new set of leaks based on attributes of the gas field
        :param time: a time object (the parameter delta_t is used)
        :return new_leaks: the new leak object
        """
        new_leaks = ecf.Emission()
        for site_dict in self.sites.values():
            site = site_dict['parameters']
            for compname, comp in site.comp_dict.items():
                n_comp = site_dict['number'] * comp['number']
                n_leaks = np.random.poisson(n_comp * comp['parameters'].emission_production_rate * time.delta_t)
                self.emission_maker(n_leaks, new_leaks, compname, n_comp, time, site)
        return new_leaks

    def met_data_maker(self):
        if self.met_data_path:
            met_dat = pd.read_csv(self.met_data_path, header=1)
            self.met['wind speed'] = met_dat['wind speed (m/s)']
            self.met['wind direction'] = met_dat['wind direction (degrees clockwise from North)']
            self.met['temperature'] = met_dat['temperature (Celsius)']
            self.met['relative humidity'] = met_dat['relative humidity (%)']
            self.met['precipitation'] = met_dat['precipitation (mm)']
            self.met['albedo'] = met_dat['albedo (-)']
            self.met['ceiling height'] = met_dat['ceiling height (m)']
            self.met['cloud cover'] = met_dat['cloud cover (%)']
            self.met['solar intensity'] = met_dat['solar intensity (direct normal irradiance--W/m^2)']

    def get_met(self, time, parameter_names, interp_modes='mean', ophrs=None):
        """
        Return the relevant meteorological condition, accounting for discrepancies between simulation time resolution
        and data time resolution
        
        :param time: time object
        :param parameter_names: specify a list of meteorological conditions to return
        :param interp_modes: can be a list of strings: mean, median, max or min
        :param ophrs: Hours to consider when interpolating met data should be of form {'begin': 5, 'end':17}
        :return met_conds: dict of meteorological conditions
        """
        if ophrs is None:
            ophrs = {'begin': 0, 'end': 24}
        hour_index = int(np.mod(time.current_time * 24, 8760))
        # if a string is passed, put it in a list with one entry
        if type(parameter_names) is str:
            parameter_names = [parameter_names]
        if type(interp_modes) is str:
            interp_modes = [interp_modes for _ in range(len(parameter_names))]
        met_conds = {}
        for ind in range(len(parameter_names)):
            parameter_name = parameter_names[ind].lower()
            interp_mode = interp_modes[ind]
            if time.delta_t <= 1/24:
                met_conds[parameter_name] = self.met[parameter_name][hour_index]
            else:
                hr = np.mod(hour_index, 24)
                start_index = hour_index - hr + int(np.max([hr, ophrs['begin']]))
                end_index = hour_index - hr + int(np.min([hr + time.delta_t * 24, hour_index + ophrs['end']]))
                relevant_metdat = self.met[parameter_name][start_index:end_index]
                if interp_mode.lower() == 'mean':
                    met_conds[parameter_name] = np.mean(relevant_metdat)
                elif interp_mode.lower() == 'max':
                    met_conds[parameter_name] = np.max(relevant_metdat)
                elif interp_mode.lower() == 'min':
                    met_conds[parameter_name] = np.min(relevant_metdat)
                elif interp_mode.lower() == 'median':
                    met_conds[parameter_name] = np.median(relevant_metdat)
                else:
                    raise ValueError("Invalid meteorological data type.")
        return met_conds

    @staticmethod
    def emission_maker(n_leaks, new_leaks, comp_name, n_comp, time, site, n_episodic=None):
        """
        Updates an Emission object with new values returned by emission_size_maker

        :param n_leaks: number of new leaks to create
        :param new_leaks: a leak object to extend
        :param comp_name: name of a component object included in site.comp_dict
        :param n_comp: the number of components to model
        :param time: a time object
        :param site: a site object
        :param n_episodic: number of episodic emissions to create
        :return: None
        """
        comp = site.comp_dict[comp_name]['parameters']
        if n_leaks > 0:
            new_leaks.extend(comp.emission_size_maker(n_leaks, comp_name, site, time, reparable=comp.base_reparable))
        if n_episodic is None:
            n_episodic = np.random.poisson(n_comp * comp.episodic_emission_per_day * time.delta_t)
        new_leaks.extend(comp.intermittent_emission_maker(n_episodic,
                                                          comp.episodic_emission_sizes,
                                                          comp.episodic_emission_duration,
                                                          time, site, comp_name))
        n_vent = 0
        if comp.vent_starts.size > 0:
            n_vent = np.random.poisson(comp.vent_duration / comp.vent_period * n_comp *
                                       min(1, time.delta_t/comp.vent_period))
        new_leaks.extend(comp.intermittent_emission_maker(n_vent, comp.vent_sizes, comp.vent_duration,
                                                          time, site, comp_name))
        return None


class Site:
    """
    A class to store the number and type of components associated with a site.
    """
    def __init__(self, name='default', comp_dict=None, prod_dat=None):
        """
        :param name: The name of the site object (a string)
        :param comp_dict: A dict of components at the site, for example:
            {'name': {'number': 650, 'parameters': Component()}}
        :param prod_dat:
        """
        self.name = name  # A string
        self.comp_dict = comp_dict
        self.prod_dat = prod_dat

        # Calculated parameters
        comp_ind = 0
        # Note: This loop ensures that comp indexes for different comp types in the same site do not overlap
        for comp_name in self.comp_dict:
            self.comp_dict[comp_name]['comp_indexes'] = [comp_ind, self.comp_dict[comp_name]['number'] + comp_ind]
            comp_ind += self.comp_dict[comp_name]['number']
        self.max_comp_ind = comp_ind
