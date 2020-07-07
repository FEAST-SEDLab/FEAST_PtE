"""
    feast_classes stores the basic classes used by FEAST.
    Additional classes are stored in DetectionModules directory and leak_objects.
"""
import numpy as np
import pandas as pd
from .emission_class_functions import leak_objects_generator as leak_obj_gen
from . import emission_class_functions as lcf
import pickle
from .simulation_functions import set_kwargs_attrs
from os.path import dirname, abspath
import os


class Time:
    """
    Instances of the time class store all time related information during a simulation
    """
    def __init__(self, delta_t=3650/4000, end_time=10 * 365, current_time=0):
        """
        Inputs:
        delta_t                 length of one timestep (days)
        end_time                length of the simulation (days)
        current_time            current time in a simulation (days)
        """
        self.n_timesteps = int(end_time/delta_t+1)
        self.end_time = end_time
        self.delta_t = delta_t
        self.current_time = current_time
        self.time_index = 0


class Component:
    """
    A class to store parameters defining a component (for example, name, leak production rate, leak size
    distribution, etc)
    """
    def __init__(self, null_repair_rate=None, **kwargs):
        self.name = 'default'
        self.component_type = "fugitive leak source"
        # Repair cost data file path
        self.repair_cost_path = 'fernandez_leak_repair_costs_2006.p'
        # Events with a constant probability of occurring and duration determined by the null repair rate
        # emission specifications may or may not be reparable before the null repair process
        self.dist_type = 'bootstrap'
        self.emission_data_path = 'production_emissions.p'
        self.emission_production_rate = 1e-5  # new emissions per component per day
        self.emission_per_comp = 'default'
        self.custom_leak_maker = None
        self.base_reparable = True
        # Permitted events with a constant probability of occurring and known duration NOT REPARABLE
        self.episodic_emission_sizes = [0]  # g/s
        self.episodic_emission_per_day = 0
        self.episodic_emission_duration = 0
        # Periodic vents: repetitive events with a known period between emissions and known duration NOT REPARABLE
        self.vent_sizes = [0]  # g/s
        self.vent_period = np.infty  # days
        self.vent_duration = 0  # days
        self.vent_starts = np.array([])
        # Update any attributes defined by kwargs
        set_kwargs_attrs(self, kwargs)
        # Distribution of leak repair costs
        rsc_path, _ = os.path.split(dirname(abspath(__file__)))
        if self.repair_cost_path in os.listdir(os.path.join(rsc_path, 'InputData', 'DataObjectInstances')):
            rsc_path = os.path.join(rsc_path, 'InputData', 'DataObjectInstances', self.repair_cost_path)
        else:
            rsc_path = self.repair_cost_path
        with open(rsc_path, 'rb') as f:
            self.repair_cost_dist = pickle.load(f)
        self.leak_size_maker, self.leak_params, self.emission_per_well, emission_per_comp \
            = leak_obj_gen(self.dist_type, self.emission_data_path, self.custom_leak_maker)
        if self.emission_per_comp == 'default':
            self.emission_per_comp = emission_per_comp
        if null_repair_rate is None:
            if self.emission_per_comp == 0:
                self.null_repair_rate = 0
            else:
                self.null_repair_rate = self.emission_production_rate / self.emission_per_comp
        else:
            self.null_repair_rate = null_repair_rate
        self.emission_maker = lcf.permitted_emission


class GasField:
    """
    GasField accommodates all data that defines a gas field at the beginning of a simulation.
    """
    def __init__(self, time=Time(), **kwargs):
        """
        Input params:
            initial_leaks     The set of leaks that exist at the beginning of the simulation
            null_repair_rate  The rate at which leaks are repaired in the Null process (repairs/leak/day)
            kwargs           All attributes defined in the kwargs section below
        """

        # -------------- Attributes that can be defined with kwargs --------------
        # Number of wells to be simulated
        self.sites = {
            'default': {'number': 100, 'parameters': Site()}
        }

        # Maximum number of wells to be surveyed with a single capital investment
        self.max_count = 6000
        # Driving distance between wells
        self.site_spacing = 700  # m
        # Initial leaks defined for the gas field
        self.initial_emissions = None
        # emissions to be created during the simulation
        self.new_leaks = None
        # emissions to be created during the simulation
        self.start_times = None
        # path to TMY data
        self.met_data_path = None
        # Update any attributes defined by kwargs
        set_kwargs_attrs(self, kwargs)

        # -------------- Calculated parameters --------------
        self.n_comps, self.n_sites = 0, 0
        site_ind = 0
        self.avg_vent = 0
        self.comp_dict = {}
        # dict to store met data
        self.met = {}

        # Define met data
        if self.met_data_path:
            self.met_data_maker(time)

        # This loop counts components for each site
        for site_dict in self.sites.values():
            site = site_dict['parameters']
            self.n_sites += site_dict['number']
            # This ensures that site indexes do not overlap between site types.
            site.site_inds = [site_ind, site_ind + site_dict['number']]
            site.production = np.random.choice(site.prod_dat, site_dict['number'])
            for compname, comp_d in site.comp_dict.items():
                comp = comp_d['parameters']
                comp_count = comp_d['number']
                if comp.vent_duration > 0:
                    comp.vent_starts = np.random.uniform(0, comp.vent_period, comp_count)
                    self.avg_vent += comp.vent_duration / comp.vent_period * comp_count * site_dict['number']
                self.n_comps += comp_d['number'] * site_dict['number']
                if compname in self.comp_dict:
                    raise ValueError("All component names must be unique.")
                self.comp_dict[compname] = comp
            site_ind += site_dict['number']
        # create initial_leaks
        cap_est = 0
        for site_dict in self.sites.values():
            cap_est += site_dict['number'] * 10
        if self.initial_emissions is None:
            self.initial_emissions = lcf.Emission(capacity=cap_est)
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
                    self.leak_maker(n_leaks, self.initial_emissions, comp_name, n_comp, time, site)
                if site.site_em_dist:
                    self.site_emissions_enforcer(site)
        self.n_sites = site_ind
        if self.new_leaks is None:
            self.new_leaks = lcf.Emission()
            for site_dict in self.sites.values():
                site = site_dict['parameters']
                for compname, comp in site.comp_dict.items():
                    n_comp = site_dict['number'] * comp['number']
                    n_leaks = np.random.poisson(n_comp * comp['parameters'].emission_production_rate * time.end_time)
                    n_episodic = np.random.poisson(n_comp * comp['parameters'].episodic_emission_per_day * time.end_time)
                    self.leak_maker(n_leaks, self.new_leaks, compname, n_comp, time, site, n_episodic=n_episodic)
            self.start_times = np.random.randint(0, time.n_timesteps, self.new_leaks.n_leaks, dtype=int)
            self.new_leaks.endtime += self.start_times * time.delta_t
        input_leaks = []
        for ind in range(time.n_timesteps):
            cond = np.where(self.start_times == ind)[0]
            input_leaks.append(lcf.Emission(flux=self.new_leaks.flux[cond],
                                            reparable=self.new_leaks.reparable[cond],
                                            site_index=self.new_leaks.site_index[cond],
                                            comp_index=self.new_leaks.comp_index[cond],
                                            endtime=self.new_leaks.endtime[cond],
                                            repair_cost=self.new_leaks.repair_cost[cond]))
        self.input_leaks = input_leaks

    # Define functions and parameters related to leaks
    def leak_size_maker(self, time):
        """
        Creates a new set of leaks based on attributes of the gas field
        :param time: a time object (the parameter delta_t is used)
        :return new_leaks: the new leak object
        """
        new_leaks = lcf.Emission()
        for site_dict in self.sites.values():
            site = site_dict['parameters']
            for compname, comp in site.comp_dict.items():
                n_comp = site_dict['number'] * comp['number']
                n_leaks = np.random.poisson(n_comp * comp['parameters'].emission_production_rate * time.delta_t)
                self.leak_maker(n_leaks, new_leaks, compname, n_comp, time, site)
        return new_leaks

    def met_data_maker(self, time):
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

    def get_met(self, time, parameter_names, interp_modes='mean', op_hrs={'begin': 0, 'end': 24}):
        """
        Return the relevant meteorological condition, accounting for discrepancies between simulation time resolution
        and data time resolution
        :param time: time object
        :param parameter_names: specify a list of meteorological conditions to return
        :param interp_mode: can be a list of strings: mean, median, max or min
        :return met_conds: dict of meteorological coditions
        """
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
                start_index = hour_index - hr + int(np.max([hr, op_hrs['begin']]))
                end_index = hour_index - hr + int(np.min([hr + time.delta_t * 24, hour_index + op_hrs['end']]))
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
    def leak_maker(n_leaks, new_leaks, comp_name, n_comp, time, site, n_episodic=None):
        """
        Updates a leak object with new values returned by leak_size_maker
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
        new_leaks.extend(comp.leak_size_maker(n_leaks, comp_name, site, time, reparable=comp.base_reparable))
        if n_episodic is None:
            n_episodic = np.random.poisson(n_comp * comp.episodic_emission_per_day * time.delta_t)
        new_leaks.extend(comp.emission_maker(n_episodic,
                                             comp.episodic_emission_sizes,
                                             comp.episodic_emission_duration,
                                             time, site, comp_name))
        n_vent = 0
        if comp.vent_starts.size > 0:
            n_vent = np.random.poisson(comp.vent_duration / comp.vent_period * n_comp *
                                       min(1, time.delta_t/comp.vent_period))
        new_leaks.extend(comp.emission_maker(n_vent, comp.vent_sizes, comp.vent_duration, time, site, comp_name))
        return None


# FinanceSettings stores all parameters relating to economic calculations
class FinanceSettings:
    def __init__(self, gas_price=2E-4, discount_rate=0.08):
        self.gas_price = gas_price  # dollars/gram (2e-4 $/g=$5/mcf methane at STP)
        self.discount_rate = discount_rate


class Results:
    """
    Class in which to save results
    """
    def __init__(self, time, gas_field, ldar_program_dict, econ_set):
        """
        Inputs:
        time                    Time object
        gas_field               GasField object
        ldar_program_dict       dict of detection methods and associated data
        econ_set                Economic settings defined for the simulation
        """
        self.time = time
        self.gas_field = gas_field
        self.ldar_program_dict = ldar_program_dict
        self.econ_settings = econ_set


class Site:
    """
    A class to store the number and type of components associated with a site.
    """
    def __init__(self, prod_dat=None, **kwargs):
        self.name = 'default'
        self.comp_dict = {'default': {'number': 650, 'parameters': Component()}}  # Based on Fort Worth data
        self.site_em_dist = False
        if prod_dat is None:
            # By default, load the production distribution from Colorado
            rsc_path, _ = os.path.split(dirname(abspath(__file__)))
            rsc_path = os.path.join(rsc_path, 'InputData', 'DataObjectInstances', 'COGCC_site_prod_2019.p')
            with open(rsc_path, 'rb') as f:
                # Column 2 contains gas production data, and that is all that is saved by default
                self.prod_dat = pickle.load(f).site_prod[:, 2]
            # By default, only use sites that recorded nonzero production
            self.prod_dat = self.prod_dat[self.prod_dat > 0]
        else:
            self.prod_dat = prod_dat
        set_kwargs_attrs(self, kwargs)

        # Calculated parameters
        comp_ind = 0
        # Note: This loop ensures that comp indexes for different comp types in the same site do not overlap
        for comp_name in self.comp_dict:
            self.comp_dict[comp_name]['comp_indexes'] = [comp_ind, self.comp_dict[comp_name]['number'] + comp_ind]
            comp_ind += self.comp_dict[comp_name]['number']
        self.max_comp_ind = comp_ind
