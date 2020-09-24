"""
Leak data, leak distribution properties, and leak objects are created in this module
"""
import pickle
import numpy as np
from os.path import dirname, abspath
import os

# Constants:
g = 9.8  # g is the strength of gravity [m/s^2]
RHO_AIR = 1225  # density of air [g/m^3]
RHO_METHANE = 681  # density of methane at atmospheric pressure [g/m^3]


class Emission:
    """
        Stores a list of leaks
    """
    def __init__(self, flux=(), capacity=0,
                 reparable=True, site_index=(), comp_index=(), endtime=np.infty, repair_cost=()):
        """
        Inputs:
        flux                leak size (g/s)
        capacity            Expected total number of leaks to be stored in this instance of Leak (allows for faster
                            extend method)
        reparable           Indicates whether an emission is reparable leak or a permanent vent
        endtime             time when emission will stop if not repaired
        """
        try:
            length_in = len(flux)
        except TypeError:
            length_in = 1

        if capacity == 0:
            self.flux = np.array(flux)
            self.site_index = np.array(site_index)
            self.comp_index = np.array(comp_index)
            if reparable is True:
                self.reparable = np.ones(length_in, dtype=np.bool)
            elif reparable is False:
                self.reparable = np.zeros(length_in, dtype=np.bool)
            else:
                self.reparable = reparable
            try:
                if len(endtime) == length_in:
                    self.endtime = np.array(endtime)
            except TypeError:
                self.endtime = np.ones(length_in) * endtime
            self.repair_cost = np.array(repair_cost)
        else:
            self.flux = np.zeros(capacity)
            self.site_index = -np.ones(capacity, dtype=int)
            self.comp_index = -np.ones(capacity, dtype=int)
            self.reparable = np.ones(capacity, dtype=np.bool)
            self.endtime = np.ones(capacity) * np.inf
            self.repair_cost = np.zeros(capacity)
            self.flux[0:length_in] = flux
            self.reparable[0:length_in] = reparable
            self.site_index[0:length_in] = site_index
            self.comp_index[0:length_in] = comp_index
            self.endtime[0:length_in] = endtime
            self.repair_cost[0:length_in] = repair_cost
        self.n_leaks = length_in

    def extend(self, leak_obj_in):
        """
        Add a new leak
        Inputs:
            leak_obj_in     a Leak object
        """
        if len(self.flux) - leak_obj_in.n_leaks - self.n_leaks >= 0:
            self.flux[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.flux
            self.reparable[self.n_leaks: self.n_leaks+leak_obj_in.n_leaks] = leak_obj_in.reparable
            self.endtime[self.n_leaks: self.n_leaks+leak_obj_in.n_leaks] = leak_obj_in.endtime
            self.site_index[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.site_index
            self.comp_index[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.comp_index
            self.repair_cost[self.n_leaks: self.n_leaks + leak_obj_in.n_leaks] = leak_obj_in.repair_cost
        else:
            self.flux = np.append(self.flux[:self.n_leaks], leak_obj_in.flux)
            self.reparable = np.append(self.reparable[:self.n_leaks], leak_obj_in.reparable)
            self.endtime = np.append(self.endtime[:self.n_leaks], leak_obj_in.endtime)
            self.site_index = np.append(self.site_index[:self.n_leaks], leak_obj_in.site_index)
            self.comp_index = np.append(self.comp_index[:self.n_leaks], leak_obj_in.comp_index)
            self.repair_cost = np.append(self.repair_cost[:self.n_leaks], leak_obj_in.repair_cost)

        self.n_leaks += leak_obj_in.n_leaks

    def delete_leaks(self, indexes_to_delete):
        """
        Delete all parameters associated with leaks at indexes 'indexes_to_delete'
        indexes_to_delete       A list of leak indexes to delete, or the string 'all'
        """
        if type(indexes_to_delete) is str:
            if indexes_to_delete == 'all':
                indexes_to_delete = list(range(0, self.n_leaks))
            else:
                raise ValueError('indexes_to_delete must be a scalar, an array or the str "all"')
        self.flux[indexes_to_delete] = 0

    def clear_zeros(self):
        """
        This function removes all leaks from the leak object that have zero flux. This is in contrast to the
        delete_leaks method, which sets the flux to zero for specified leak indexes. These methods are implemented
        separately for computation efficiency: delete_leaks is called frequently, and only changes the flux value of
        a few indeces. clear_zeros is called less frequently and creates a new copy of the leak attribute arrays with th
        e 0 flux leaks omitted.
        :return:
        """
        indexes_to_delete = np.argwhere(self.flux[:self.n_leaks] == 0)
        self.flux = np.delete(self.flux, indexes_to_delete)
        self.reparable = np.delete(self.reparable, indexes_to_delete)
        self.endtime = np.delete(self.endtime, indexes_to_delete)
        self.site_index = np.delete(self.site_index, indexes_to_delete)
        self.comp_index = np.delete(self.comp_index, indexes_to_delete)
        self.repair_cost = np.delete(self.repair_cost, indexes_to_delete)
        try:
            self.n_leaks -= len(indexes_to_delete)
        except TypeError:
            self.n_leaks -= 1

    def sort_by_site(self):
        """
        sorts all of the leak attributes based on the site they are associated with.
        :return:
        """
        sortorder = np.argsort(self.site_index[:self.n_leaks])
        self.flux[:self.n_leaks] = self.flux[:self.n_leaks][sortorder]
        self.reparable[:self.n_leaks] = self.reparable[:self.n_leaks][sortorder]
        self.endtime[:self.n_leaks] = self.endtime[:self.n_leaks][sortorder]
        self.site_index[:self.n_leaks] = self.site_index[:self.n_leaks][sortorder]
        self.comp_index[:self.n_leaks] = self.comp_index[:self.n_leaks][sortorder]
        self.repair_cost[:self.n_leaks] = self.repair_cost[:self.n_leaks][sortorder]


def bootstrap_emission_maker(n_leaks_in, comp_name, site, time, capacity=0, reparable=True):
    """
    Create leaks using a bootstrap method
    n_leaks_in      number of leaks to generate
    comp_name       key to a Component object in site.comp_dict
    site            a Site object
    time            a time object
    capacity        the size of array to make to store leaks (if zero, the arrays are given a size equal to n_leaks_in).
    """
    comp = site.comp_dict[comp_name]['parameters']
    leak_params = comp.emission_params
    detection_methods = list(leak_params.leak_sizes.keys())
    flux = []
    round_err, leaks_per_well, counter = [], [], -1
    # Calculate leaks per well identified with each detection method stored in leak_params
    for method in detection_methods:
        n_leaks = len(leak_params.leak_sizes[method])
        n_wells = leak_params.well_counts[method]
        leaks_per_well.append(n_leaks/n_wells)
    # Generate the appropriate number of leaks from the distribution associated with each detection method
    for method in detection_methods:
        counter += 1
        n_leaks_key = leaks_per_well[counter] / sum(leaks_per_well) * n_leaks_in
        flux.extend(np.random.choice(leak_params.leak_sizes[method], int(n_leaks_key)))
        round_err.append(n_leaks_key % 1)
    # Add leaks omitted due to inability to add fractional leaks
    # The "round" function in the following line is intended to eliminate floating point errors.
    chooser = np.random.uniform(0, sum(round_err), round(sum(round_err)))
    error_intervals = np.cumsum(round_err)
    for choose in chooser:
        ind = 0
        # Add a leak from the appropriate detection method
        while choose > error_intervals[ind]:
            ind += 1
        flux.append(np.random.choice(leak_params.leak_sizes[detection_methods[ind]]))
    flux = np.array(flux)
    np.random.shuffle(flux)
    site_indexes = np.random.randint(site.site_inds[0], site.site_inds[1], len(flux))
    comp_indexes = comp_indexes_fcn(site, comp_name, len(flux))
    if site.comp_dict[comp_name]['parameters'].null_repair_rate > 0:
        end_times = time.current_time + \
                    np.random.exponential(1 / site.comp_dict[comp_name]['parameters'].null_repair_rate, len(flux))
    else:
        end_times = np.inf
    repair_costs = np.random.choice(comp.repair_cost_dist.repair_costs, len(flux))
    return Emission(flux=flux, reparable=reparable, capacity=capacity,
                    site_index=site_indexes, comp_index=comp_indexes, endtime=end_times, repair_cost=repair_costs)


def comp_indexes_fcn(site, comp_name, n_inds):
    """
    Returns an array of indexes to associate with new emissions
    :param site: a EmissionSimModules.simulation_classes.Site object
    :param comp_name: name of a component contained in Site.comp_dict
    :param n_inds: Integer of indexes to generate
    :return: An array of indexes in the range specified for the relevant component
    """
    low_ind = site.comp_dict[comp_name]['comp_indexes'][0]
    high_ind = site.comp_dict[comp_name]['comp_indexes'][1]
    return np.random.randint(low_ind, high_ind, n_inds)


def emission_objects_generator(dist_type, emission_data_path, custom_emission_maker=None):
    """
    leak_objects is a parent function that will be called to initialize gas fields
    Inputs:
        dist_type           Type of leak distribution to be used
        leak_data_path      Path to a leak data file
    """
    rsc_path = emission_data_path
    with open(rsc_path, 'rb') as f:
        emission_params = pickle.load(f)

    # Define leak params and emission_size_maker based on the leak distribution type
    if dist_type == 'bootstrap':
        emission_size_maker = bootstrap_emission_maker
    elif dist_type.lower() == 'custom':
        if custom_emission_maker is None:
            raise ValueError("custom_emission_maker must be defined for a custom emission distribution type")
        emission_size_maker = custom_emission_maker
    else:
        raise NameError('emission distribution type unsupported in GasField')

    # Number of leaking components at each well (Poisson distribution)
    detection_types, em_per_well, em_per_comp = emission_params.leak_sizes.keys(), 0, 0
    # This sums the leaks found per well by every method in emission_data
    for key in detection_types:
        em_per_well += len(emission_params.leak_sizes[key]) / emission_params.well_counts[key]
        em_per_comp += len(emission_params.leak_sizes[key]) / \
            (emission_params.comp_counts[key] * emission_params.well_counts[key])
    return emission_size_maker, emission_params, em_per_well, em_per_comp


def permitted_emission(n_emit, sizes, duration, time, site, comp_name):
    """
    Creates a leak object specifying new permitted emissions
    :param n_emit: number of emissions to create
    :param sizes: a list of leak sizes from which to specify the emission rate
    :param duration: a float defining the duration of the emission
    :param time: a Time object
    :param site: a Site object
    :param comp_name: Name of the component to be considered from within site.comp_dict
    :return: A Leak object
    """
    flux = np.random.choice(sizes, n_emit)
    reparable = False
    endtime = time.current_time + duration
    site_indexes = np.random.randint(site.site_inds[0], site.site_inds[1], len(flux))
    comp_indexes = comp_indexes_fcn(site, comp_name, len(flux))
    repair_cost = np.zeros(len(flux))
    return Emission(flux=flux, reparable=reparable, endtime=endtime,
                    site_index=site_indexes, comp_index=comp_indexes, repair_cost=repair_cost)
