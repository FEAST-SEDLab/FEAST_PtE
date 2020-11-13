"""
A class for storing emission properties and functions for modifying emission proporeties throughout a simulation are
defined in this module.
"""
import pickle
import numpy as np
import pandas as pd


class Emission:
    """
    Stores all properties of all emissions that exist at a particular instant in a simulation.
    """
    def __init__(self, flux=(), reparable=True, site_index=(), comp_index=(), start_time=0, end_time=np.infty,
                 repair_cost=(), emission_id=None):
        """
        :param flux: An array of emission rates (array of floats--gram/second)
        :param reparable: An array of True/False values to indicate whether or not an emission is reparable
        :param site_index: An array indicating the index of the site that contains every emission
        :param comp_index: An array indicating the index of the component that is the source of each emission
        :param start_time: An array specifying the time when every emission begins
        :param end_time: An array specifying the time when every emission will end (days)
        :param emission_id:
        :param repair_cost: An array storing the cost of repairing every emission ($)
        """
        try:
            length_in = len(flux)
        except TypeError:
            length_in = 1

        if reparable is True:
            rep_array = np.ones(length_in, dtype=np.bool)
        elif reparable is False:
            rep_array = np.zeros(length_in, dtype=np.bool)
        else:
            rep_array = np.array(reparable)
        if emission_id is None:
            emission_id = np.linspace(0, length_in - 1, length_in, dtype=int)
        try:
            if len(end_time) == length_in:
                end_time = np.array(end_time)
        except TypeError:
            end_time = np.ones(length_in) * end_time
        self.emissions = pd.DataFrame({
            'flux': np.array(flux),
            'site_index': np.array(site_index),
            'comp_index': np.array(comp_index),
            'reparable': rep_array,
            'end_time': end_time,
            'repair_cost': np.array(repair_cost),
            'start_time': np.array(start_time)
        }, index=np.array(emission_id))
        self.emissions.index.name = 'emission_id'

    def get_current_emissions(self, time):
        """
        Returns all emissions that exist at time.current_time
        :param time: a Time object
        :return: a DataFrame of current emissions
        """
        cond = (self.emissions['start_time'] <= time.current_time) & (self.emissions['end_time'] > time.current_time)
        return self.emissions.loc[cond]

    def get_emissions_in_range(self, t0, t1, reparable=None):
        """
        Returns all emissions that existed between t0 and t1
        :param t0: beginning of interval (days)
        :param t1: end of interval (days)
        :param reparable: boolean condition. If set, only returns emissions with a matching reparable property
        :return: a DataFrame of all emissions that existed at any time in the interval t0:t1
        """
        cond = (self.emissions['start_time'] < t1) & (self.emissions['end_time'] >= t0)
        if reparable is not None:
            cond = cond & (self.emissions['reparable'] == reparable)
        return self.emissions[cond]

    def em_rate_in_range(self, t0, t1, reparable=None):
        """
        Returns the sum of emissions that existed between t0 and t1 integrated over the time period
        :param t0: beginning of interval (days)
        :param t1: end of interval (days)
        :param reparable: boolean condition. If set, only returns emissions with a matching reparable property
        :return: Average emission rate between t1 and t0 (g/s)
        """
        em = self.get_emissions_in_range(t0, t1, reparable=reparable)
        st = em['start_time'].to_numpy()
        st[st < t0] = t0
        et = em['end_time'].to_numpy()
        et[et > t1] = t1
        duration = et - st
        return np.sum(duration * em.flux) / (t1 - t0)

    def extend(self, *args):
        """
        Extends the existing emissions data frame with all of the entries in args
        :param args: a list of Emission objects
        :return:
        """
        emission_list = [self.emissions]
        [emission_list.append(a.emissions) for a in args]
        self.emissions = pd.concat(emission_list)


def bootstrap_emission_maker(n_em_in, comp_name, site, time, start_time=None, reparable=True):
    """
    Create leaks using a bootstrap method.

    :param n_em_in: number of leaks to generate
    :param comp_name: key to a Component object in site.comp_dict
    :param site: a Site object
    :param time: a Time object
    :param start_time: the times at which each emission begins
    :param reparable: Specifies whether emissions should be reparable or not (boolean)
    """
    if start_time is None:
        start_time = np.ones(n_em_in) * time.current_time
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
        n_leaks_key = leaks_per_well[counter] / sum(leaks_per_well) * n_em_in
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
    return Emission(flux=flux, reparable=reparable, start_time=start_time,
                    site_index=site_indexes, comp_index=comp_indexes, end_time=end_times, repair_cost=repair_costs)


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
    emission_objects_generator is a parent function that will be called to initialize gas fields

    :param dist_type: Type of leak distribution to be used
    :param leak_data_path: Path to a leak data file
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


def permitted_emission(n_emit, sizes, duration, time, site, comp_name, start_time):
    """
    Creates an emission object specifying new permitted emissions

    :param n_emit: number of emissions to create
    :param sizes: a list of leak sizes from which to specify the emission rate
    :param duration: a float defining the duration of the emission
    :param time: a Time object
    :param site: a Site object
    :param comp_name: Name of the component to be considered from within site.comp_dict
    :param start_times: array of times at which emissions start
    :return: an Emission object
    """
    flux = np.random.choice(sizes, n_emit)
    start_time = np.random.uniform(0, time.end_time, n_emit)
    reparable = False
    endtime = time.current_time + duration
    site_indexes = np.random.randint(site.site_inds[0], site.site_inds[1], len(flux))
    comp_indexes = comp_indexes_fcn(site, comp_name, len(flux))
    repair_cost = np.zeros(len(flux))
    return Emission(flux=flux, reparable=reparable, end_time=endtime,
                    site_index=site_indexes, comp_index=comp_indexes, repair_cost=repair_cost, start_time=start_time)
