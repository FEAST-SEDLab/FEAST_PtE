from . import EmissionSimModules
import pandas as pd
import numpy as np


def gascomp_reader(path_to_gas_comp, feast_delta_t, duration, comps_per_site, rep_cost_path,
                   met_data_path=None, met_start_hr=None):
    """
    Read a MEET 1.0 GasComp result file and return a feast Time object and GasField object specifying all of the
    emissions to occur in a FEAST simulation

    :param path_to_gas_comp: path to a gas_comp file
    :param feast_delta_t: time resolution of the FEAST simulation (days). FEAST will represent emissions with the \
    precision of the original MEET simulation, but will invoke LDAR events as though they occur instantaneously with \
    the time resolution specified here.
    :param duration: duration of the simulation to run in FEAST (days)
    :param comps_per_site: A dict of form {loc: number of components}. This will specify the number of components to \
    be inspected for leaks at every site. Must have a key for every location ID in the gas_comp file.
    :param rep_cost_path: Path to a repair cost data file. The FEAST will randomly assign these repair costs to leaks \
    simulated in MEET.
    :param met_data_path: Path to a meteorological data file (if left as None, FEAST will run without met data)
    :param met_start_hr: FEAST will rotate the meteorological data to begin at the hour specified by an integer here.
    :return time_obj: a FEAST Time object to specify a simulation
    :return gas_field: a FEAST GasField object
    """
    gas_comp, meet_delta_t = load_gas_comp_file(path_to_gas_comp)
    dat_dict = gas_comp_to_dict(gas_comp, meet_delta_t)
    time_obj, gas_field = gc_dat_to_gas_field(feast_delta_t, duration, comps_per_site, rep_cost_path,
                                              met_data_path,
                                              met_start_hr, **dat_dict)
    return time_obj, gas_field


def gc_dat_to_gas_field(feast_delta_t, duration, comps_per_site, rep_cost_path, met_data_path=None,
                        met_start_hr=None, size=0, loc=0, start=0, stop=0, emtype=0, em_id=None):
    """
    Ports data from a data_dict returned by the funcation gas_comp_to_dict into a FEAST GasField

    :param feast_delta_t: FEAST time resolution (days)
    :param duration: FEAST simulation duration (days)
    :param comps_per_site: A dict of form {loc: number of components}. This will specify the number of components to \
    be inspected for leaks at every site
    :param rep_cost_path: Path to a repair cost data file. The FEAST will randomly assign these repair costs to leaks \
    simulated in MEET
    :param met_data_path: Path to a meteorological data file (if left as None, FEAST will run without met data)
    :param met_start_hr: FEAST will rotate the meteorological data to begin at the hour specified by an integer here.
    :param size: an array of emission sizes (g/s)
    :param loc: an array of site ids
    :param start: an array of emission start times (days)
    :param stop: an array of emission stop times (days)
    :param emtype: an array of emission types (reparable or no)
    :param em_id: an array of emission IDs
    :return: a FEAST Time object
    :return: a FEAST GasField object
    """
    n_site = len(np.unique(loc))
    if len(comps_per_site) != n_site:
        raise ValueError("comps_per_site must be a dict or array like with len() == n_site ")
    if duration > np.max(stop):
        raise ValueError("Requested FEAST simulation duration is longer than the MEET simulation.")
    timeobj = EmissionSimModules.simulation_classes.Time(delta_t=feast_delta_t, end_time=duration)
    reparable_cond = emtype == 'FUGITIVE'
    site_dict = {}
    comp_ind = np.zeros(len(size))
    rep_costs_in = np.zeros(len(size))
    for loc_id in np.unique(loc):
        comp_general = EmissionSimModules.infrastructure_classes.Component(
            name='General ' + str(loc_id),
            repair_cost_path=rep_cost_path
        )
        basicpad = EmissionSimModules.infrastructure_classes.Site(
            name='basic pad ' + str(loc),
            comp_dict={
                # ------ The number of components is proportional to the number of wells, ind2
                'General ' + str(loc_id): {'number': comps_per_site[loc_id], 'parameters': comp_general},
            }
        )
        loc_inds = np.where(loc == loc_id)[0]
        comp_ind[loc_inds] = np.random.randint(0, comps_per_site[loc_id], len(loc_inds))
        site_dict['basic pad ' + str(loc_id)] = {'number': 1, 'parameters': basicpad}
        rep_costs_in[loc_inds] = np.random.choice(comp_general.repair_cost_dist.repair_costs, len(loc_inds))
    em = EmissionSimModules.emission_class_functions.Emission(
        flux=size,  # g/s
        reparable=reparable_cond,
        site_index=loc,
        end_time=stop,
        repair_cost=rep_costs_in,
        start_time=start,
        comp_index=comp_ind,
        emission_id=em_id
    )
    gas_field = EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj,
        emissions=em
    )
    if met_data_path:
        gas_field.met_data_path = met_data_path
        if met_start_hr:
            gas_field.met_data_maker(start_hr=met_start_hr)
        else:
            gas_field.met_data_maker(start_hr=np.random.randint(0, 8760))
    return timeobj, gas_field


def load_gas_comp_file(path_to_gas_comp):
    """
    Loads a gas_comp.csv file generated by MEET 1.0

    :param path_to_gas_comp: a path (str) to a gas_comp file
    :return: a DataFrame of the gas_comp file
    :return: the time resolution of meet simulation (seconds)
    """
    dat = pd.read_csv(path_to_gas_comp)
    cond = (dat.timestamp != 0) | ((dat.timestamp == 0) & (dat.METHANE > 0))
    dat = dat[cond]
    dat = dat.set_index(['name', 'emissionCategory'])
    dat = dat.sort_index()
    time_res_dict = {
        'SECONDS': 1,
        'MINUTES': 60,
        'HOURS': 3600,
        'DAYS': 3600 * 24
    }
    delta_t = time_res_dict[dat['emissionInterval'][0]]
    return dat, delta_t


def gas_comp_to_dict(gc, delta_t):
    """
    converts a gas_comp DataFrame to a dict specifying the start time, top time, emission rate, location,
    emission type and emission ID for every emitter in gas_comp. If an emission changes its emission rate,
    it is recorded as stopping and restarting at the time of the change. LDAR programs will treat the emission
    correctly due to the emission ID that persists after the change in emission rate.

    :param gc: a gas_comp DataFrame as created by the function load_gas_comp_file
    :param delta_t: The time resolution of the gas_comp file (seconds)
    :return: the dict containing emitter data from the gas_comp file
    """
    dat_dict = {
        'loc': [],
        'start': [],
        'stop': [],
        'size': [],
        'emtype': [],
        'em_id': []
    }
    em_counter = -1
    for emitter in gc.index.levels[0]:
        for em_type in gc.loc[(emitter)].index.unique(level='emissionCategory'):
            em_counter += 1
            tseries = gc.loc[(emitter, em_type)].timestamp
            tzero = ((tseries.iloc[1:] - tseries.iloc[:-1]) != delta_t).values
            methane_rate = gc.loc[(emitter, em_type), 'METHANE']
            tchange = ((methane_rate.iloc[1:] - methane_rate.iloc[:-1]) != 0).values
            tstops = list(np.sort(tseries[:-1][tzero | tchange]))
            tstarts = [tseries[0] - delta_t]
            tstarts.extend(list(np.sort(tseries[1:][tzero | tchange]) - delta_t))
            tstops.append(np.max(tseries))
            dat_dict['stop'].extend(tstops)
            dat_dict['start'].extend(tstarts)
            dat_dict['emtype'].extend([em_type] * len(tstarts))
            dat_dict['em_id'].extend(np.ones(len(tstarts)) * em_counter)
            dat_dict['loc'].extend(np.ones(len(tstarts)) * gc.loc[(emitter, em_type), 'location'][0])
            cond = gc.loc[(emitter, em_type), 'timestamp'].isin(np.array(tstarts) + delta_t)
            size_em = methane_rate[cond].values / delta_t * 1000  # g/s
            dat_dict['size'].extend(size_em)
    for k, v in dat_dict.items():
        dat_dict[k] = np.array(v)
    dat_dict['start'] = dat_dict['start'] / 3600 / 24  # days
    dat_dict['stop'] = dat_dict['stop'] / 3600 / 24  # days

    return dat_dict
