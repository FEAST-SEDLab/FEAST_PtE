import numpy as np
from pickle import load
from os import listdir
from os.path import isfile, join


def results_analysis(directory):
    """
    Process many realizations of a single scenario stored in a directory
    Inputs:
        directory       A directory of results files all generated under the same scenario
    Return:
        null_npv          array of null-NPV of each LDAR program in each realization [k$/well]
        emissions_timeseries  Array of emissions in each LDAR program in each realization at each time step
        costs                 Array of costs associated with each LDAR program (no discounting, all costs summed)
        techs           list of detection program names
    """
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    n_realizations = len(files)
    sample = load(open(directory + '/' + files[0], 'rb'))
    n_prog = len(sample.ldar_program_dict)
    # Initialize an array to store a time series of emissions for every LDAR program in every realization
    emissions_timeseries = np.zeros([n_prog, sample.time.n_timesteps, n_realizations])
    costs = np.zeros([n_prog, n_realizations])
    progs = list(sample.ldar_program_dict.keys())
    null_npv = dict()
    npv_keys = ['Repair', 'Finding', 'Gas', 'Total']
    for key in npv_keys:
        null_npv[key] = np.zeros([n_prog-1, n_realizations])
    # iterate through each realization
    for jindex in range(0, len(files)):
        sample = load(open(directory + '/' + files[jindex], 'rb'))
        # iterate through each LDAR program
        for index in range(0, len(progs)):
            emissions_timeseries[index, :, jindex] = sample.ldar_program_dict[progs[index]].emissions_timeseries
            costs[index, jindex] = np.sum(sample.ldar_program_dict[progs[index]].find_cost) + \
                np.sum(sample.ldar_program_dict[progs[index]].repair_cost)
        # iterate through each category of value
        null_npv_temp = npv_calculator(directory + '/' + files[jindex])
        for key in null_npv.keys():
            # no_repair_npv[key][:, jindex] = no_repair_npv_temp[key]
            null_npv[key][:, jindex] = null_npv_temp[key]
    return null_npv, emissions_timeseries, costs, progs


def npv_calculator(filepath):
    """
        Calculates the net present value (NPV) of each LDAR program in the results file
        Inputs:
            filepath        path to a results file
        Return:
            null_npv          NPV of each LDAR program compared to a scenario with only the Null LDAR program [k$/well]
    """
    sample = load(open(filepath, 'rb'))
    # site_count = sample.gas_field.site_count
    progs = list(sample.ldar_program_dict.keys())
    # The 'Null' LDAR program is special here because null_npv is calculated with respect to it
    if 'Null' not in progs:
        raise NameError('tech_dict must contain a "Null" detection method to use npv_calculator')
    null_emissions = np.array(sample.ldar_program_dict['Null'].emissions_timeseries)
    time = np.linspace(0, sample.time.n_timesteps * sample.time.delta_t, sample.time.n_timesteps)
    discount_array = (1 + sample.econ_settings.discount_rate)**(time / 365)
    # Initialize all arrays
    gas_value_n = np.zeros(len(progs) - 1)
    find_cost = np.zeros([len(progs), sample.time.n_timesteps])
    find_cost_null = np.zeros([len(progs) - 1, sample.time.n_timesteps])
    repair_cost_null = np.zeros([len(progs) - 1, sample.time.n_timesteps])
    # Store the cost of repairs in the null module
    null_repair_cost = sample.ldar_program_dict['Null'].repair_cost / discount_array
    # Calculate costs associated with each LDAR program
    null_correct = 0
    for index in range(0, len(progs)):
        find_cost[index, :] = sample.ldar_program_dict[progs[index]].find_cost / discount_array
        if progs[index] != 'Null':
            ind = index - null_correct
            repair_cost_null[ind, :] = sample.ldar_program_dict[progs[index]].repair_cost / discount_array - \
                null_repair_cost
            gas_value_n[ind] = sum((null_emissions - sample.ldar_program_dict[progs[index]].emissions_timeseries) /
                                   discount_array) * sample.time.delta_t * 24 * 3600 * sample.econ_settings.gas_price
            find_cost_null[ind, :] = sample.ldar_program_dict[progs[index]].find_cost / discount_array
        else:
            null_correct += 1
    # consolidate costs into totals
    find_n = np.sum(find_cost_null, 1)
    n_repair = np.sum(repair_cost_null, 1)
    null_npv = {'Repair': n_repair, 'Finding': find_n, 'Gas': gas_value_n}
    tot_n = gas_value_n - n_repair - find_n
    null_npv['Total'] = tot_n
    return null_npv
