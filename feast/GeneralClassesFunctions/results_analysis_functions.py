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
        no_repair_npv     array of no-repair-NPVs of each LDAR program in each realization [k$/well]
        null_npv          array of null-NPV of each LDAR program in each realization [k$/well]
        leaks_found       Cumulative lists of leaks found by each LDAR program over all realizations
        emissions_timeseries  Array of emissions in each LDAR program in each realization at each time step
    """
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    n_realizations = len(files)
    sample = load(open(directory + '/' + files[0], 'rb'))
    n_techs = len(sample.tech_dict)
    # Initialize an array to store a time series of emissions for every LDAR program in every realization
    emissions_timeseries = np.zeros([n_techs + 1, sample.time.n_timesteps, n_realizations])
    leaks_found = [[] for _ in range(n_techs)]
    techs = list(sample.tech_dict.keys())
    no_repair_npv = dict()
    null_npv = dict()
    npv_keys = ['Capital', 'Maintenance', 'Repair', 'Finding', 'Gas', 'Total']
    for key in npv_keys:
        no_repair_npv[key] = np.zeros([n_techs, n_realizations])
        null_npv[key] = np.zeros([n_techs-1, n_realizations])
    # iterate through each realization
    for jindex in range(0, len(files)):
        sample = load(open(directory + '/' + files[jindex], 'rb'))
        # iterate through each LDAR program
        for index in range(0, len(techs)):
            emissions_timeseries[index, :, jindex] = sample.tech_dict[techs[index]].emissions
            if techs[index] == 'Null':
                leaks_found[index].append(sample.tech_dict[techs[index]].null_repaired)
            else:
                leaks_found[index].append(sample.tech_dict[techs[index]].leaks_found)
        # iterate through each category of value
        no_repair_npv_temp, null_npv_temp = npv_calculator(directory + '/' + files[jindex])
        for key in null_npv.keys():
            # no_repair_npv[key][:, jindex] = no_repair_npv_temp[key]
            null_npv[key][:, jindex] = null_npv_temp[key]
    return no_repair_npv, null_npv, leaks_found, emissions_timeseries


def npv_calculator(filepath):
    """
        Calculates the net present value (NPV) of each LDAR program in the results file
        Inputs:
            filepath        path to a results file
        Return:
            no_repair_npv     NPV of each LDAR program compared to a scenario with no repair [k$/well]
            null_npv          NPV of each LDAR program compared to a scenario with only the Null LDAR program [k$/well]
    """
    sample = load(open(filepath, 'rb'))
    # site_count = sample.gas_field.site_count
    techs = list(sample.tech_dict.keys())
    # The 'Null' LDAR program is special here because null_npv is calculated with respect to it
    if 'Null' not in techs:
        raise NameError('tech_dict must contain a "Null" detection method to use npv_calculator')
    null_emissions = np.array(sample.tech_dict['Null'].emissions)
    time = np.linspace(0, sample.time.n_timesteps * sample.time.delta_t, sample.time.n_timesteps)
    discount_array = (1 + sample.econ_settings.discount_rate)**(time / 365)
    # Initialize all arrays
    gas_value_n = np.zeros(len(techs)-1)
    gas_value_no_repair = np.zeros(len(techs))
    find_cost = np.zeros([len(techs), sample.time.n_timesteps])
    find_cost_null = np.zeros([len(techs)-1, sample.time.n_timesteps])
    repair_cost_null = np.zeros([len(techs)-1, sample.time.n_timesteps])
    repair_cost_no_repair = np.zeros([len(techs), sample.time.n_timesteps])
    capital_cost = np.zeros([len(techs), sample.time.n_timesteps])
    capital_cost_null = np.zeros([len(techs)-1, sample.time.n_timesteps])
    maintenance_cost = np.zeros([len(techs), sample.time.n_timesteps])
    maintenance_cost_null = np.zeros([len(techs)-1, sample.time.n_timesteps])
    # Store the cost of repairs in the null module
    null_repair_cost = sample.tech_dict['Null'].repair_cost / discount_array
    # Calculate costs associated with each LDAR program
    null_correct = 0
    for index in range(0, len(techs)):
        find_cost[index, :] = sample.tech_dict[techs[index]].find_cost / discount_array
        maintenance_cost[index, :] = sample.tech_dict[techs[index]].maintenance / discount_array
        capital_cost[index, :] = sample.tech_dict[techs[index]].capital / discount_array
        # gas_value_no_repair[index] = sum((np.array(sample.no_repair_emissions)-sample.tech_dict[techs[index]].emissions) /
        #                                  discount_array) * sample.time.delta_t*24*3600*sample.econ_settings.gas_price
        # repair_cost_no_repair[index, :] = sample.tech_dict[techs[index]].repair_cost / discount_array
        if techs[index] != 'Null':
            ind = index - null_correct
            repair_cost_null[ind, :] = sample.tech_dict[techs[index]].repair_cost / discount_array - null_repair_cost
            gas_value_n[ind] = sum((null_emissions-sample.tech_dict[techs[index]].emissions) / discount_array) * \
                               sample.time.delta_t*24*3600*sample.econ_settings.gas_price
            capital_cost_null[ind, :] = sample.tech_dict[techs[index]].capital / discount_array
            find_cost_null[ind, :] = sample.tech_dict[techs[index]].find_cost / discount_array
            maintenance_cost_null[ind, :] = sample.tech_dict[techs[index]].maintenance / discount_array
        else:
            null_correct += 1
    # consolidate costs into totals
    cap, find, main = np.sum(capital_cost, 1), np.sum(find_cost, 1), np.sum(maintenance_cost, 1)
    no_repair = np.sum(repair_cost_no_repair, 1)
    no_repair_npv = {'Capital': cap, 'Maintenance': main, 'Repair': no_repair, 'Finding': find}
    cap_n, find_n = np.sum(capital_cost_null, 1), np.sum(find_cost_null, 1)
    main_n = np.sum(maintenance_cost_null, 1)
    n_repair = np.sum(repair_cost_null, 1)
    null_npv = {'Capital': cap_n, 'Maintenance': main_n, 'Repair': n_repair, 'Finding': find_n, 'Gas': gas_value_n}
    tot_n = -cap_n - main_n - n_repair - find_n + gas_value_n
    null_npv['Total'] = tot_n
    no_repair_npv['Gas'] = gas_value_no_repair
    no_repair_npv['Total'] = -cap - main - no_repair - find + gas_value_no_repair
    # adjust units
    # for key in no_repair_npv.keys():
    #     no_repair_npv[key] /= site_count*1000
    # for key in null_npv.keys():
    #     null_npv[key] /= site_count*1000
    return no_repair_npv, null_npv
