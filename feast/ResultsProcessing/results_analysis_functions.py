import numpy as np
from pickle import load
from os import listdir
from os.path import isfile, join


def results_analysis(directory, discount_rate, gas_price):
    """
    Process many realizations of a single scenario stored in a directory

    :param directory:       A directory of results files all generated under the same scenario
    :return:
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
        null_npv_temp = npv_calculator(directory + '/' + files[jindex], discount_rate, gas_price)
        for key in null_npv.keys():
            # no_repair_npv[key][:, jindex] = no_repair_npv_temp[key]
            null_npv[key][:, jindex] = null_npv_temp[key]
    return null_npv, emissions_timeseries, costs, progs


def npv_calculator(filepath, discount_rate, gas_price):
    """
        Calculates the net present value (NPV) of each LDAR program in the results file

        :param filepath:        path to a results file
        :param discount_rate: The discount rate of future cash flows (should be between 0 and 1)
        :param gas_price: The value to assign to mitigated gas losses ($/gram)
        :return:
            null_npv          NPV of each LDAR program compared to a scenario with only the Null LDAR program [k$/well]
    """
    sample = load(open(filepath, 'rb'))
    # The 'Null' LDAR program is special here because null_npv is calculated with respect to it
    if 'Null' not in sample.ldar_program_dict:
        raise NameError('tech_dict must contain a "Null" detection method to use npv_calculator')
    null_emissions = np.array(sample.ldar_program_dict['Null'].emissions_timeseries)
    time = np.linspace(0, sample.time.n_timesteps * sample.time.delta_t, sample.time.n_timesteps)
    discount_array = (1 + discount_rate)**(time / 365)
    # Initialize all arrays
    gas_value_n = np.zeros(len(sample.ldar_program_dict) - 1)
    find_cost = np.zeros([len(sample.ldar_program_dict), sample.time.n_timesteps])
    find_cost_null = np.zeros([len(sample.ldar_program_dict) - 1, sample.time.n_timesteps])
    repair_cost_null = np.zeros(len(sample.ldar_program_dict) - 1)
    # Store the cost of repairs in the null module
    em = sample.ldar_program_dict['Null'].emissions.emissions
    em_tc = em.loc[em.end_time < time[-1], ['end_time', 'repair_cost']]
    null_repair_cost = np.sum(em_tc['repair_cost'] / ((1 + discount_rate)**(em_tc['end_time'] / 365)))
    # Calculate costs associated with each LDAR program
    null_correct = 0
    index = 0
    for pr_name, program in sample.ldar_program_dict.items():
        find_cost[index, :] = program.find_cost / discount_array
        if pr_name != 'Null':
            em = program.emissions.emissions
            em_tc = em.loc[em.end_time < time[-1], ['end_time', 'repair_cost']]
            repair_cost_null[index] = np.sum(em_tc['repair_cost'] / ((1 + discount_rate)**(em_tc['end_time'] / 365)))
            repair_cost_null[index] -= null_repair_cost
            gas_value_n[index] = sum((null_emissions - program.emissions_timeseries) / discount_array) * \
                sample.time.delta_t * 24 * 3600 * gas_price
            find_cost_null[index, :] = (program.find_cost - sample.ldar_program_dict['Null'].find_cost) / discount_array
            index += 1
    # consolidate costs into totals
    find_n = np.sum(find_cost_null, 1)
    null_npv = {'Repair': repair_cost_null, 'Finding': find_n, 'Gas': gas_value_n}
    tot_n = gas_value_n - repair_cost_null - find_n
    null_npv['Total'] = tot_n
    return null_npv
