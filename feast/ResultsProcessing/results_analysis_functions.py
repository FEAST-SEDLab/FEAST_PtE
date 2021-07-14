import numpy as np
from pickle import load
from os import listdir
from os.path import isfile, join
import json
import pandas as pd

def results_analysis(directory, discount_rate, gas_price):
    """
    Process many realizations of a single scenario stored in a directory

    :param directory:       A directory of results files all generated under the same scenario
    :param discount_rate:   Discount rate to apply in net present value (NPV) calculations
    :param gas_price:       value of gas saved ($/g)
    :return:
        null_npv          array of null-NPV of each LDAR program in each realization [k$/well]
        emissions_timeseries  Array of emissions in each LDAR program in each realization at each time step
        techs           list of detection program names
    """
    files = [f for f in listdir(directory) if isfile(f"{directory}/{f}") if '.json' not in f if 'Frame.p' not in f]
    n_realizations = len(files)
    sample = load(open(join(directory, files[0]), 'rb'))
    n_prog = len(sample.ldar_program_dict)
    # Initialize an array to store a time series of emissions for every LDAR program in every realization
    emissions_timeseries = np.zeros([n_prog, sample.time.n_timesteps, n_realizations])
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
        # iterate through each category of value
        null_npv_temp = npv_calculator(directory + '/' + files[jindex], discount_rate, gas_price)
        for key in null_npv.keys():
            # no_repair_npv[key][:, jindex] = no_repair_npv_temp[key]
            null_npv[key][:, jindex] = null_npv_temp[key]
    return null_npv, emissions_timeseries, progs


def results_analysis_json(directory, discount_rate, gas_price):
    """
    Process many JSON realizations of a single scenario stored in a directory

    :param directory:       A directory of results files all generated under the same scenario
    :param discount_rate:   Discount rate to apply in net present value (NPV) calculations
    :param gas_price:       value of gas saved ($/g)
    :return:
        null_npv          array of null-NPV of each LDAR program in each realization [k$/well]
        emissions_timeseries  Array of emissions in each LDAR program in each realization at each time step
        costs                 Array of costs associated with each LDAR program (no discounting, all costs summed)
        techs           list of detection program names
    """
    not_prog_keys = ['time']
    files = [f for f in listdir(directory) if isfile(join(directory, f)) and '.json' in f]
    n_realizations = len(files)
    with open(f"{directory}/{files[0]}",'r') as jfile:
        data = jfile.read()
    sample = json.loads(data)
    n_prog = len([i for i in sample.keys() if i not in not_prog_keys])
    # Initialize an array to store a time series of emissions for every LDAR program in every realization
    emissions_timeseries = np.zeros([n_prog, sample['time']['n_timesteps'], n_realizations])
    costs = np.zeros([n_prog, n_realizations])
    progs = [i for i in sample.keys() if i not in not_prog_keys]
    null_npv = dict()
    npv_keys = ['Repair', 'Finding', 'Gas', 'Total']
    for key in npv_keys:
        null_npv[key] = np.zeros([n_prog-1, n_realizations])
    # iterate through each realization
    for jindex in range(0, len(files)):
        with open(f"{directory}/{files[jindex]}", 'r') as jfile:
            data = jfile.read()
        sample = json.loads(data)
        # iterate through each LDAR program
        for index in range(0, len(progs)):
            emissions_timeseries[index, :, jindex] = sample[progs[index]]['emission timeseries']
        # iterate through each category of value
        null_npv_temp = npv_calculator_json(directory + '/' + files[jindex], discount_rate, gas_price)
        for key in null_npv.keys():
            # no_repair_npv[key][:, jindex] = no_repair_npv_temp[key]
            null_npv[key][:, jindex] = null_npv_temp[key]
    return null_npv, emissions_timeseries, progs


def npv_calculator_json(filepath, discount_rate, gas_price):
    """
        Calculates the net present value (NPV) of each LDAR program in the  JSON results file

        :param filepath:        path to a results file
        :param discount_rate: The discount rate of future cash flows (should be between 0 and 1)
        :param gas_price: The value to assign to mitigated gas losses ($/gram)
        :return:
            null_npv          NPV of each LDAR program compared to a scenario with only the Null LDAR program [k$/well]
    """
    not_prog_keys = ['time']
    with open(filepath, 'r') as jfile:
        data = jfile.read()
    sample = json.loads(data)
    # The 'Null' LDAR program is special here because null_npv is calculated with respect to it
    if 'Null' not in [i for i in sample.keys() if i not in not_prog_keys]:
        raise NameError('tech_dict must contain a "Null" detection method to use npv_calculator')
    null_emissions = np.array(sample['Null']['emission timeseries'])
    # Initialize all arrays
    gas_value_n = np.zeros(len([i for i in sample.keys() if i not in not_prog_keys]) - 1)
    find_cost = np.zeros(len([i for i in sample.keys() if i not in not_prog_keys]) - 1)
    repair_cost_null = np.zeros(len([i for i in sample.keys() if i not in not_prog_keys]) - 1)
    # Store the cost of repairs in the null module
    em = pd.read_pickle(f"{filepath[:-5]}_emissionsDataFrame.p")
    em_tc = em.loc[em.end_time < sample['time']['end_time'], ['end_time', 'repair_cost']]
    cond = em_tc['end_time'] < sample['time']['end_time']
    null_repair_cost = np.sum(em_tc['repair_cost'][cond] / ((1 + discount_rate)**(em_tc['end_time'][cond] / 365)))
    # Calculate costs associated with each LDAR program
    time_array = np.linspace(0, sample['time']['end_time'] - sample['time']['delta_t'], sample['time']['n_timesteps'])
    discount_array = (1 + discount_rate)**(time_array / 365)
    index = 0
    for pr_name, program in sample.items():
        if pr_name not in ['Null', 'time']:
            for prog in sample[pr_name].keys():
                if prog not in ['emission timeseries', 'vent timeseries', 'repair cost', 'ogi repair']:
                    if len(sample[pr_name][prog]['deployment costs']['time_value']) > 0:
                        t = np.array(sample[pr_name][prog]['deployment costs']['time_value'])[:, 0]
                        v = np.array(sample[pr_name][prog]['deployment costs']['time_value'])[:, 1]
                        find_cost[index] += np.sum(v / (1 + discount_rate) ** (t / 365))
            repair_cost_null[index] = np.sum(em_tc['repair_cost'] / ((1 + discount_rate)**(em_tc['end_time'] / 365)))
            repair_cost_null[index] -= null_repair_cost
            gas_value_n[index] = sum((null_emissions - sample[pr_name]['emission timeseries']) / discount_array) * \
                sample['time']['delta_t'] * 24 * 3600 * gas_price
            index += 1
    # consolidate costs into totals

    null_npv = {'Repair': repair_cost_null, 'Finding': find_cost, 'Gas': gas_value_n}
    tot_n = gas_value_n - repair_cost_null - find_cost
    null_npv['Total'] = tot_n
    return null_npv


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
    time = sample.time
    # Initialize all arrays
    gas_value_n = np.zeros(len(sample.ldar_program_dict) - 1)
    find_cost = np.zeros(len(sample.ldar_program_dict) - 1)
    repair_cost_null = np.zeros(len(sample.ldar_program_dict) - 1)
    # Store the cost of repairs in the null module
    em = sample.ldar_program_dict['Null'].emissions.emissions
    em_tc = em.loc[em.end_time < time.end_time, ['end_time', 'repair_cost']]
    cond = em_tc['end_time'] < time.end_time
    null_repair_cost = np.sum(em_tc['repair_cost'][cond] / ((1 + discount_rate)**(em_tc['end_time'][cond] / 365)))
    # Calculate costs associated with each LDAR program
    time_array = np.linspace(0, time.end_time - time.delta_t, time.n_timesteps)
    discount_array = (1 + discount_rate)**(time_array / 365)
    index = 0
    for pr_name, program in sample.ldar_program_dict.items():
        if pr_name != 'Null':
            for dm in program.tech_dict.values():
                if len(dm.deployment_cost.time_value) > 0:
                    t = np.array(dm.deployment_cost.time_value)[:, 0]
                    v = np.array(dm.deployment_cost.time_value)[:, 1]
                    find_cost[index] += np.sum(v / (1 + discount_rate)**(t / 365))
            em = program.emissions.emissions
            em_tc = em.loc[em.end_time < time.end_time, ['end_time', 'repair_cost']]
            repair_cost_null[index] = np.sum(em_tc['repair_cost'] / ((1 + discount_rate)**(em_tc['end_time'] / 365)))
            repair_cost_null[index] -= null_repair_cost
            gas_value_n[index] = sum((null_emissions - program.emissions_timeseries) / discount_array) * \
                sample.time.delta_t * 24 * 3600 * gas_price
            index += 1
    # consolidate costs into totals

    null_npv = {'Repair': repair_cost_null, 'Finding': find_cost, 'Gas': gas_value_n}
    tot_n = gas_value_n - repair_cost_null - find_cost
    null_npv['Total'] = tot_n
    return null_npv
