import matplotlib.pyplot as plt
from pickle import load
from os import listdir, getcwd
from os.path import isfile, join
import numpy as np
from feast.ResultsProcessing import results_analysis_functions
from matplotlib import rc, rcParams
import json

# The color set defined here is the Tableau 10 color set.
color_set = np.array([
    [78, 121, 167],
    [242, 142,  43],
    [225,  87,  89],
    [118, 183, 178],
    [89, 161,  79],
    [237, 201, 72],
    [176, 122, 161],
    [255, 157, 167],
    [156, 117,  95],
    [186, 176, 172]
]) / 255


def plot_fixer(fig=None, ax=None, fsize=18, color=(0, 0, 0), tight_layout=True, line_width=4, fontweight='bold'):
    rc('font', weight=fontweight)
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    # fig.canvas.draw()
    if fig is None:
        fig = plt.gcf()
    if ax is None:
        ax_list = fig.get_axes()
    elif type(ax) is not list:
        ax_list = [ax]
    else:
        ax_list = ax
    for ax in ax_list:
        ax.tick_params(axis='both', which='major', labelsize=fsize, color=color)
        lw = line_width
        for ln in ax.lines:
            ln.set_linewidth(lw)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fsize)
            item.set_color(color)
        try:
            ax.legend_.prop.set_size(fsize)
            for legobj in ax.legend_.legendHandles:
                legobj.set_linewidth(4.0)
        except AttributeError:
            pass
        [i.set_linewidth(2) for i in ax.spines.values()]
        [i.set_color(color) for i in ax.spines.values()]
        ax.xaxis.label.set_fontweight(fontweight)
        ax.yaxis.label.set_fontweight(fontweight)
    if tight_layout:
        plt.tight_layout()


def time_series(results_file, line_width=2):
    """
    Display a time series of emissions from each detection method in a results file

    :param results_file:    path to a results file
    :param line_width:      width at which to plot lines
    """
    results = load(open(results_file, 'rb'))
    tech_dict = results.ldar_program_dict
    fig = plt.figure()
    ax = fig.add_subplot(111)
    counter = -1
    for tech in tech_dict.keys():
        lab = tech
        if tech.lower() == 'null':
            lab = 'No LDAR'
        counter += 1
        tech_dict[tech].line = ax.plot(np.array(range(0, results.time.n_timesteps)) * results.time.delta_t / 365,
                                       np.array(tech_dict[tech].emissions_timeseries)/tech_dict[
                                           'Null'].emissions_timeseries[0], label=lab,
                                       color=color_set[counter], linewidth=line_width)
        avg_emissions = np.mean(np.array(tech_dict[tech].emissions_timeseries)/tech_dict[tech].emissions_timeseries[0])
        ax.plot([0, results.time.end_time / 365],
                [avg_emissions, avg_emissions],
                '--', label=lab + ' Average', color=color_set[counter], linewidth=line_width)
    plt.xlabel('Time [years]', fontweight='bold')
    plt.ylabel('Fraction of initial emissions', fontweight='bold')
    plt.legend(bbox_to_anchor=(1, 0.72))


def time_series_allmc(results_dir, mc_line_width=0.5, mean_line_width=1):
    """
    Display a time series of emissions from each detection method for all Monte Carlo simulations

    :param results_dir:    directory containing the results files for a simulation
    :param mc_line_width:      width at which to plot lines for all the Monte Carlo simulations
    :param mean_line_width:    line width for the mean of all MCs Line-plot
    """
    files = [f for f in listdir(results_dir) if isfile(join(results_dir, f))
             if '.json' not in f if 'Frame.p' not in f if 'plot.pdf' not in f]

    res0 = load(open(f"{results_dir}/{files[0]}", 'rb'))
    techs = res0.ldar_program_dict.keys()
    ti = np.array([range(0, res0.time.n_timesteps)]) * res0.time.delta_t / 365
    tseries = np.repeat(ti, len(files) + 1, axis=0)
    em_tseries_dict = {}
    em_mean_dict = {}

    for tech in techs:
        em_tseries = np.array([])
        if tech not in em_mean_dict.keys():
            em_mean_dict[tech] = np.array([])
        for f in files:
            res = load(open(f"{results_dir}/{f}", 'rb'))
            if tech not in em_tseries_dict.keys():
                em_tseries_dict[tech] = None
            if len(em_tseries) == 0:
                # Normalize program emissions to the Null Scenario
                em_tseries = np.divide(
                    res.ldar_program_dict[tech].emissions_timeseries,
                    res.ldar_program_dict['Null'].emissions_timeseries[0]
                )
            if len(em_tseries) > 0:
                # Normalize program emissions to the Null Scenario
                em_tseries = np.vstack((em_tseries,
                                        np.divide(
                                            res.ldar_program_dict[tech].emissions_timeseries,
                                            res.ldar_program_dict['Null'].emissions_timeseries[0]
                                        )
                                        ))
        em_tseries_dict[tech] = em_tseries
        for em in em_tseries_dict[tech].T:
            em_mean_dict[tech] = np.append(em_mean_dict[tech], np.mean(em))
        # Plot Normalized LDAR Program Emissions
        if tech == 'Null':
            title = 'No_LDAR'
        else:
            title = tech
        fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [3, 1]})
        fig.suptitle(f"{title} timeseries", fontsize=16)
        plt.ion()
        # Line plot of all MCs and mean for given LDAR program (tech)
        ax1.plot(tseries, em_tseries_dict[tech], color='gray', linewidth=mc_line_width, label=f"{title}")
        ax1.plot(tseries[0], em_mean_dict[tech], color='red', linewidth=mean_line_width, label=f"{title} Mean")
        ax1.set_xlabel('Time (years)')
        ax1.set_ylabel('Fraction of initial emissions')

        # Histogram plot
        weights = np.ones(len(em_tseries_dict[tech].flatten()))/len(files)  # normalize counts to number of MCs
        ax2.hist(
            em_tseries_dict[tech].flatten(),
            bins=50, color='gray', orientation='horizontal', weights=weights)
        ax2.hist(
            em_mean_dict[tech],
            bins=50, color='red', orientation='horizontal')
        plt.ioff()
        plt.savefig(f"{results_dir}/{title}_timeseries_plot.pdf", dpi=600, format='pdf')
        plt.close(fig)

        # Boxplots
        box_data = {prog: techvals.flatten() for (prog, techvals) in em_tseries_dict.items()}
        fig1, ax3 = plt.subplots(figsize=(20, 10))
        plt.ion()
        ax3.boxplot(box_data.values())
        ax3.set_xticklabels([i for i in box_data.keys()])
        ax3.set_xlabel('LDAR Programs', fontsize=20)
        ax3.set_ylabel('Fraction of Initial Emissions', fontsize=20)
        ax3.tick_params(axis='both', which='major', labelsize=18)
        plt.ioff()
        plt.savefig(f"{results_dir}/LDAR_Programs_Boxplot.pdf", dpi=600, format='pdf')
        plt.close(fig1)

def time_series_allmc_json(results_dir, mc_line_width=0.5, mean_line_width=1):
    """
    Display a time series of emissions from each detection method for all Monte Carlo simulations saved as JSON files

    :param results_dir:    directory containing the results files for a simulation
    :param mc_line_width:      width at which to plot lines for all the Monte Carlo simulations
    :param mean_line_width:    line width for the mean of all MCs Line-plot
    """
    not_prog_keys = ['time']
    files = [f for f in listdir(results_dir) if isfile(join(results_dir, f))
             if '.json' in f]
    with open(f"{results_dir}/{files[0]}", 'r') as jfile:
        data = jfile.read()
    res0 = json.loads(data)
    techs = [i for i in res0.keys() if i not in not_prog_keys]
    ti = np.array([range(0, res0['time']['n_timesteps'])]) * res0['time']['n_timesteps'] / 365
    tseries = np.repeat(ti, len(files) + 1, axis=0)
    em_tseries_dict = {}
    em_mean_dict = {}

    for tech in techs:
        em_tseries = np.array([])
        if tech not in em_mean_dict.keys():
            em_mean_dict[tech] = np.array([])
        for f in files:
            with open(f"{results_dir}/{f}", 'r') as jfile:
                data = jfile.read()
            res = json.loads(data)
            if tech not in em_tseries_dict.keys():
                em_tseries_dict[tech] = None
            if len(em_tseries) == 0:
                # Normalize program emissions to the Null Scenario
                em_tseries = np.divide(
                    res[tech]['emission timeseries'],
                    res['Null']['emission timeseries'][0]
                )
            if len(em_tseries) > 0:
                # Normalize program emissions to the Null Scenario
                em_tseries = np.vstack((em_tseries,
                                        np.divide(
                                            res[tech]['emission timeseries'],
                                            res['Null']['emission timeseries'][0]
                                        )
                                        ))
        em_tseries_dict[tech] = em_tseries
        for em in em_tseries_dict[tech].T:
            em_mean_dict[tech] = np.append(em_mean_dict[tech], np.mean(em))
        # Plot Normalized LDAR Program Emissions
        if tech == 'Null':
            title = 'No_LDAR'
        else:
            title = tech
        fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [3, 1]})
        fig.suptitle(f"{title} timeseries", fontsize=16)
        plt.ion()
        # Line plot of all MCs and mean for given LDAR program (tech)
        ax1.plot(tseries, em_tseries_dict[tech], color='gray', linewidth=mc_line_width, label=f"{title}")
        ax1.plot(tseries[0], em_mean_dict[tech], color='red', linewidth=mean_line_width, label=f"{title} Mean")
        ax1.set_xlabel('Time (years)')
        ax1.set_ylabel('Fraction of initial emissions')

        # Histogram plot
        weights = np.ones(len(em_tseries_dict[tech].flatten()))/len(files)  # normalize counts to number of MCs
        ax2.hist(
            em_tseries_dict[tech].flatten(),
            bins=50, color='gray', orientation='horizontal', weights=weights)
        ax2.hist(
            em_mean_dict[tech],
            bins=50, color='red', orientation='horizontal')
        plt.ioff()
        plt.savefig(f"{results_dir}/{title}_timeseries_plot.pdf", dpi=600, format='pdf')
        plt.close(fig)

    # Boxplots
    box_data = {prog: techvals.flatten() for (prog, techvals) in em_tseries_dict.items()}
    fig1, ax3 = plt.subplots(figsize=(20,10))
    plt.ion()
    ax3.boxplot(box_data.values())
    ax3.set_xticklabels(box_data.keys())
    ax3.set_xlabel('LDAR Programs', fontsize=20)
    ax3.set_ylabel('Fraction of Initial Emissions', fontsize=20)
    ax3.tick_params(axis='both', which='major', labelsize=18)
    plt.ioff()
    plt.savefig(f"{results_dir}/LDAR_Programs_Boxplot.pdf", dpi=600, format='pdf')
    plt.close(fig1)



def abatement_cost_plotter(directory, gwp=34, discount_rate=0, gas_price=0):
    """
    Generates a box plot of the cost of abatement from pickle files
    gwp defaults to 34, which is the value provided in the IPCC 5th assessment report including climate-carbon feedbacks
    (see Table 8.7, page 714 in Chapter 8 of Climate Change 2013: The Physical Science Basis.)

    :param directory: A directory containing one or more realizations of a scenario
    :param gwp: global warming potential of methane
    :param discount_rate: discount rate to use in NPV calculations
    :param gas_price: value of saved gas ($/g)
    :return:
    """
    npv, emissions, techs = results_analysis_functions.results_analysis(directory, discount_rate, gas_price)
    files = [f for f in listdir(directory) if isfile(join(directory, f)) if '.json' not in f if 'Frame.p' not in f]
    with open(directory + '/' + files[0], 'rb') as f:
        sample = load(f)
    emissions = np.sum(emissions * sample.time.delta_t * 3600 * 24 / 1e6, axis=1)  # metric tonnes
    em_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    cost_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    for ind in range(em_abate.shape[0]):
        em_abate[ind, :] = emissions[-1, :] - emissions[ind, :]
        cost_abate[ind, :] = -npv['Total'][ind, :]
    abatement_cost = cost_abate / em_abate / gwp
    medianprops = dict(color='k')
    boxprops = dict(linewidth=4)
    boxplot = plt.boxplot(np.transpose(abatement_cost), medianprops=medianprops,
                          boxprops=boxprops, patch_artist=True)
    ind = 1
    for bx in boxplot['boxes']:
        bx.set(facecolor=color_set[ind])
        ind += 1
    ax = plt.gca()
    ax.set_xticklabels(techs[:len(techs) - 1])
    ax.set_ylabel('Mitigation cost\n(USD/metric ton CO$_2$ eq.)')
    ax.set_xlabel('LDAR program')
    plot_fixer()


def abatement_cost_plotter_json(directory, gwp=34, discount_rate=0, gas_price=0):
    """
    Generates a box plot of the cost of abatement from json files
    gwp defaults to 34, which is the value provided in the IPCC 5th assessment report including climate-carbon feedbacks
    (see Table 8.7, page 714 in Chapter 8 of Climate Change 2013: The Physical Science Basis.)

    :param directory: A directory containing one or more realizations of a scenario
    :param gwp: global warming potential of methane
    :param discount_rate: discount rate to use in NPV calculations
    :param gas_price: value of saved gas ($/g)
    :return:
    """
    npv, emissions, techs = results_analysis_functions.results_analysis(directory, discount_rate, gas_price)
    files = [f for f in listdir(directory) if isfile(join(directory, f)) if '.json' in f]
    with open(f"{directory}/{files[0]}", 'r') as jfile:
        data = jfile.read()
    sample = json.loads(data)
    emissions = np.sum(emissions * sample['time']['delta_t'] * 3600 * 24 / 1e6, axis=1)  # metric tonnes
    em_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    cost_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    for ind in range(em_abate.shape[0]):
        em_abate[ind, :] = emissions[-1, :] - emissions[ind, :]
        cost_abate[ind, :] = -npv['Total'][ind, :]
    abatement_cost = cost_abate / em_abate / gwp
    medianprops = dict(color='k')
    boxprops = dict(linewidth=4)
    boxplot = plt.boxplot(np.transpose(abatement_cost), medianprops=medianprops,
                          boxprops=boxprops, patch_artist=True)
    ind = 1
    for bx in boxplot['boxes']:
        bx.set(facecolor=color_set[ind])
        ind += 1
    ax = plt.gca()
    ax.set_xticklabels(techs[:len(techs) - 1])
    ax.set_ylabel('Mitigation cost\n(USD/metric ton CO$_2$ eq.)')
    ax.set_xlabel('LDAR program')
    plot_fixer()
