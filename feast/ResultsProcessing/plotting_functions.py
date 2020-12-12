import matplotlib.pyplot as plt
from pickle import load
from os import listdir
from os.path import isfile, join
import numpy as np
from feast.ResultsProcessing import results_analysis_functions
from matplotlib import rc, rcParams

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


def time_series(results_file, line_width=6):
    """
    Display a time series of emissions from each detection method in a results file

    :param results_file:    path to a results file
    :param line_width:      width at which to plot lines
    """
    results = load(open(results_file, 'rb'))
    tech_dict = results.tech_dict
    fig = plt.figure()
    ax = fig.add_subplot(111)
    counter = -1
    for tech in tech_dict.keys():
        lab = tech
        if tech.lower() == 'null':
            lab = 'No LDAR'
        counter += 1
        tech_dict[tech].line = ax.plot(np.array(range(0, results.time.n_timesteps)) * results.time.delta_t / 365,
                                       np.array(tech_dict[tech].emissions)/tech_dict[tech].emissions[0], label=lab,
                                       color=color_set[counter], linewidth=line_width)
        avg_emissions = np.mean(np.array(tech_dict[tech].emissions)/tech_dict[tech].emissions[0])
        ax.plot([0, results.time.end_time / 365],
                [avg_emissions, avg_emissions],
                '--', label=lab + ' Average', color=color_set[counter], linewidth=line_width)
    plt.xlabel('Time [years]', fontweight='bold')
    plt.ylabel('Fraction of initial emissions', fontweight='bold')
    plt.legend(bbox_to_anchor=(1, 0.72))


def abatement_cost_plotter(directory, gwp=34):
    """
    Generates a box plot of the cost of abatement
    gwp defaults to 34, which is the value provided in the IPCC 5th assessment report including climate-carbon feedbacks
    (see Table 8.7, page 714 in Chapter 8 of Climate Change 2013: The Physical Science Basis.)

    :param directory: A directory containing one or more realizations of a scenario
    :param gwp: global warming potential of methane
    :return:
    """
    _, emissions, costs, techs = results_analysis_functions.results_analysis(directory)
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    with open(directory + '/' + files[0], 'rb') as f:
        sample = load(f)
    emissions = np.sum(emissions * sample.time.delta_t * 3600 * 24 / 1e6, axis=1)  # metric tonnes
    em_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    cost_abate = np.zeros([emissions.shape[0] - 1, emissions.shape[1]])
    for ind in range(em_abate.shape[0]):
        em_abate[ind, :] = emissions[-1, :] - emissions[ind, :]
        cost_abate[ind, :] = costs[ind, :] - costs[-1, :]
    abatement_cost = cost_abate / em_abate / gwp
    medianprops = dict(color='k')
    boxplot = plt.boxplot(np.transpose(abatement_cost), medianprops=medianprops, patch_artist=True)
    ind = 1
    for bx in boxplot['boxes']:
        bx.set(facecolor=color_set[ind])
        ind += 1
    ax = plt.gca()
    ax.set_xticklabels(techs[:2])
    ax.set_ylabel('Mitigation cost\n(USD/metric ton CO$_2$ eq.)')
    ax.set_xlabel('LDAR program')
    plot_fixer()


def summary_plotter(directory, n_wells=None, ylabel=None):
    """
    The NPV for each realization stored in 'directory is calculated and displayed in a stacked bar chart. Each component
    of the NPV is displayed separately in the chart.

    :param directory:    path to a directory containing results files
    :param n_wells:      if set to a number, then the NPV will be reported on a per well basis
    :param ylabel:       yaxis label
    """
    # load data
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    sample = load(open(directory + '/' + files[0], 'rb'))
    npv_in, emissions, costs, techs = results_analysis_functions.results_analysis(directory)
    # Normalize to the number of wells
    if n_wells is not None:
        for cost_type in npv_in:
            npv_in[cost_type] = npv_in[cost_type] / (1000 * n_wells)
        ylabel = "NPV (k$/well)"
    n_realizations = len(npv_in['Total'][0, :])
    tech_keys = []
    # Get a list of indexes that do not have data related to the null detection method
    for key in sample.tech_dict.keys():
        if key != 'Null':
            tech_keys.append(key)
    n_tech = len(sample.tech_dict)
    npv_std = np.std(npv_in['Total'], 1)
    totals = np.average(npv_in['Total'], 1)
    npv = dict()
    for key in npv_in.keys():
        npv[key] = np.sum(npv_in[key], 1)/n_realizations

    # create figure
    plt.figure(figsize=(10, 5))
    width = 0.35
    ind = np.linspace(0, n_tech-2, n_tech-1)
    counter = 0
    neg_bottoms = [0]*(n_tech-1)
    pos_bottoms = [0]*(n_tech-1)
    key_list = list(npv.keys())
    key_list.sort()
    for key in key_list:
        if key == 'Total' or key == 'Gas':
            continue
        else:
            for k in ind:
                m = int(k)
                if npv[key][m] >= 0:
                    if m == 0:
                        plt.bar(m, -npv[key][m], width, bottom=neg_bottoms[m], color=color_set[counter], label=key)
                    else:
                        plt.bar(m, -npv[key][m], width, bottom=neg_bottoms[m], color=color_set[counter])
                    neg_bottoms[m] -= npv[key][m]
                else:
                    if m == 0:
                        plt.bar(m, -npv[key][m], width, bottom=pos_bottoms[m], color=color_set[counter], label=key)
                    else:
                        plt.bar(m, -npv[key][m], width, bottom=pos_bottoms[m], color=color_set[counter])
                    pos_bottoms[m] -= npv[key][m]
        counter += 1
    plt.bar(ind, npv['Gas'], width, bottom=pos_bottoms, color=color_set[counter], label='Gas')
    plt.errorbar(ind, totals, yerr=npv_std, linewidth=8, fmt='o', ms=1, label='Total NPV', ecolor='yellow', capsize=8)
    if ylabel is None:
        plt.ylabel("NPV ($)")
    else:
        plt.ylabel(ylabel, fontsize=18)
    plt.xticks(np.linspace(0, len(tech_keys) - 1, len(tech_keys)), tech_keys, fontsize=18)
    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1.4, 1), fontsize=18)
    plot_fixer()
