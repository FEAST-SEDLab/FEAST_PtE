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


def abatement_cost_plotter(directory, gwp=34, discount_rate=0, gas_price=0):
    """
    Generates a box plot of the cost of abatement
    gwp defaults to 34, which is the value provided in the IPCC 5th assessment report including climate-carbon feedbacks
    (see Table 8.7, page 714 in Chapter 8 of Climate Change 2013: The Physical Science Basis.)

    :param directory: A directory containing one or more realizations of a scenario
    :param gwp: global warming potential of methane
    :param discount_rate: discount rate to use in NPV calculations
    :param gas_price: value of saved gas ($/g)
    :return:
    """
    npv, emissions, techs = results_analysis_functions.results_analysis(directory, discount_rate, gas_price)
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
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
