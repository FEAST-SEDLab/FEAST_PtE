import matplotlib.pyplot as plt
from pickle import load
from os import listdir
from os.path import isfile, join
import numpy as np
from . import results_analysis_functions
from matplotlib import rc, rcParams

# Stanford bright color palette (https://identity.stanford.edu/overview/color)
# color_set = [[140, 21, 21], [0, 152, 219], [0, 155, 118],  [178, 111, 22], [83, 40, 79],  [0, 0, 0], [77, 79, 83]]

color_set = np.array([
    [0, 61, 50],
    [193, 176, 92],
    [72, 24, 21],
    [48, 155, 255],
    [22, 104, 9],
    [0, 0, 0],
    [32, 25, 55]
]) / 255


def plot_fixer(fig=None, ax=None, fsize=18, color=[0, 0, 0], tight_layout=True):
    rc('font', weight='bold')
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
        lw = 4
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
        ax.xaxis.label.set_fontweight('bold')
        ax.yaxis.label.set_fontweight('bold')
    if tight_layout:
        plt.tight_layout()


# for i in range(0, len(color_set)):
#     for j in range(0, len(color_set[i])):
#         color_set[i][j] /= 255


def display_save(display, save, file_out):
    """
    display_save can display and/or save a plot
    Inputs
        display         choose whether to display the plot or not
        save            choose whether to save the plot or not
        file_out        define a path at which to save the file
    """
    if save:
        if file_out is None:
            file_out = 'Figures/time_series' + str(len(listdir('Figures'))) + '.pdf'
        plt.savefig(file_out, bbox_inches="tight")

    if display:
        plt.show()
    else:
        plt.close()
    plt.ioff()


def display_settings(interactive=True, inline=True):
    """
    If running in IPYTHON, display_settings uses magic commands to plot interactively
    interactive         determines whether the plot is interactive
    """
    if interactive:
        # display settings depend on the python environment being used
        try:
            __IPYTHON__
            ipy = get_ipython()
            if inline is True:
                try:
                    ipy.magic("matplotlib inline")
                except:
                    ipy.magic("matplotlib")
            else:
                ipy.magic("matplotlib")
        except:
            plt.ion()


def time_series(results_file, display=True, save=False, file_out=None, interactive=True):
    """
    Display a time series of emissions from each detection method in a results file
    Inputs:
        results_file    path to a results file
        display         choose whether to display the plot or not
        save            choose whether to save the plot or not
        file_out        define a path at which to save the file
        interactive     choose whether to make the plot interactive
    """
    display_settings(interactive)
    results = load(open(results_file, 'rb'))
    tech_dict = results.tech_dict
    fig = plt.figure()
    ax = fig.add_subplot(111)
    counter = -1
    for tech in tech_dict.keys():
        counter += 1
        tech_dict[tech].line = ax.plot(np.array(range(0, results.time.n_timesteps)) * results.time.delta_t,
                                       np.array(tech_dict[tech].emissions)/tech_dict[tech].emissions[0], label=tech,
                                       color=color_set[counter])
    plt.xlabel('Time [days]', fontweight='bold')
    plt.ylabel('Fraction of initial emissions', fontweight='bold')
    plt.legend()
    display_save(display, save, file_out)
    return ax


def summary_plotter(directory, display=True, save=False, file_out=None, interactive=True, n_wells=None, ylabel=None):
    """
    The NPV for each realization stored in 'directory is calculated and displayed in a stacked bar chart. Each component
    of the NPV is displayed separately in the chart.
    Inputs:
        directory    path to a directory containing results files
        display      boolean value that determines whether or not the plot is displayed
        save         boolean value that determines whether the plot is saved
        file_out     file_out path to a file in which to save the figure if save is True
        interactive  determines whether the plot is produced interactively
    """
    # load data
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    sample = load(open(directory + '/' + files[0], 'rb'))
    _, npv_in, _, emissions = results_analysis_functions.results_analysis(directory)
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
    display_settings(interactive)
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
    npv_std /= np.sqrt(n_realizations)
    plt.errorbar(ind, totals, yerr=npv_std, linewidth=8, fmt='o', ms=1, label='Total NPV', ecolor='yellow')
    if ylabel is None:
        plt.ylabel("NPV ($)")
    else:
        plt.ylabel(ylabel, fontsize=18)
    plt.xticks(np.linspace(0, len(tech_keys) - 1, len(tech_keys)), tech_keys, fontsize=18)
    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1.4, 1), fontsize=18)
    plot_fixer()
    display_save(display, save, file_out)


def hist_plotter(directory, display=True, save=False, file_out=None, interactive=True, inline=False, bins=None):
    """
    Plots histogram of leaks found by each technology based on results in the folder 'Directory'
    Inputs:
        directory     Path to a results folder to be analyzed
        display      boolean value that determines whether or not the plot is displayed
        save         boolean value that determines whether the plot is saved
        file_out     file_out path to a file in which to save the figure if save is True
        interactive  determines whether the plot is produced interactively
        inline       determines whether the plot is displayed inline if possible
        bins         defines the bins used in the histograms
    Return: None
    """

    display_settings(interactive, inline=inline)
    _, _, found, _ = results_analysis_functions.results_analysis(directory)
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    leaks_made = []
    sample = load(open(directory + '/' + files[0], 'rb'))
    techs = sample.tech_dict.keys()
    tech_leaks_found = {}
    well_count = 0
    # Load data from files
    for tech in techs:
        tech_leaks_found[tech] = []
    for file in files:
        with open(directory + '/' + file, 'rb') as f:
            sample = load(f)
            leaks_made.extend(sample.gas_field.initial_leaks.flux)
            well_count += sample.gas_field.site_count
            for leak_list in sample.new_leaks:
                leaks_made.extend(leak_list.flux)
            for tech in techs:
                tech_leaks_found[tech].extend(sample.tech_dict[tech].leaks_found)
    leaks_made = np.array(leaks_made)
    n_bins = 30
    if bins is None:
        bins = np.linspace(0, max(leaks_made), n_bins+1)
    cumsum_points = 100
    counts, bins, patches = plt.hist(leaks_made, bins, normed=0, alpha=0.75)
    plt.close()
    # Calculative the cumulative fraction of emissions above a range of leak fluxes
    cumsum = np.zeros(cumsum_points)
    cumsumx = np.linspace(0, max(bins), cumsum_points)
    for ind in reversed(range(0, cumsum_points-1)):
        cumsum[ind] = sum(leaks_made[leaks_made > cumsumx[ind]])
    cumsum /= cumsum[0]
    centers = 0.5*(bins[1:] + bins[:-1])
    ind = 0
    delta = bins[1] - bins[0]
    # Create a separate plot for each technology
    for tech in techs:
        if tech.lower() == 'null':
            continue
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        counts_temp, bins_temp, patches = ax1.hist(tech_leaks_found[tech], bins, normed=0, alpha=0.75)
        ax1.cla()
        counts_temp /= well_count
        ax1.bar(bins_temp[0:n_bins], counts_temp, width=delta, color=color_set[ind+1], label="Leaks found")
        ax1.plot(centers, counts/well_count, 'o', color=color_set[0], label="Leaks generated")
        ax2.plot(cumsumx, cumsum, label="Cumulative emissions fraction")
        plt.title(tech)
        ax1.set_xlabel("Leak flux (g/s)")
        ax1.set_ylabel("Leaks found per well")
        ax2.set_ylabel("Cumulative emissions fraction")
        ax1.set_ylim(0, 2)
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles = handles1
        handles.extend(handles2)
        labels = labels1
        labels.extend(labels2)
        fig.legend(handles, labels, bbox_to_anchor=(0.9, 0.7), fontsize=16)
        ind += 1
        plot_fixer()
        display_save(display, save, file_out)
