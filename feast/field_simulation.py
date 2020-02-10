from .GeneralClassesFunctions import simulation_classes
from . import DetectionModules as Dm
from .GeneralClassesFunctions.simulation_functions import save_results
import numpy as np


def check_timestep(gas_field, time):
    """
    Prints a warning if time.delta_t is greater than the duration of some permitted emissions
    :param gas_field: a gas_field object
    :param time: a time object
    :return:
    """
    for sitedict in gas_field.sites.values():
        site = sitedict['parameters']
        for comp_temp in site.comp_dict.values():
            comp = comp_temp['parameters']
            if 0 < comp.episodic_emission_duration < time.delta_t:
                print("Warning: Episodic emission duration in site '{}', component '{}' is less than the simulation "
                      "time "
                      "step.\n"
                      "Episodic emissions will be treated as though they have a duration of one time step.".format(
                          site.name, comp.name))


def field_simulation(gas_field=None, input_leaks=None, dir_out='Results', time=None,
                     econ_set=None, tech_dict=None, detection_techs=None, display_status=True,
                     save_norepair_scenario=False):
    """
    field_simulation generates a single realization of scenario. The scenario is defined by the input values.
    gas_field           a GasField object
    input_leaks         a list of leaks to be generated at each timestep
    dir_out             directory name in which to save results
    time                a Time object
    econ_set            a FinanceSettings object
    tech_dict           a dict of detection technology objects
    detection_techs     a list of detection technology identifying strings
    save_norepair_scenario    if True, save the total emissions time series in a no-repair scenario
    """
    # -------------- Define settings --------------
    # time defines parameters related to time in the model. Time units are days.
    if time is None:
        time = simulation_classes.Time()

    if gas_field is None:
        gas_field = simulation_classes.GasField(time)

    # Note: econ_settings are not used during the simulation, but are saved for use in post simulation data processing
    if econ_set is None:
        econ_set = simulation_classes.FinanceSettings()

    if detection_techs is None:
        detection_techs = ['null.Null', 'tech_detect.TechDetect']

    new_leaks, no_repair_emissions = [], []
    if tech_dict is None:
        tech_dict = dict()
        for tech in detection_techs:
            [tech_mod, tech_class] = tech.split(sep='.')
            tech_dict[tech_class] = eval('Dm.' + tech_mod.lower() + '.' + tech_class + '(time, gas_field)')
    elif 'Null' not in tech_dict:
        tech_dict['Null'] = Dm.null.Null(time, gas_field)
    # -------------- Run the simulation --------------
    # check_timestep(gas_field, time)
    for time.time_index in range(0, time.n_timesteps):
        if display_status and time.current_time % int(time.end_time/10) < time.delta_t:
            print("Currently evaluating time step " + str(time.time_index) + " of " + str(time.n_timesteps))
        if gas_field.input_leaks is None:
            new_leaks.append(gas_field.leak_size_maker(time))
        else:
            new_leaks.append(gas_field.input_leaks[time.time_index])
        if save_norepair_scenario:
            no_repair_emissions.append(no_repair_emissions[-1] + np.sum(new_leaks.flux))
        # Loop through each LDAR program:
        for tech_obj in tech_dict.values():
            # The extra factor here accounts for emissions that end part way through a timestep
            if len(tech_obj.leaks.flux) > 0:
                timestep_fraction = (tech_obj.leaks.endtime - time.current_time + time.delta_t) / time.delta_t
                emissions = tech_obj.leaks.flux * np.min([np.ones(len(tech_obj.leaks.endtime)),
                                                               timestep_fraction], axis=0)
                tech_obj.emissions.append(np.sum(emissions))
                tech_obj.vents.append(np.sum(emissions[np.invert(tech_obj.leaks.reparable)]))
            else:
                tech_obj.emissions.append(0)
                tech_obj.vents.append(0)
            tech_obj.detection(time, gas_field)
            tech_obj.leaks.extend(new_leaks[-1])
            if time.time_index % 1000 == 0:
                tech_obj.leaks.clear_zeros()
        time.current_time += time.delta_t

    # -------------- Save results --------------
    results = simulation_classes.Results(time, gas_field, tech_dict, no_repair_emissions, econ_set,
                                         new_leaks)
    save_results(dir_out, results)
