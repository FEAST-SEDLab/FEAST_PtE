import feast.EmissionSimModules.infrastructure_classes
from .EmissionSimModules import simulation_classes
from . import DetectionModules as Dm
from .EmissionSimModules.simulation_functions import save_results
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
                      "time step.\n"
                      "Episodic emissions will be treated as though they have a duration of one time step.".format(
                          site.name, comp.name))


def field_simulation(gas_field, dir_out='Results', time=None, ldar_program_dict=None,
                     econ_set=None, display_status=True):
    """
    field_simulation generates a single realization of scenario. The scenario is defined by the input values.
    gas_field           a GasField object
    dir_out             directory name in which to save results
    time                a Time object
    econ_set            a FinanceSettings object
    tech_dict           a dict of detection technology objects
    detection_techs     a list of detection technology identifying strings
    """
    # -------------- Define settings --------------
    # time defines parameters related to time in the model. Time units are days.
    if time is None:
        time = simulation_classes.Time()

    # Note: econ_settings are not used during the simulation, but are saved for use in post simulation data processing
    if econ_set is None:
        econ_set = simulation_classes.FinanceSettings()

    if ldar_program_dict is None:
        ldar_program_dict = {'Null': Dm.ldar_program.LDARProgram(time, gas_field)}

    elif 'Null' not in ldar_program_dict:
        ldar_program_dict['Null'] = Dm.ldar_program.LDARProgram(time, gas_field)
    # -------------- Run the simulation --------------
    # check_timestep(gas_field, time)
    for time.time_index in range(0, time.n_timesteps):
        if display_status and time.current_time % int(time.end_time/10) < time.delta_t:
            print("The evaluation is {:0.0f}% complete" .format(100 * time.time_index / time.n_timesteps))
        # Loop through each LDAR program:
        for ldar_program in ldar_program_dict.values():
            # The extra factor here accounts for emissions that end part way through a timestep
            if len(ldar_program.emissions.flux) > 0:
                timestep_fraction = (ldar_program.emissions.endtime[:ldar_program.emissions.n_leaks] -
                                     time.current_time + time.delta_t) / time.delta_t
                emissions = ldar_program.emissions.flux[:ldar_program.emissions.n_leaks] * \
                    np.min([np.ones(ldar_program.emissions.n_leaks), timestep_fraction], axis=0)
                ldar_program.emissions_timeseries.append(np.sum(emissions))
                ldar_program.vents_timeseries.append(np.sum(emissions[np.invert(ldar_program.emissions.reparable
                                                                            [:ldar_program.emissions.n_leaks])]))
            else:
                ldar_program.emissions_timeseries.append(0)
                ldar_program.vents_timeseries.append(0)
            ldar_program.action(time, gas_field)
            ldar_program.emissions.extend(gas_field.input_emissions[time.time_index])
            if time.time_index % 1000 == 0:
                ldar_program.emissions.clear_zeros()
        time.current_time += time.delta_t

    # -------------- Save results --------------
    results = simulation_classes.Results(time, gas_field, ldar_program_dict, econ_set)
    save_results(dir_out, results)
