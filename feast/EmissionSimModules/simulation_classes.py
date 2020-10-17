"""
    simulation_classes stores the classes used to represent time, results and financial settings in simulations.
"""
import os
import pickle
from .. import DetectionModules as Dm
import numpy as np


class Time:
    """
    Instances of the time class store all time related information during a simulation
    """
    def __init__(self, delta_t=1, end_time=365, current_time=0):
        """
        :param delta_t: length of one timestep (days)
        :param end_time: length of the simulation (days)
        :param current_time: current time in a simulation (days)
        """
        self.n_timesteps = int(end_time/delta_t+1)
        self.end_time = end_time
        self.delta_t = delta_t
        self.current_time = current_time
        self.time_index = 0


class Scenario:
    """
    A class to store all data specifying a scenario and the methods to run and save a realization
    """
    def __init__(self, time, gas_field, ldar_program_dict):
        """
        :param time: Time object
        :param gas_field: GasField object
        :param ldar_program_dict: dict of detection methods and associated data
        """
        self.time = time
        self.gas_field = gas_field
        self.ldar_program_dict = ldar_program_dict

    def run(self, dir_out="Results", display_status=True):
        """
        run generates a single realization of a scenario.

        :param dir_out: path to a directory in which to save results (string)
        :param display_status: if True, display a status update whenever 10% of the time steps are completed
        :return: None
        """
        # -------------- Define settings --------------
        # time defines parameters related to time in the model. Time units are days.
        if 'Null' not in self.ldar_program_dict:
            self.ldar_program_dict['Null'] = Dm.ldar_program.LDARProgram(self.time, self.gas_field, tech_dict={})
        # -------------- Run the simulation --------------
        # check_timestep(gas_field, time)
        for self.time.time_index in range(0, self.time.n_timesteps):
            if display_status and self.time.current_time % (self.time.end_time / 10) < self.time.delta_t:
                print("The evaluation is {:0.0f}% complete".format(100 * self.time.time_index / self.time.n_timesteps))
            # Loop through each LDAR program:
            for ldar_program in self.ldar_program_dict.values():
                # The extra factor here accounts for emissions that end part way through the present timestep.
                if len(ldar_program.emissions.flux) > 0:
                    timestep_fraction = (ldar_program.emissions.endtime[:ldar_program.emissions.n_em] -
                                         self.time.current_time + self.time.delta_t) / self.time.delta_t
                    emissions = ldar_program.emissions.flux[:ldar_program.emissions.n_em] * \
                                np.min([np.ones(ldar_program.emissions.n_em), timestep_fraction], axis=0)
                    ldar_program.emissions_timeseries.append(np.sum(emissions))
                    vented_em_indexes = np.invert(ldar_program.emissions.reparable[:ldar_program.emissions.n_em])
                    ldar_program.vents_timeseries.append(np.sum(emissions[vented_em_indexes]))
                else:
                    ldar_program.emissions_timeseries.append(0)
                    ldar_program.vents_timeseries.append(0)
                ldar_program.action(self.time, self.gas_field)
                ldar_program.emissions.extend(self.gas_field.input_emissions[self.time.time_index])
                if self.time.time_index % 1000 == 0:
                    ldar_program.emissions.clear_zeros()
            self.time.current_time += self.time.delta_t

        # -------------- Save results --------------
        self.save(dir_out)

    def save(self, dir_out):
        """
        Save results to a file

        :param dir_out: Name of directory in which to save output file.
        """

        if not os.path.exists(dir_out):
            os.makedirs(dir_out)
        n_realization = len(os.listdir(dir_out))
        file_out = dir_out + '/realization' + str(n_realization) + '.p'
        pickle.dump(self, open(file_out, 'wb'))

    def check_timestep(self):
        """
        Prints a warning if time.delta_t is greater than the duration of some permitted emissions

        :param gas_field: a GasField object
        :param time: a Time object
        :return: None
        """
        for sitedict in self.gas_field.sites.values():
            site = sitedict['parameters']
            for comp_temp in site.comp_dict.values():
                comp = comp_temp['parameters']
                if 0 < comp.episodic_emission_duration < self.time.delta_t:
                    print(
                        "Warning: Episodic emission duration in site '{}', component '{}' is less than the simulation "
                        "time step.\n"
                        "Episodic emissions will be treated as though they have a duration of one time step.".format(
                            site.name, comp.name))

