import copy
import numpy as np
from ..EmissionSimModules.simulation_functions import set_kwargs_attrs
from scipy import interpolate as interp


class DetectionMethod:
    """
    DetectionMethod is an abstract super class that defines the form required for all detection methods
    """
    def __init__(self, time, **kwargs):
        """
        Inputs:
            time         a time object (Defined in simulation_classes)
            kwargs       optional input dict that will override default parameters
        """
        self.find_cost = np.zeros(time.n_timesteps)
        self.repair_cost = np.zeros(time.n_timesteps)
        self.op_envelope = {}
        self.detection_variables = {}  # Dict with format {name: interpolation mode}
        self.hull = None
        self.site_queue = []
        # Set all attributes defined in kwargs, regardless of whether they already exist
        set_kwargs_attrs(self, kwargs, only_existing=True)
        if type(self.detection_variables) is not dict:
            raise TypeError("Detection_variables must be a dict of form {name: interpolation mode,}")

    def check_time(self, time):
        """
        Determines whether or not the detection method is active during the present time step
        :param time:
        :return:
        """
        oktime = self.ophrs['begin'] <= np.mod(time.current_time, 1) * 24 < self.ophrs['end']
        # accounts for a delta_t that is greater than the daily working hours
        timesize = time.delta_t * 24 > self.ophrs['end'] - self.ophrs['begin']
        return oktime or timesize

    def check_op_envelope(self, gas_field, time, site_index=None):
        """
        Returns false if an operating envelope condition fails
        The method supports 8 classes of operating envelope conditions:
        1.) A meteorological condition based on min-max values that apply to the whole field (eg. temperature)
        2.) A meteorological condition based on min-max values that are site-specific (eg. wind direction)
        3.) A meteorological condition based on a fail list that applies to the whole field (eg. precipitation type)
        4.) A meteorological condition based on a fail list that is site specific (possible but not expected)
        5.) A site condition based on min-max values that apply to the whole field (eg site production)
        6.) A site condition based on min-max values that are site-specific (possible but not expected)
        7.) A site condition based on a fail list that applies to the whole field (eg. site type)
        8.) A site condition based on a fail list that is site specific (possible but not expected).
        :param gas_field:
        :param time:
        :param site_index:
        :return:
        """
        status = 'field pass'
        # iterate across all operating envelope conditions
        for name, params in self.op_envelope.items():
            # Load a meteorological condition or a site attribute, depending on the class of envelope parameter
            if params['class'] in [1, 2, 3, 4]:
                if 'interp_mode' in params:
                    im = params['interp_mode']
                else:
                    im = 'mean'
                condition = gas_field.get_met(time, name, interp_modes=im, ophrs=self.ophrs)[name]
            else:
                site = gas_field.sites[self.find_site_name(gas_field, site_index)]
                condition = site.op_env_params[name]
            if params['class'] == 1 and not self.check_min_max_condition(condition, params):
                return 'field fail'
            elif params['class'] in [2, 6]:
                status = 'site pass'
                site_params = {'min': params['min'][site_index], 'max': params['max'][site_index]}
                if not self.check_min_max_condition(condition, site_params):
                    status = 'site fail'
            elif params['class'] == 3 and condition in params['enum_fail_list']:
                return 'field fail'
            elif params['class'] in [4, 8]:
                status = 'site pass'
                if condition in params['enum_fail_list'][site_index]:
                    status = 'site fail'
            elif params['class'] == 5:
                status = 'site pass'
                if self.check_min_max_condition(condition, params):
                    status = 'site fail'
            elif params['class'] == 7:
                status = 'site pass'
                if condition in params['enum_fail_list']:
                    status = 'site fail'
        return status

    def choose_sites(self, gas_field, time, n_sites, clear_sites=True):
        """
        Identifies sites to survey at this time step
        :param gas_field:
        :param time:
        :param n_sites: Max number of sites to survey at this time step
        :param clear_sites: If true, clear sites selected from the queue. If False, leave sites in the queue
                            (eg. continuous monitors)
        :return:
        """
        site_inds = []
        queue_ind = 0
        while len(site_inds) < n_sites and queue_ind < len(self.site_queue):
            op_env = self.check_op_envelope(gas_field, time, self.site_queue[queue_ind])
            if op_env == 'field pass':
                # This case applies if there are no site-specific operating envelope conditions.
                site_inds = self.site_queue[:n_sites]
                if clear_sites:
                    del self.site_queue[:n_sites]
            elif op_env == 'site pass':
                # This case applies if the operating envelope is satisifed for this site,
                # but may fail for a different site.
                if clear_sites:
                    site_inds.append(self.site_queue.pop(queue_ind))
                else:
                    site_inds.append(self.site_queue[queue_ind])
                    queue_ind += 1
            elif op_env == 'field fail':
                # This case applies if the operating envelope fails for the entire field at this time step.
                return site_inds
            else:
                # This condition applies if the site fails the operating envelope but other sites may pass.
                queue_ind += 1
        return site_inds

    @staticmethod
    def check_min_max_condition(condition, params):
        """
        Checks a min-max condition defined by params. Supports float, integer, list and array based min max conditions.
        If the min-max condition is specified as a min float/integer and max float/integer, the numbers are placed in
        min and max lists each with length 1. The function returns True if the condition is between the min and max
        values, False otherwise. If the min and max values are array-like, the function returns true if the condition is
        between any pair of min-max values.
        :param condition: condition to check (must be a number)
        :param params: a dict with 'min' and 'max' keys. The min and max values can be numbers or array-like.
        :return:
        """
        try:
            _ = len(params['min'])
        except TypeError:
            params['min'] = [params['min']]
            params['max'] = [params['max']]
        condition_allowed = False
        for ind in range(len(params['min'])):
            if params['min'][ind] <= condition <= params['max'][ind]:
                condition_allowed = True
        return condition_allowed

    @staticmethod
    def find_site_name(gas_field, site_index):
        """
        Determines the key for a site based on  its index
        :param gas_field:
        :param site_index:
        :return:
        """
        for sitename, site in gas_field.sites.items():
            if site['parameters'].site_inds[0] <= site_index < site['parameters'].site_inds[1]:
                return sitename
        return -1

    @staticmethod
    def find_comp_name(gas_field, sitename, comp_index):
        """
        Determines the key for a component based on its index and site
        :param gas_field:
        :param sitename:
        :param comp_index
        :return:
        """
        for compname, comp in gas_field.sites[sitename]['parameters'].comp_dict.items():
            if comp.comp_inds[0] <= comp_index < comp.comp_inds[1]:
                return compname
        return -1

    def get_current_conditions(self, time, gas_field, emissions, em_indexes):
        """
        Extracts conditions specified in self.detection_variables
        :param time: simulation time object
        :param gas_field: simulation gas field object
        :param emissions: emissions object
        :param em_indexes: indexes of emissions to consider
        :return conditions: an array (n_emissions, n_variables) of conditions for use in the PoD calculation
        """
        conditions = np.zeros([len(em_indexes), len(self.detection_variables)])
        index = 0
        for v, im in self.detection_variables.items():
            if v in gas_field.met:
                conditions[:, index] = gas_field.get_met(time, v, interp_modes=im, ophrs=self.ophrs)[v]
            else:
                conditions[:, index] = emissions.__getattribute__(v)[em_indexes]
            index += 1
        return conditions

    @staticmethod
    def empirical_interpolator(test_conditions, test_results, sim_conditions):
        """
        calculates the probabiity of detection analytically
        :param test_conditions: conditions to be interpolated from
        :param test_results: results associated with each condition listed in test_conditions
        :param sim_conditions: Nxk array of current conditions, where N is the number of emissions to consider,
            and k is the number of conditions
        :return:
        """
        probs = interp.griddata(test_conditions, test_results, sim_conditions)
        # griddata returns NaN for vars outside the convex hull of the interpolation data points when using default
        # linear interpolation. The following code sets those NaN values (outside the convex hull)
        # to the nearest interpolation point.
        cond = np.where(np.isnan(probs))[0]
        probs[cond] = interp.griddata(test_conditions, test_results, sim_conditions[cond], method='nearest')
        return np.ndarray.flatten(probs)

    def extend_site_queue(self, site_inds):
        """
        Add new sites to the site_queue if they are not already in the queue
        :param site_inds: List of indexes to add to the queue
        :return:
        """
        for si in site_inds:
            if si not in self.site_queue:
                self.site_queue.append(si)
