"""
result_classes defines classes that are used to store event counts and continuous variable data for saving.
"""

import numpy as np


class ResultAggregate:
    """
    A super class designed to store aggregate results during a simulation.
    Time and value pairs are stored in a list.
    """
    def __init__(self, units=None, time_value=None):
        self.time_value = time_value or []
        self.units = units

    def get_vals(self, t_start=0, t_end=np.infty):
        """
        Returns all values associated with times between t_start and t_end.

        :param t_start: Time to begin the sum
        :param t_end: time to end the sum
        :return: All values associated with time
        """
        if len(self.time_value) == 0:
            return self.time_value
        tv = np.array(self.time_value)
        vals = tv[(tv[:, 0] >= t_start) & (tv[:, 0] < t_end), 1]
        return vals

    def append_entry(self, time_value):
        """
        Add a new entry to the ResultAggregate object

        :param time_value: an ordered pair following this pattern: [time, value]
        :return: None
        """
        self.time_value.append(time_value)
        return None


class ResultDiscrete(ResultAggregate):
    """
    Designed to store discrete values associated with specific times, as opposed to continuous rates that persist
    between consecutive data points. For example, the number of sites surveyed can be recorded as a discrete data type.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_cumulative_vals(self, t_start=0, t_end=np.infty):
        """
        Returns a cumulative sum of the attribute "value"

        :param t_start: Time to begin the cumulative sum
        :param t_end: time to end the cumulative sum
        :return: Array of times in between t_start and t_end, cumulative sum of the attribute "value"
        """
        ti = np.array(self.time_value)[:, 0]
        ti = ti[(ti >= t_start) & (ti < t_end)]
        return ti, np.cumsum(self.get_vals(t_start, t_end))

    def get_sum_val(self, t_start=0, t_end=np.infty):
        """
        Returns the sum of values between t_start and t_end

        :param t_start: Time to begin the sum
        :param t_end: time to end the sum
        :return: sum of values between t_start and t_end
        """
        return np.sum(self.get_vals(t_start, t_end))


class ResultContinuous(ResultAggregate):
    """
    Designed to store continuous rates that endure between consecutive time recordings as opposed to discrete
    variables that occur at a specific time. For example, emission rate can be recorded as a continuous data type.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_time_integrated(self, start_time=0, end_time=None, unit_factor=1):
        """
        Calculates the integral of value over the time period start_time:end_time

        :param start_time: Beginning of the integration period
        :param end_time: End of the integration period
        :param unit_factor: A factor that may be used to ensure that the units of value are consistent with the units of
            time. For example, if time is measured in days and emissions are measured in g/s, a conversion factor of
            3600 * 24 should be used to convert gram/second*days to grams.
        :return: The integrated value
        """
        tv = np.array(self.time_value)
        ti, vals = tv[:, 0], tv[:, 1]
        cond = (ti >= start_time) & (ti < end_time)
        ti = ti[cond]
        ti = np.concatenate((np.array([start_time]), ti, np.array([end_time])))
        delta_t = ti[1:] - ti[:-1]
        vals = vals[cond]
        index_at_starttime = np.where(np.array(self.time_value)[:, 0] <= start_time)[0]
        if len(index_at_starttime) == 0:
            raise ValueError("ResultContinuous.value is undefined at the specified start time.")

        initial_val = self.time_value[np.max(index_at_starttime)][1]
        vals = np.concatenate((np.array([initial_val]), vals))
        # compute the time interval between consecutive values
        integral = np.sum(vals * delta_t) * unit_factor
        return integral