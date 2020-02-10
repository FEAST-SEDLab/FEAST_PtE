"""
This module stores various methods than can be called to simulate site level emissions
"""
import numpy as np


def alvarez_colorado_em_dist(prod):
    """
        Generates an emission for each site in prod according to a lognormal distribution.
        See Eq. 1 in the Supporting Materials for Alvarez et al 2018.
        :input prod: an array of site production values. Must all be >0.
        :return em: an array of site emissions [g/s]
    """
    if np.min(prod) <= 0:
        raise ValueError('All production values must be greater than zero to use alvarez_colorado_em_dist')
    a = 2.6
    b = -2.2
    c = 0.2
    theta1 = 0.6
    theta2 = 1.4
    sigma = 1.3
    xp = np.log(prod)
    mup = a + b * xp ** theta1 + c * xp ** theta2
    em = np.random.lognormal(mup, sigma)
    return em * 1000 / 3600


def start_at_zero(prod):
    """
        Sets all site emissions to 0
        :input prod: an array of site production values. Must all be >0.
        :return em: an array of site emissions [g/s]
    """
    em = np.zeros(len(prod))
    return em
