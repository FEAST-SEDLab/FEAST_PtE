"""
These functions can be called by DetectionModules to aid in calculations
"""


def replacement_cap(time, tech):
    """
    Enters the replacement value of the technology into tech.capital at the end of its lifetime.
    Inputs:
        time        Object defining time parameters for the simulation
        tech        LDAR technology object
    """
    lifetimes = 0
    while lifetimes < time.end_time:
        tech.capital[max(1, round(lifetimes / time.delta_t))] = tech.capital_0
        lifetimes += tech.lifetime
