import numpy as np
import copy
import feast
import os
import pickle
import time as ti
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.emission_class_functions as ecf
from Tests.test_helper import basic_gas_field


def test_gas_field():
    gf = basic_gas_field()
    if gf.n_sites != 100:
        raise ValueError("gas_field.__init__() is defining gf.n_sites incorrectly.")
    if gf.n_comps != 10000:
        raise ValueError("gas_field.__init__() is not defining gf.n_comps correctly")
    if gf.sites['basic pad']['parameters'].site_inds != [0, 100]:
        raise ValueError("gas_field.__init__() is not defining gf.site_inds correctly")
    comp = gf.sites['basic pad']['parameters'].comp_dict['Fugitive']['parameters']
    if gf.new_leaks.n_leaks > comp.emission_production_rate * 2 * gf.n_comps * 5 or gf.new_leaks.n_leaks == 0:
        raise ValueError("gas_field.__init__() is not generating new leaks as expected")

    np.random.seed(0)
    n_sites = 100
    site_dict = {}
    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_data_path='production_emissions.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive': {'number': 100, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1, end_time=2)
    initial_leaks = feast.GeneralClassesFunctions.emission_class_functions.Emission(
        flux=np.ones(100), site_index=np.random.randint(0, n_sites, 100),
        comp_index=np.random.randint(0, 100, 100), endtime=np.infty, repair_cost=np.ones(100) * 2
    )
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj,
        initial_emissions=initial_leaks,
        met_data_path='TMY-DataExample.csv'
    )
    if gas_field.met['solar intensity'][10] != 226:
        raise ValueError("gas_field.met_data_maker is not loading TMY data correctly.")

    timeobj.current_time = 10 / 24
    timeobj.delta_t = 1/24
    met_dat = gas_field.get_met(timeobj, ['wind speed', 'solar intensity'])
    if met_dat['wind speed'] != 2.6 or met_dat['solar intensity'] != 226:
        raise ValueError("gas_field.get_met not returning the right values")
    timeobj.current_time = 365
    timeobj.delta_t = 1
    met_dat = gas_field.get_met(timeobj, ['wind speed', 'solar Intensity'],
                                interp_modes=['max', 'min'],
                                op_hrs={'begin': 900, 'end': 1700})
    if met_dat['solar intensity'] != 0:
        raise ValueError("gas_field.get_met not returning the correct values when interp mode is min")

    if met_dat['wind speed'] != 4.1:
        raise ValueError("gas_field.get_met not returning the correct values when max is specified")




test_gas_field()
