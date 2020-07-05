import numpy as np
import feast
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.emission_class_functions as lcf
import feast.GeneralClassesFunctions.results_analysis_functions as raf
from feast import DetectionModules as Dm
import os
import pickle


def basic_gas_field():
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
        initial_emissions=initial_leaks
    )
    return gas_field
