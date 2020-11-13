from feast import MEET_1_importer as meet
import numpy as np


def test_load_gas_comp_file():
    dat, delta_t = meet.load_gas_comp_file('gas_composition.csv')
    if delta_t != 3600:
        raise ValueError("load_gas_comp_file() is not returning the correct time resolution")


def test_gas_comp_to_dict():
    gas_comp, meet_delta_t = meet.load_gas_comp_file('gas_composition.csv')
    dat_dict = meet.gas_comp_to_dict(gas_comp, meet_delta_t)
    if len(np.unique(dat_dict['loc'])) != 10:
        raise ValueError('gas_comp_to_dict is not registering the correct number of sites')
    if len(np.unique(dat_dict['em_id'])) != 293:
        raise ValueError('gas_comp_to_dict is not extracting the correct number of emitters')


def test_gc_dat_to_gas_field():
    gas_comp, meet_delta_t = meet.load_gas_comp_file('gas_composition.csv')
    dat_dict = meet.gas_comp_to_dict(gas_comp, meet_delta_t)
    comps_per_site = {}
    for ind in np.unique(dat_dict['loc']):
        comps_per_site[ind] = 1000
    rep_cost_path = '../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p'
    time_obj, gas_field = meet.gc_dat_to_gas_field(1 / 24, 1, comps_per_site, rep_cost_path, **dat_dict)
    if time_obj.delta_t != 1/24 or time_obj.end_time != 1:
        raise ValueError("gc_dat_to_gas_field() is not returning the correct Time object")
    if gas_field.n_sites != 10:
        raise ValueError("gc_dat_to_gas_field() is not a gas_field with the correct n_sites")
    if len(np.unique(gas_field.emissions.emissions.index)) != 293:
        raise ValueError('gc_dat_to_gas_field() is not returning the right number of emission IDs')
    em = gas_field.emissions.emissions
    max_em = em.loc[em['flux'].idxmax(), :].iloc[0]
    if np.abs(max_em.flux - 2758.8) > 0.1:
        raise ValueError('gc_dat_to_gas_field() is not returning the correct emission rates')
    if max_em.start_time != 9 / 24:
        raise ValueError('gc_dat_to_gas_field() is not returning the correct start times')
    if max_em.end_time != 10 / 24:
        raise ValueError('gc_dat_to_gas_field() is not returning the correct end times')
    if max_em.reparable:
        raise ValueError('gc_dat_to_gas_field() is not returning the correct emission types')
    if max_em.site_index != 6:
        raise ValueError('gc_dat_to_gas_field() is not returning the correct site indexes')


test_load_gas_comp_file()
test_gas_comp_to_dict()
test_gc_dat_to_gas_field()

print("Successfully completed meet_reader tests.")