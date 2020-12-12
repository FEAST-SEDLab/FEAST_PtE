import numpy as np
import feast
import feast.EmissionSimModules.infrastructure_classes as ic
from feast.EmissionSimModules import simulation_classes as sc
import feast.EmissionSimModules.emission_class_functions as lcf
import os
import pickle
from Tests.test_helper import basic_gas_field


def test_component():
    comp = ic.Component(
        emission_production_rate=1e-5,
        dist_type='bootstrap',
        emission_data_path='../ExampleData/DataObjectInstances/production_emissions.p',
        repair_cost_path='../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        base_reparable=True
    )
    if comp.emission_production_rate != 1e-5:
        raise ValueError("Component() is not storing the emission production rate correctly")
    if np.abs(comp.emission_per_comp - 0.0025) > 0.0001:
        raise ValueError("Component() is not computing emission_per_comp from the input emission data correctly")
    if np.abs(comp.null_repair_rate - comp.emission_production_rate / comp.emission_per_comp) > 1e-5:
        raise ValueError("Component() is not computing the correct null_repair_rate for the steady state assumption")


def test_gas_field():
    gf = basic_gas_field()
    if gf.n_sites != 100:
        raise ValueError("gas_field.__init__() is defining gf.n_sites incorrectly.")
    if gf.n_comps != 10000:
        raise ValueError("gas_field.__init__() is not defining gf.n_comps correctly")
    if gf.sites['basic pad']['parameters'].site_inds != [0, 100]:
        raise ValueError("gas_field.__init__() is not defining gf.site_inds correctly")
    comp = gf.sites['basic pad']['parameters'].comp_dict['Fugitive']['parameters']
    new_emissions = gf.emissions.emissions[gf.emissions.emissions.start_time > 0]
    if len(new_emissions.flux) > comp.emission_production_rate * 2 * gf.n_comps * 5 or \
            len(new_emissions.flux) == 0:
        raise ValueError("gas_field.__init__() is not generating new leaks as expected")

    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=2)
    if gf.met['solar intensity'][10] != 226:
        raise ValueError("gas_field.met_data_maker is not loading TMY data correctly.")

    timeobj.current_time = 10 / 24
    timeobj.delta_t = 1/24
    met_dat = gf.get_met(timeobj, ['wind speed', 'solar intensity'])
    if met_dat['wind speed'] != 2.6 or met_dat['solar intensity'] != 226:
        raise ValueError("gas_field.get_met not returning the right values")
    timeobj.current_time = 365
    timeobj.delta_t = 1
    met_dat = gf.get_met(timeobj, ['wind speed', 'solar Intensity'],
                                interp_modes=['max', 'min'],
                                ophrs={'begin': 9, 'end': 17})
    if met_dat['solar intensity'] != 0:
        raise ValueError("gas_field.get_met not returning the correct values when interp mode is min")

    if met_dat['wind speed'] != 4.1:
        raise ValueError("gas_field.get_met not returning the correct values when max is specified")
    wd = gf.get_met(timeobj, ['wind direction'], interp_modes=['mean'], ophrs={'begin': 10, 'end': 11})
    gf.met_data_maker(100)
    wd2 = gf.get_met(timeobj, ['wind direction'], interp_modes=['mean'], ophrs={'begin': 10, 'end': 11})
    if wd2 == wd:
        raise ValueError("gas_field.met_data_maker is not changing the start hour correctly")



def test_gasfield_leak_maker():
    gf = basic_gas_field()
    new_leaks = lcf.Emission()
    time = sc.Time()
    time.current_time = 10
    np.random.seed(0)
    gf.emission_maker(100, new_leaks, 'Fugitive', 0, time, gf.sites['basic pad']['parameters'])
    for f in new_leaks.emissions.flux:
        if f not in gf.comp_dict['Fugitive'].emission_params.leak_sizes['All']:
            raise ValueError("unexpected flux value")
    if np.abs(np.mean(new_leaks.emissions.end_time) - 1 / gf.comp_dict['Fugitive'].null_repair_rate) > 100:
        # Note: this is a probabilistic test, but in 1e7 iterations the test never failed randomly.
        raise ValueError("unexpected endtime distribution")
    return None


def test_bootstrap_emission_maker():
    comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
        name='Fugitive emitters',
        emission_data_path='../ExampleData/DataObjectInstances/production_emissions.p',
        repair_cost_path='../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
        name='basic pad',
        comp_dict={
            'Fugitive': {'number': 450, 'parameters': comp_fug},
        }
    )
    basicpad.site_inds = [0, 10]
    time = sc.Time()
    time.current_time = 10
    n_leaks_in = 1000
    np.random.seed(0)
    leak = lcf.bootstrap_emission_maker(n_leaks_in, 'Fugitive', basicpad, time)
    mean_leak = np.mean(comp_fug.emission_params.leak_sizes['All'])
    if len(leak.emissions.flux) != n_leaks_in:
        raise ValueError("bootstrap_emission_maker is not creating the correct number of emissions")
    if np.abs(np.sum(leak.emissions.flux) / n_leaks_in - mean_leak) / mean_leak > 0.1:
        raise ValueError("bootstrap_emission_maker is not returning the expected mean emissions size")


def test_gasfield_emission_size_maker():
    comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
        name='Fugitive emitters',
        emission_data_path='../ExampleData/DataObjectInstances/production_emissions.p',
        repair_cost_path='../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 450 * 2, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1 / 24, end_time=1)
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    _ = gas_field.emission_size_maker(timeobj)
    return


def test_emission_obj():
    dat_test = feast.input_data_classes.LeakData(notes='test')
    leak_data = {'All': np.array([1, 4, 2, 3, 2, 4])}
    well_counts = {'All': 1}  # Number of wells in the study
    comp_counts = {'All': 600}  # Assumed components per well
    dat_test.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)
    file_out = 'temp_dat.p'
    with open(file_out, 'wb') as f:
        pickle.dump(dat_test, f)
    comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
        name='Fugitive emitters',
        emission_data_path='temp_dat.p',
        repair_cost_path='../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        emission_per_comp=0.00231,
        emission_production_rate=5.4 / 650 / 365
    )
    lsm, _, _, _ = feast.EmissionSimModules.emission_class_functions.emission_objects_generator('bootstrap', file_out)
    basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
        name='basic pad',
        comp_dict={
            # ------ The number of components is proportional to the number of wells, ind2
            'Fugitive': {'number': 600, 'parameters': comp_fug}
        },
    )
    basicpad.site_inds = [0, 1]
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1 / 24, end_time=3 * 365)
    leak = lsm(6, 'Fugitive', basicpad, timeobj)
    for f in leak.emissions.flux:
        if f not in leak_data['All']:
            raise ValueError('emission_objects_generator is not using the correct input emission data file')
    os.remove(file_out)


def test_emission_class():
    leak_specs = {
        'flux': [1, 2, 3, 4, 0, 3, 2],
        'site_index': [1, 2, 3, 4, 1, 2, 0],
        'comp_index': [1, 356, 20, 478, 233, 5, 530],
        'repair_cost': [1000, 1000, 1000, 1000, 1000, 1000, 30],
        'reparable': [True, False, True, True, True, True, False],
        'end_time': [100, np.inf, 1000, 300, 456, 762, 3],
        'start_time': [0, 0, 0, 0, 0, 0, 0],
    }
    leak = feast.EmissionSimModules.emission_class_functions.Emission(**leak_specs)
    if len(leak.emissions.flux) != 7:
        raise ValueError("emission_class_functions.Emission is not initializing the array correctly")
    if leak.emissions.comp_index[1] != 356:
        raise ValueError("emission_class_functions.Emission is not initializing to input values correctly")
    time = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=1001, current_time=400)
    em = leak.get_current_emissions(time)
    if np.any(em.flux != np.array([2, 3, 0, 3])):
        raise ValueError("emission_class_functions.Emission.get_current_emissions is not returning the correct "
                         "emissions.")
    leak2 = feast.EmissionSimModules.emission_class_functions.Emission(**leak_specs)
    leak.extend(leak2)
    if len(leak.emissions.comp_index) != 14:
        raise ValueError("emission_class_functions.Emission.extend is not concatenating emission data frames correctly")

    em = leak2.get_emissions_in_range(300, 764)
    if np.any(em.flux != np.array([2, 3, 4, 0, 3])):
        raise ValueError("emission_class_functions.Emission.get_emissions_in_range is not returning the correct "
                         "emissions")
    em = leak2.get_emissions_in_range(300, 764, reparable=True)
    if np.any(em.flux != np.array([3, 4, 0, 3])):
        raise ValueError("emission_class_functions.Emission.get_emissions_in_range is not returning the correct "
                         "emissions when a reparable condition is set")
    if leak2.em_rate_in_range(0, 1) != np.sum(leak2.emissions.flux):
        raise ValueError("emission_class_functions.Emission.sum_emissions_in_range is not returning the correct value")

    if leak2.em_rate_in_range(99, 101) != (np.sum(np.array(leak_specs['flux'][1:6])) * 2 + leak_specs['flux'][0]) / 2:
        raise ValueError("emission_class_functions.Emission.sum_emissions_in_range is no returning the correct value")

    leak2.emissions.end_time = np.ones(len(leak2.emissions.flux))
    if leak2.em_rate_in_range(99, 101) != 0:
        raise ValueError("emission_class_functions.Emission.sum_emissions_in_range is failing to return zero")


test_component()

test_gas_field()

test_gasfield_leak_maker()

test_bootstrap_emission_maker()

test_gasfield_emission_size_maker()

test_emission_obj()

test_emission_class()

print("Successfully completed emission tests.")
