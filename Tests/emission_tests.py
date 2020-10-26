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
    if gf.new_emissions.n_em > comp.emission_production_rate * 2 * gf.n_comps * 5 or gf.new_emissions.n_em == 0:
        raise ValueError("gas_field.__init__() is not generating new leaks as expected")

    np.random.seed(0)
    n_sites = 100
    site_dict = {}
    comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
        name='Fugitive emitters',
        emission_data_path='../ExampleData/DataObjectInstances/production_emissions.p',
        repair_cost_path='../ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive': {'number': 100, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=2)
    initial_leaks = feast.EmissionSimModules.emission_class_functions.Emission(
        flux=np.ones(100), site_index=np.random.randint(0, n_sites, 100),
        comp_index=np.random.randint(0, 100, 100), endtime=np.infty, repair_cost=np.ones(100) * 2
    )
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
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
                                ophrs={'begin': 9, 'end': 17})
    if met_dat['solar intensity'] != 0:
        raise ValueError("gas_field.get_met not returning the correct values when interp mode is min")

    if met_dat['wind speed'] != 4.1:
        raise ValueError("gas_field.get_met not returning the correct values when max is specified")


def test_gasfield_leak_maker():
    gf = basic_gas_field()
    new_leaks = lcf.Emission()
    time = sc.Time()
    time.current_time = 10
    np.random.seed(0)
    gf.emission_maker(100, new_leaks, 'Fugitive', 0, time, gf.sites['basic pad']['parameters'])
    for f in new_leaks.flux:
        if f not in gf.comp_dict['Fugitive'].emission_params.leak_sizes['All']:
            raise ValueError("unexpected flux value")
    if np.abs(np.mean(new_leaks.endtime) - 1 / gf.comp_dict['Fugitive'].null_repair_rate) > 100:
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
    if len(leak.flux) != n_leaks_in:
        raise ValueError("bootstrap_emission_maker is not creating the correct number of emissions")
    if np.abs(np.sum(leak.flux) / n_leaks_in - mean_leak) / mean_leak > 0.1:
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
    for f in leak.flux:
        if f not in leak_data['All']:
            raise ValueError('emission_objects_generator is not using the correct input emission data file')
    os.remove(file_out)


def test_emission_class():
    fl = [1, 2, 3, 4, 0, 3, 2]
    si = [1, 2, 3, 4, 1, 2, 0]
    ci = [1, 356, 20, 478, 233, 5, 530]
    rc = [1000, 1000, 1000, 1000, 1000, 1000, 30]
    rep = [True, False, True, True, True, True, False]
    et = [100, np.inf, 1000, 300, 456, 762, 3]
    leak = feast.EmissionSimModules.emission_class_functions.Emission(
        flux=fl,
        site_index=si,
        comp_index=ci,
        repair_cost=rc,
        reparable=rep,
        endtime=et,
        capacity=100,
    )
    if leak.n_em != 7:
        raise ValueError("leak_class_function.Leak is not initializing the number of leaks correctly")
    leak2 = feast.EmissionSimModules.emission_class_functions.Emission(
        flux=np.ones(100),
        site_index=np.ones(100),
        comp_index=np.ones(100),
        repair_cost=np.ones(100),
        reparable=True,
    )
    if len(leak2.flux) != 100:
        raise ValueError("leak_class_function.Leak is not initializing the array capacity correctly")
    leak.extend(leak2)
    if leak.n_em != 107:
        raise ValueError("leak_class_function.Leak.extend is not updating n_leaks correctly")
    if np.min(leak.flux) != 0:
        raise ValueError("leak_class_function.Leak.extend is not preserving the minimum leak flux")
    if np.sum(leak.flux == 0) > 1:
        raise ValueError("leak_class_function.Leak.extend is leaving too many zero values")
    if np.any(leak.flux != np.append(fl, leak2.flux)):
        raise ValueError("leak_class_function.Leak.extend is not updating flux correctly")
    if np.any(leak.site_index != np.append(si, leak2.site_index)):
        raise ValueError("leak_class_function.Leak.extend is not updating site index correctly")
    if np.any(leak.comp_index != np.append(ci, leak2.comp_index)):
        raise ValueError("leak_class_function.Leak.extend is not updating component index correctly")
    if np.any(leak.repair_cost != np.append(rc, leak2.repair_cost)):
        raise ValueError("leak_class_function.Leak.extend is not updating repair cost correctly")
    if np.any(leak.reparable != np.append(rep, leak2.reparable)):
        raise ValueError("leak_class_function.Leak.extend is not updating reparable correctly")
    if np.any(leak.endtime != np.append(et, leak2.endtime)):
        raise ValueError("leak_class_function.Leak.extend is not updating end time correctly")
    leak.delete_leaks(2)
    if leak.flux[2] != 0:
        raise ValueError("leak_class_function.Leak.delete is not setting flux to zero")
    if leak.site_index[2] != si[2]:
        raise ValueError("leak_class_function.Leak.delete is changing the site index of an emission when it should not")
    if leak.n_em != 107:
        raise ValueError("leak_class_function.Leak.delete is changing n_leaks when it should not")
    leak.clear_zeros()
    if leak.n_em != 105:  # Started with 107, but one with 0 flux, then removed one with leak.delete
        raise ValueError("leak_class_function.Leak.clear_zeros is not updating n_leaks correctly")
    if len(leak.flux) != 105:
        raise ValueError("leak_class_function.Leak.clear_zeros is not deleting elements correctly")
    if np.min(leak.flux[:leak.n_em]) <= 0:
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements flux")
    if len(leak.flux) != len(leak.site_index):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from site_index")
    if len(leak.flux) != len(leak.comp_index):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from comp_index")
    if len(leak.flux) != len(leak.reparable):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from reparable")
    if len(leak.flux) != len(leak.endtime):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from endtime")
    if len(leak.flux) != len(leak.repair_cost):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from repair_cost")
    leak.sort_by_site()
    if leak.flux[0] != 2:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting flux values correctly")
    if leak.comp_index[0] != ci[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting component index values correctly")
    if leak.repair_cost[0] != rc[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting repair cost values correctly")
    if leak.reparable[0] != rep[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting reparable values correctly")
    if leak.endtime[0] != et[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting end time values correctly")


test_component()

test_gas_field()

test_gasfield_leak_maker()

test_bootstrap_emission_maker()

test_gasfield_emission_size_maker()

test_emission_obj()

test_emission_class()

print("Successfully completed emission tests.")
