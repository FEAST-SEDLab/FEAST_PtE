import numpy as np
import feast
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.leak_class_functions as lcf
import feast.GeneralClassesFunctions.results_analysis_functions as raf
from feast import DetectionModules as Dm
import os
import pickle


def test_gasfield_leak_maker():
    gf = sc.GasField()
    new_leaks = lcf.Leak()
    time = sc.Time()
    time.current_time = 10
    np.random.seed(0)
    gf.leak_maker(10, new_leaks, 'default', 0, time, gf.sites['default']['parameters'])
    expected = np.array([0.02689507, 0.00775865, 0.00228167, 0.01050579, 0.00157465,
                         0.0149844, 0.04953858, 0.49138137, 0.09214868, 0.04406508])
    if np.max(np.abs(expected - new_leaks.flux)) > 1e-8:
        raise ValueError("unexpected flux value")
    expected_endtimes = np.array([
        391.04863847, 165.94777145, 225.15696612, 14.85666044,
        256.14930494, 252.46640114, 255.67999985, 746.85338952,
        303.1982008, 124.06931841
    ])
    if np.max(np.abs(expected_endtimes - new_leaks.endtime)) > 1e-6:
        raise ValueError("unexpected endtime")
    return None


def test_null_repair():
    gf = sc.GasField()
    timeobj = sc.Time(delta_t=1 / 24, end_time=50)
    timeobj.current_time = 10
    np.random.seed(0)
    flux = np.random.uniform(0.1, 1, 100)
    end_time = np.random.exponential(10, 100)
    gf.initial_leaks = lcf.Leak(flux=flux, endtime=end_time)
    null = Dm.null.Null(timeobj, gf)
    null.null_detection(timeobj, gf)
    if np.sum(null.leaks.flux > 0) != 39:
        raise ValueError("unexpected number of leaks persisting")
    return None


def test_alvarez_colorado_em_dist():
    np.random.seed(0)
    em_out = feast.GeneralClassesFunctions.site_emission_methods.alvarez_colorado_em_dist(np.array([11, 121]))
    if np.max(np.abs(em_out - np.array([1.77656241, 0.13525278]))) > 1e-6:
        raise ValueError("Unexpected emission from test_alvarez_colorado_em_dist")


def test_emissions_enforcer_low():
    pad = sc.Site(site_em_dist=feast.GeneralClassesFunctions.site_emission_methods.alvarez_colorado_em_dist)
    gf = sc.GasField(sites={'site0': {'number': 1, 'parameters': pad}})
    gf.initial_leaks = lcf.Leak(capacity=1)
    if np.sum(gf.initial_leaks.flux) != 0:
        raise ValueError("lcf.Leak intializing non-zero emissions")
    gf.site_emissions_enforcer(gf.sites['site0']['parameters'])
    if np.sum(gf.initial_leaks.flux) <= 0:
        raise ValueError("site_emissions_enforcer is not generating emissions when required")
    if gf.initial_leaks.site_index != 0:
        raise ValueError("site_emissions_enforcer does not set the site index correctly")
    if gf.initial_leaks.comp_index != -1:
        raise ValueError("site_emissions_enforcer does not set the component index correctly")


def test_emissions_enforcer_high():
    pad = sc.Site(site_em_dist=feast.GeneralClassesFunctions.site_emission_methods.alvarez_colorado_em_dist)
    gf = sc.GasField(sites={'site0': {'number': 1, 'parameters': pad}})
    gf.initial_leaks = lcf.Leak(capacity=10,
                                flux=np.array([1000, 100, 10, 1, 0.1]),
                                site_index=np.array([0, 0, 0, 0, 0]),
                                comp_index=np.array([1, 2, 3, 4, 5]))
    if np.sum(gf.initial_leaks.flux) != 1111.1:
        raise ValueError("lcf.Leak initializing unexpected emissions")
    np.random.seed(0)
    gf.site_emissions_enforcer(gf.sites['site0']['parameters'])
    if np.max(gf.initial_leaks.flux[1:]) > 0 or gf.initial_leaks.flux[0] >= 1000:
        raise(ValueError("emissions_enforcer is not removing emissions as expected"))


def test_emissions_enforcer_no_repairable():
    pad = sc.Site(site_em_dist=feast.GeneralClassesFunctions.site_emission_methods.alvarez_colorado_em_dist)
    gf = sc.GasField(sites={'site0': {'number': 1, 'parameters': pad}})
    gf.initial_leaks = lcf.Leak(capacity=10,
                                flux=np.array([1000, 100, 10, 1, 0.1]),
                                site_index=np.array([0, 0, 0, 0, 0]),
                                comp_index=np.array([1, 2, 3, 4, 5]),
                                reparable=np.array([False, False, False, False, False]))
    if np.sum(gf.initial_leaks.flux) != 1111.1:
        raise ValueError("lcf.Leak initializing unexpected emissions")
    np.random.seed(0)
    gf.site_emissions_enforcer(gf.sites['site0']['parameters'])


def test_bootstrap_leak_maker():
    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_data_path='production_emissions.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
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
    leak = lcf.bootstrap_leak_maker(n_leaks_in, 'Fugitive', basicpad, time)
    mean_leak = np.mean(comp_fug.leak_params.leak_sizes['All'])
    if len(leak.flux) != n_leaks_in:
        raise ValueError("bootstrap_leak_maker is not creating the correct number of emissions")
    if np.abs(np.sum(leak.flux) / n_leaks_in - mean_leak) / mean_leak > 0.1:
        raise ValueError("bootstrap_leak_maker is not returning the expected mean emissions size")


def test_gasfield_leak_size_maker():
    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_data_path='production_emissions.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 450 * 2, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1 / 24, end_time=1)
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    _ = gas_field.leak_size_maker(timeobj)
    return


def test_field_simulation():
    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_data_path='production_emissions.p',
        emission_per_comp=0.0026,
        emission_production_rate=5.4 / 650 / 365
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 450 * 2, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1 / 24, end_time=100)
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    feast.field_simulation.field_simulation(
            time=timeobj, gas_field=gas_field,
            tech_dict={}, dir_out='ResultsTemp', display_status=True
        )
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


def test_leak_obj():
    dat_test = feast.InputData.input_data_classes.LeakData(notes='test')
    leak_data = {'All': np.array([1, 4, 2, 3, 2, 4])}
    well_counts = {'All': 1}  # Number of wells in the study
    comp_counts = {'All': 600}  # Assumed components per well
    dat_test.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)
    file_out = 'temp_dat.p'
    pickle.dump(dat_test, open(file_out, 'wb'))
    comp_fug = feast.GeneralClassesFunctions.simulation_classes.Component(
        name='Fugitive emitters',
        emission_data_path='temp_dat.p',
        emission_per_comp=0.00231,
        emission_production_rate=5.4 / 650 / 365
    )
    lsm, _, _, _ = feast.GeneralClassesFunctions.leak_class_functions.leak_objects_generator('bootstrap', file_out)
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        name='basic pad',
        comp_dict={
            # ------ The number of components is proportional to the number of wells, ind2
            'Fugitive': {'number': 600, 'parameters': comp_fug}
        },
    )
    basicpad.site_inds = [0, 1]
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1 / 24, end_time=3 * 365)
    leak = lsm(6, 'Fugitive', basicpad, timeobj)
    for f in leak.flux:
        if f not in leak_data['All']:
            raise ValueError('leak_objects_generator is not using the correct input emission data file')
    os.remove(file_out)


def test_results_analysis():
    for ind in range(3):
        feast.field_simulation.field_simulation(dir_out='ResultsTemp', display_status=False)
    null_npv, emissions = raf.results_analysis('ResultsTemp')
    if len(null_npv.keys()) != 6:
        raise ValueError("results analysis function returning the wrong number of keys")
    if null_npv['Finding'].shape != (1, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or LDAR programs in "
                         "null npv")
    if emissions.shape != (2, 4001, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or "
                         "LDAR programs in emissions")
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


test_gasfield_leak_maker()

test_null_repair()

test_alvarez_colorado_em_dist()

test_emissions_enforcer_low()

test_emissions_enforcer_high()

test_emissions_enforcer_no_repairable()

test_bootstrap_leak_maker()

test_gasfield_leak_size_maker()

test_field_simulation()

test_leak_obj()

test_results_analysis()

print("Successfully completed all tests")
