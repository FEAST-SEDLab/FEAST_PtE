import numpy as np
import feast
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.leak_class_functions as lcf
from feast import DetectionModules as Dm
import os


def test_gasfield_leak_maker():
    gf = sc.GasField()
    new_leaks = lcf.Leak()
    time = sc.Time()
    time.current_time = 10
    np.random.seed(0)
    gf.leak_maker(10, new_leaks, 'default', 0, time, gf.sites['default']['parameters'])
    expected = np.array([0.001951722536, 0.002732411550, 0.002112817094, 0.000154898614, 0.050311069805, 0.003423259368,
                         0.001375499692, 0.022832055693, 0.018364779668, 0.036927829561])
    if np.max(np.abs(expected - new_leaks.flux)) > 1e-8:
        raise ValueError("unexpected flux value")
    expected_endtimes = np.array([
        1464.853673030148, 605.412672264341, 831.474926906344, 28.542856664766, 949.804487527753, 935.743064465414,
        948.012668427203, 2823.325523307107, 1129.438403135151, 445.519642662429
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
    if np.max(np.abs(gf.initial_leaks.flux - np.array([0.68284285, 0, 0, 0, 0, 0, 0, 0, 0, 0]))) > 1e-6:
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
        leak_data_path='production_emissions.p',
        leaks_per_comp=0.0026,
        leak_production_rate=5.4 / 650 / 365
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
        leak_data_path='production_emissions.p',
        leaks_per_comp=0.0026,
        leak_production_rate=5.4 / 650 / 365
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
        leak_data_path='production_emissions.p',
        leaks_per_comp=0.0026,
        leak_production_rate=5.4 / 650 / 365
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
    os.remove('ResultsTemp/realization0.p')
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

print("Successfully completed all tests")