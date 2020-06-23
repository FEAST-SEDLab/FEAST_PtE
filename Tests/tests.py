import numpy as np
import feast
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.emission_class_functions as lcf
import feast.GeneralClassesFunctions.results_analysis_functions as raf
from feast import DetectionModules as Dm
import os
import pickle
import time as ti


def test_gasfield_leak_maker():
    gf = sc.GasField()
    new_leaks = lcf.Emission()
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
    repair_cost = np.zeros(100)
    gf.initial_emissions = lcf.Emission(flux=flux, endtime=end_time, repair_cost=repair_cost)
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
            tech_dict={}, dir_out='ResultsTemp', display_status=False
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
    with open(file_out, 'wb') as f:
        pickle.dump(dat_test, f)
    comp_fug = feast.GeneralClassesFunctions.simulation_classes.Component(
        name='Fugitive emitters',
        emission_data_path='temp_dat.p',
        emission_per_comp=0.00231,
        emission_production_rate=5.4 / 650 / 365
    )
    lsm, _, _, _ = feast.GeneralClassesFunctions.emission_class_functions.leak_objects_generator('bootstrap', file_out)
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
    null_npv, emissions, costs, techs = raf.results_analysis('ResultsTemp')
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


def test_npv_calculator():
    file_out = 'temp_emissions.p'
    rep_file_out = 'temp_rep_costs.p'
    em_array = np.ones(100)
    emissions = feast.InputData.input_data_classes.LeakData()
    # The dict structure allows for multiple types of detection methods used in the study
    leak_data = {'All': em_array}
    well_counts = {'All': 100}  # Number of wells in the study
    comp_counts = {'All': 100}  # Assumed components per well
    emissions.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)
    with open(file_out, 'wb') as f:
        pickle.dump(emissions, f)
    repair_out = feast.InputData.input_data_classes.RepairData()
    repair_out.define_data(repair_costs=np.ones(5) * 2)
    with open(rep_file_out, 'wb') as f:
        pickle.dump(repair_out, f)

    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_data_path=file_out,
        emission_per_comp=0.1,
        repair_cost_path=rep_file_out,
        emission_production_rate=0
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 100, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1, end_time=2)
    initial_leaks = feast.GeneralClassesFunctions.emission_class_functions.Emission(
        flux=np.ones(1000), site_index=np.random.randint(0, n_sites, 1000),
        comp_index=np.random.randint(0, 100, 1000), endtime=np.infty, repair_cost=np.ones(1000) * 2
    )
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj,
        initial_leaks=initial_leaks
    )
    tech_dict = {'TechDetect': feast.DetectionModules.tech_detect.TechDetect(timeobj, gas_field, survey_speed=10000)}
    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field, dir_out='ResultsTemp', display_status=False, tech_dict=tech_dict
    )
    npv = feast.GeneralClassesFunctions.results_analysis_functions.npv_calculator('ResultsTemp/realization0.p')
    if npv['Repair'] != 2000:
        raise ValueError("npv_calculator not returning the expected repair cost")
    if np.abs(npv['Gas'] - 34560) > 100:  # 34560 = gas_value * 1000 g/s * 2 days * 3600 * 24 sec/day--no discount rate.
        raise ValueError("npv_calculator not returning the expected value of gas saved")
    os.remove(file_out)
    os.remove(rep_file_out)
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    try:
        os.rmdir('ResultsTemp')
    except PermissionError:
        # If there is an automated syncing process, a short pause may be necessary before removing "ResultsTemp"
        ti.sleep(5)
        os.rmdir('ResultsTemp')


def test_tiered_detect_find_cost():
    comp_fug = sc.Component(
        name='Fugitive emitters',
        emission_per_comp=0.1,
        emission_production_rate=0
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
        # Simulates two wells, one tank, total components=11302
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 100, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1, end_time=2)
    initial_leaks = feast.GeneralClassesFunctions.emission_class_functions.Emission(
        flux=np.ones(100), site_index=np.ones(100),
        comp_index=np.random.randint(0, 100, 100), endtime=np.infty, repair_cost=np.ones(100) * 2
    )
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj,
        initial_leaks=initial_leaks
    )
    tech_dict = {'TieredDetect': feast.DetectionModules.tiered_detect.TieredDetect(
                     timeobj,
                     gas_field,
                     sites_per_day=1000,
                     site_cost=1,
                     secondary_comps_hr=100,
                     labor=100,
                     ophrs={'begin': 800, "end": 1700}
                 )}
    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field, dir_out='ResultsTemp', display_status=False, tech_dict=tech_dict
    )
    npv = feast.GeneralClassesFunctions.results_analysis_functions.npv_calculator('ResultsTemp/realization0.p')
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    try:
        os.rmdir('ResultsTemp')
    except PermissionError:
        # If there is an automated syncing process, a short pause may be necessary before removing "ResultsTemp"
        ti.sleep(5)
        os.rmdir('ResultsTemp')
    if npv['Finding'][0] != 200:
        raise ValueError("Finding costs not calculated correctly for tiered detect method.")


def test_leak_class():
    fl = [1, 2, 3, 4, 0, 3, 2]
    si = [1, 2, 3, 4, 1, 2, 0]
    ci = [1, 356, 20, 478, 233, 5, 530]
    rc = [1000, 1000, 1000, 1000, 1000, 1000, 30]
    ld = [0, 0, 0, 1, 0, 0, 1]
    rep = [True, False, True, True, True, True, False]
    et = [100, np.inf, 1000, 300, 456, 762, 3]
    leak = feast.GeneralClassesFunctions.emission_class_functions.Emission(
        flux=fl,
        site_index=si,
        comp_index=ci,
        repair_cost=rc,
        leaks_detected=ld,
        reparable=rep,
        endtime=et,
        capacity=100,
    )
    if leak.n_leaks != 7:
        raise ValueError("leak_class_function.Leak is not initializing the number of leaks correctly")
    leak2 = feast.GeneralClassesFunctions.emission_class_functions.Emission(
        flux=np.ones(100),
        site_index=np.ones(100),
        comp_index=np.ones(100),
        repair_cost=np.ones(100),
        reparable=True,
    )
    if len(leak2.flux) != 100:
        raise ValueError("leak_class_function.Leak is not initializing the array capacity correctly")
    leak.extend(leak2)
    if leak.n_leaks != 107:
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
    if np.any(leak.leaks_detected != np.append(ld, leak2.leaks_detected)):
        raise ValueError("leak_class_function.Leak.extend is not updating leaks detected correctly")
    if np.any(leak.reparable != np.append(rep, leak2.reparable)):
        raise ValueError("leak_class_function.Leak.extend is not updating reparable correctly")
    if np.any(leak.endtime != np.append(et, leak2.endtime)):
        raise ValueError("leak_class_function.Leak.extend is not updating end time correctly")
    leak.delete_leaks(2)
    if leak.flux[2] != 0:
        raise ValueError("leak_class_function.Leak.delete is not setting flux to zero")
    if leak.site_index[2] != si[2]:
        raise ValueError("leak_class_function.Leak.delete is changing the site index of an emission when it should not")
    if leak.n_leaks != 107:
        raise ValueError("leak_class_function.Leak.delete is changing n_leaks when it should not")
    leak.clear_zeros()
    if leak.n_leaks != 105:  # Started with 107, but one with 0 flux, then removed one with leak.delete
        raise ValueError("leak_class_function.Leak.clear_zeros is not updating n_leaks correctly")
    if len(leak.flux) != 105:
        raise ValueError("leak_class_function.Leak.clear_zeros is not deleting elements correctly")
    if np.min(leak.flux[:leak.n_leaks]) <= 0:
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
    if len(leak.flux) != len(leak.leaks_detected):
        raise ValueError("leak_class_function.Leak.clear_zeros is not removing the correct elements from "
                         "leaks_detected")
    leak.sort_by_site()
    if leak.flux[0] != 2:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting flux values correctly")
    if leak.comp_index[0] != ci[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting component index values correctly")
    if leak.repair_cost[0] != rc[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting repair cost values correctly")
    if leak.reparable[0] != rep[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting reparable values correctly")
    if leak.leaks_detected[0] != ld[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting leaks_detected values correctly")
    if leak.endtime[0] != et[-1]:
        raise ValueError("leak_class_function.Leak.sort_by_site is not sorting end time values correctly")


test_gasfield_leak_maker()

test_null_repair()

test_alvarez_colorado_em_dist()

test_bootstrap_leak_maker()

test_gasfield_leak_size_maker()

test_field_simulation()

test_leak_obj()

test_results_analysis()

test_npv_calculator()

test_tiered_detect_find_cost()

test_leak_class()

print("Successfully completed all tests")
