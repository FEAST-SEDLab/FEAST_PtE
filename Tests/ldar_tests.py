import numpy as np
import copy
import feast
import os
import pickle
import time as ti
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.emission_class_functions as ecf
from feast import DetectionModules as Dm
from Tests.test_helper import basic_gas_field


def test_repair():
    n_em = 10
    flux, capacity, reparable, repair_cost = np.zeros(n_em), 10, np.zeros(n_em, dtype=bool), np.ones(n_em)
    site_index, comp_index = np.zeros(n_em), np.zeros(n_em)
    flux[5:] = 1
    reparable[3: 8] = True
    endtime = np.linspace(0, 9, n_em)
    emission = ecf.Emission(flux, capacity, reparable, site_index, comp_index, repair_cost=repair_cost, endtime=endtime)
    time = sc.Time()
    time.time_index = 3
    time.current_time = 4.5
    detected = np.array([3, 4, 5, 6, 7, 8, 9], dtype=int)
    repair_proc = Dm.repair.Repair(repair_delay=1)
    repair_proc.action(emit_inds=detected)
    repair_proc.repair(time, emission)
    expected = np.array([0., 1., 2., 3., 4., 5., 5.5, 5.5, 8., 9.])
    for ind in range(len(emission.endtime)):
        if emission.endtime[ind] != expected[ind]:
            raise ValueError("DetectionModules.repair.Repair is not adjusting "
                             "emission endtimes correctly at index {:0.0f}".format(ind))


def test_check_time():
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    tech = Dm.comp_detect.CompDetect(
        time,
        survey_interval=50,
        dispatch_object=rep
    )
    if not tech.check_time(time):
        raise ValueError("Check time returning False when it should be true at time 0")
    time = sc.Time(delta_t=0.1, end_time=10, current_time=0)
    tech = Dm.comp_detect.CompDetect(
        time,
        survey_interval=50,
        dispatch_object=rep
    )
    if tech.check_time(time):
        raise ValueError("check_time returning True when it should return False at time 0")


def test_comp_detect():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    find_cost = np.zeros(time.n_timesteps)
    rep = Dm.repair.Repair(repair_delay=0)
    tech = Dm.comp_detect.CompDetect(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    emissions = gas_field.initial_emissions
    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
    tech.detect(time, gas_field, emissions, find_cost)
    if np.max(emissions.site_index[rep.to_repair]) != 13:
        raise ValueError("tech.detect repairing emissions at incorrect sites")
    expected_detected_inds = [54, 55, 75, 66, 87, 62, 58, 84, 5, 25, 43, 71, 13, 79, 97]
    for ind in range(len(rep.to_repair)):
        if expected_detected_inds[ind] != rep.to_repair[ind]:
            raise ValueError("tech.detect not detecting the correct emissions")
    if len(expected_detected_inds) != len(rep.to_repair):
        raise ValueError("tech.detect not detecting the correct emissoins")
    rep.repair(time, emissions)
    if np.max(emissions.endtime[np.array(expected_detected_inds)]) != 0:
        raise ValueError("rep.repair not adjusting the end times correctly")

def test_site_detect():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    find_cost = np.zeros(time.n_timesteps)
    # Test __init__
    tech = Dm.site_detect.SiteDetect(
        time,
        survey_interval=50,
        sites_per_day=100,
        ophrs={'begin': 800, 'end': 1700},
        site_cost=100,
        dispatch_object=Dm.comp_detect.CompDetect(time)
    )
    emissions = gas_field.initial_emissions
    np.random.seed(0)
    # test detect_prob_curve
    detect = tech.detect_prob_curve([0, 1, 2], emissions)
    if detect != np.array([0]):
        raise ValueError("site_detect.detect_prob_curve not returning expected sites.")
    # test sites_surveyed with empty queue
    sites_surveyed = tech.sites_surveyed(time, find_cost)
    if sites_surveyed != []:
        raise ValueError("sites_surveyed returning sites when it should not")
    if find_cost[0] > 0:
        raise ValueError("sites_surveyed updating find_cost when it should not")
    # test detect
    np.random.seed(0)
    tech.sites_to_survey = [0, 1, 2]
    tech.detect(time, gas_field, emissions, np.zeros(time.n_timesteps))
    if tech.dispatch_object.sites_to_survey != [0]:
        raise ValueError("site_detect.detect not updating dispatch object sites to survey correctly")
    # test action and sites_surveyed with
    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
    if tech.sites_to_survey != list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)):
        raise ValueError("action is not updating sites_to_survey as expected")
    # test sites_surveyed with full queue
    sites_surveyed = tech.sites_surveyed(time, find_cost)
    if (sites_surveyed != np.linspace(0, 99, 100, dtype=int)).any():
        raise ValueError("sites_surveyed not identifying the correct sites")
    if find_cost[0] != 10000:
        raise ValueError("sites_surveyed not updating find_cost as expected.")


def test_ldar_program():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    ogi = Dm.comp_detect.CompDetect(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    ogi_no_survey = Dm.comp_detect.CompDetect(
        time,
        survey_interval=None,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    plane_survey = Dm.site_detect.SiteDetect(
        time,
        survey_interval=50,
        sites_per_day=200,
        site_cost=100,
        mu=0.1,
        dispatch_object=ogi_no_survey
    )
    # test __init__
    ogi_survey = Dm.ldar_program.LDARProgram(
        time, copy.deepcopy(gas_field), {'ogi': ogi}, rep
    )
    if len(ogi_survey.find_cost) != 11:
        raise ValueError("find_cost not set to the correct length")
    if np.sum(ogi_survey.emissions.flux) != 100:
        raise ValueError("Unexpected emission rate in LDAR program initialization")

    # test end_emissions
    ogi_survey.emissions.endtime[0] = 0
    ogi_survey.end_emissions(time)
    if ogi_survey.emissions.flux[0] != 0:
        raise ValueError("ldar_program.end_emissions not zeroing emissions as expected")
    # test action
    ogi_survey.action(time, gas_field)
    if np.sum(ogi_survey.emissions.flux) != 84:
        raise ValueError("Unexpected emission rate after LDAR program action")
    # test combined program
    tech_dict = {
        'plane': plane_survey,
        'ogi': ogi_no_survey
    }
    tiered_survey = Dm.ldar_program.LDARProgram(
        time, gas_field, tech_dict, rep
    )
    # test action
    tiered_survey.action(time, gas_field)
    if np.sum(tiered_survey.emissions.flux) != 80:
        raise ValueError("Unexpected emission rate after LDAR program action with tiered survey")


def test_field_simulation():
    gas_field = basic_gas_field()
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1, end_time=2)
    rep = Dm.repair.Repair(repair_delay=0)
    ogi = Dm.comp_detect.CompDetect(
        timeobj,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    ogi_no_survey = Dm.comp_detect.CompDetect(
        timeobj,
        survey_interval=None,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    plane_survey = Dm.site_detect.SiteDetect(
        timeobj,
        survey_interval=50,
        sites_per_day=200,
        site_cost=100,
        mu=0.1,
        dispatch_object=ogi_no_survey
    )
    ogi_survey = Dm.ldar_program.LDARProgram(
        timeobj, copy.deepcopy(gas_field), {'ogi': ogi}, rep
    )
    tech_dict = {
        'plane': plane_survey,
        'ogi': ogi_no_survey
    }
    tiered_survey = Dm.ldar_program.LDARProgram(
        timeobj, gas_field, tech_dict, rep
    )
    feast.field_simulation.field_simulation(
            time=timeobj, gas_field=gas_field,
            ldar_program_dict={'tiered': tiered_survey, 'ogi': ogi_survey},
            dir_out='ResultsTemp', display_status=False
        )

    with open('ResultsTemp/realization0.p', 'rb') as f:
        res = pickle.load(f)
    if res.ldar_program_dict['tiered'].emissions_timeseries[-1] >= \
            res.ldar_program_dict['ogi'].emissions_timeseries[-1]:
        raise ValueError("field_simulation is not returning emission reductions as expected")
    if res.ldar_program_dict['ogi'].emissions_timeseries[-1] >= res.ldar_program_dict['Null'].emissions_timeseries[-1]:
        raise ValueError("field_simulation is not returning emission reductions as expected")

    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    try:
        os.rmdir('ResultsTemp')
    except PermissionError:
        # If there is an automated syncing process, a short pause may be necessary before removing "ResultsTemp"
        ti.sleep(5)
        os.rmdir('ResultsTemp')


test_repair()
test_comp_detect()
test_check_time()
test_site_detect()
test_ldar_program()
test_field_simulation()


print("Successfully completed LDAR tests.")
