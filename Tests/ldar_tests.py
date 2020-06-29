import numpy as np
import feast
from feast.GeneralClassesFunctions import simulation_classes as sc
import feast.GeneralClassesFunctions.emission_class_functions as ecf
import feast.GeneralClassesFunctions.results_analysis_functions as raf
from feast import DetectionModules as Dm
import os
import pickle
import time as ti

with open('test_helper.py', 'r') as f:
    exec(f.read())


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


def test_comp_detect_check_time():
    gas_field = basic_gas_field()
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


def test_ldar_program():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    find_cost = np.zeros(time.n_timesteps)
    rep = Dm.repair.Repair(repair_delay=0)
    ogi = Dm.comp_detect.CompDetect(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 800, 'end': 1700},
        labor=100,
        dispatch_object=rep
    )
    ogi_survey = Dm.ldar_program.LDARProgram(
        time, gas_field, [ogi], rep
    )

    if len(ogi_survey.find_cost) != 11:
        raise ValueError("find_cost not set to the correct length")

    if np.sum(ogi_survey.emissions.flux) != 100:
        raise ValueError("Unexpected emission rate in LDAR program initialization")

    ogi_survey.action(time, gas_field)

    if np.sum(ogi_survey.emissions.flux) != 85:
        raise ValueError("Unexpected emission rate after LDAR program action")

# feast.field_simulation.field_simulation()

# test_repair()
# test_comp_detect()
# test_comp_detect_check_time()
# test_ldar_program()

feast.field_simulation.field_simulation()

print("Successfully completed LDAR tests.")
