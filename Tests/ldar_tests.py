import numpy as np
import copy
import feast
import os
import pickle
import pandas as pd
from feast.EmissionSimModules import simulation_classes as sc
import feast.EmissionSimModules.emission_class_functions as ecf
from feast import DetectionModules as Dm
from Tests.test_helper import basic_gas_field
from Tests.test_helper import ex_prob_detect_arrays


def test_repair():
    n_em = 10
    flux, reparable, repair_cost = np.zeros(n_em), np.zeros(n_em, dtype=bool), np.ones(n_em)
    site_index, comp_index = np.zeros(n_em), np.zeros(n_em)
    flux[5:] = 1
    reparable[3: 8] = True
    endtime = np.linspace(0, 9, n_em)
    emission = ecf.Emission(flux, reparable, site_index, comp_index, repair_cost=repair_cost,
                            end_time=endtime)
    time = sc.Time()
    time.time_index = 3
    time.current_time = 4.5
    detected = np.array([3, 4, 5, 6, 7, 8, 9], dtype=int)
    repair_proc = Dm.repair.Repair(repair_delay=1)
    repair_proc.action(emit_inds=detected)
    repair_proc.repair(time, emission)
    expected = np.array([0., 1., 2., 3., 4., 5., 5.5, 5.5, 8., 9.])
    for ind in range(len(emission.emissions.end_time)):
        if emission.emissions.end_time[ind] != expected[ind]:
            raise ValueError("DetectionModules.repair.Repair is not adjusting "
                             "emission endtimes correctly at index {:0.0f}".format(ind))
    if repair_proc.repair_count.get_sum_val(0, 10) != 2:
        raise ValueError("Repair.repair_count is not updated correctly.")
    if repair_proc.repair_cost.get_sum_val(0, 10) != 2:
        raise ValueError("Repair.repai_cost is not updated correctly")


def test_check_time():
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_speed=150,
        labor=100,
        site_queue=[0, 1, 2, 3],
        detection_probability_points=[1, 2, 3],
        detection_probabilities=[0, 0.5, 1],
        detection_variables={'flux': 'mean'},
        survey_interval=50,
        dispatch_object=rep,
        ophrs={'begin': 8, 'end': 17}
    )
    if not tech.check_time(time):
        raise ValueError("Check time returning False when it should be true at time 0")
    time = sc.Time(delta_t=0.1, end_time=10, current_time=0)
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_speed=150,
        labor=100,
        site_queue=[0, 1, 2, 3],
        detection_probability_points=[1, 2, 3],
        detection_probabilities=[0, 0.5, 1],
        detection_variables={'flux': 'mean'},
        survey_interval=50,
        dispatch_object=rep,
        ophrs={'begin': 8, 'end': 17}
    )
    if tech.check_time(time):
        raise ValueError("check_time returning True when it should return False at time 0")


def test_comp_survey():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    tech = Dm.comp_survey.CompSurvey(
        time,
        site_queue=list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)),
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    emissions = gas_field.emissions.get_current_emissions(time)
    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
    tech.detect(time, gas_field, emissions)
    if tech.deployment_count.get_sum_val(0, 1) != 14:
        raise ValueError("deployment count is not updated correctly in comp_survey")
    if np.abs(tech.deployment_cost.get_sum_val(0, 1) - 9 * 100) > 1e-5:
        raise ValueError("CompSurvey is not calculating deployment costs correctly")
    if np.max(emissions.site_index[emissions.index.isin(rep.to_repair)]) != 13:
        raise ValueError("tech.detect repairing emissions at incorrect sites")
    expected_detected_ids = [54, 55, 75, 66, 87, 62, 58, 84, 5, 25, 43, 71, 13, 79, 97]
    for ind in range(len(rep.to_repair)):
        if expected_detected_ids[ind] != rep.to_repair[ind]:
            raise ValueError("tech.detect not detecting the correct emissions")
    if len(expected_detected_ids) != len(rep.to_repair):
        raise ValueError("tech.detect not detecting the correct emissoins")
    rep.repair(time, gas_field.emissions)
    expected_repaired = gas_field.emissions.emissions.index.isin(expected_detected_ids)
    if np.max(gas_field.emissions.emissions.end_time[expected_repaired]) != 0:
        raise ValueError("rep.repair not adjusting the end times correctly")


def test_comp_survey_emitters_surveyed():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    rep = Dm.repair.Repair(repair_delay=0)
    wind_dirs_mins = np.zeros(gas_field.n_sites)
    wind_dirs_maxs = np.ones(gas_field.n_sites) * 90
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        op_envelope={
            'wind speed': {'class': 1, 'min': 1, 'max': 10},
            'wind direction': {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
        },
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        site_queue=[]
    )
    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
    emissions = gas_field.emissions.get_current_emissions(time)
    emitter_inds = tech.emitters_surveyed(time, gas_field, emissions)
    if emitter_inds:
        # emitter_inds is expected to be []
        raise ValueError("CompSurvey.emitters_surveyed is not returning expected emitter indexes")
    wind_dirs_maxs[11] = 200
    emitter_inds = tech.emitters_surveyed(time, gas_field, emissions)
    if emitter_inds != [71]:
        # The wind direction op envelope was updated to pass at site 11 only. Site 11 has one emission at index 71.
        raise ValueError("CompSurvey.emitters_surveyed is not returning expected indexes")


def test_site_survey():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.474)) / (1.36 * np.sqrt(2))) for f
                                  in points])
    rep = Dm.repair.Repair(repair_delay=0)
    comp_temp = Dm.comp_survey.CompSurvey(
        time,
        site_queue=[],
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    tech = Dm.site_survey.SiteSurvey(
        time,
        survey_interval=50,
        sites_per_day=100,
        ophrs={'begin': 8, 'end': 17},
        site_cost=100,
        dispatch_object=comp_temp,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    emissions = gas_field.emissions.get_current_emissions(time)
    np.random.seed(0)
    # test detect_prob_curve
    detect = tech.detect_prob_curve(time, gas_field, [0, 1, 2], emissions)
    if detect != np.array([0]):
        raise ValueError("SiteSurvey.detect_prob_curve not returning expected sites.")
    # test sites_surveyed with empty queue
    sites_surveyed = tech.sites_surveyed(gas_field, time)
    if sites_surveyed:
        raise ValueError("sites_surveyed returning sites when it should not")
    # test detect
    np.random.seed(0)
    tech.site_queue = [0, 1, 2]
    tech.detect(time, gas_field, emissions)
    if tech.deployment_count.get_sum_val(0, 1) != 3:
        raise ValueError("SiteSurvey is not counting deployments correctly")
    if tech.detection_count.get_sum_val(0, 1) != 1:
        raise ValueError("SiteSurvey is not counting detections correctly.")
    if tech.deployment_cost.get_sum_val(0, 1) != 300:
        raise ValueError("SiteSurvey is not calculating survey costs correctly")
    if tech.dispatch_object.site_queue != [0]:
        raise ValueError("site_detect.detect not updating dispatch object sites to survey correctly")
    # test action and sites_surveyed with
    tech.action(list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)))
    if tech.site_queue != list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)):
        raise ValueError("action is not updating site_queue as expected")
    # test sites_surveyed with full queue
    sites_surveyed = tech.sites_surveyed(gas_field, time)
    if (sites_surveyed != np.linspace(0, 99, 100, dtype=int)).any():
        raise ValueError("sites_surveyed not identifying the correct sites")


def test_sitedetect_sites_surveyed():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    wind_dirs_mins = np.zeros(gas_field.n_sites)
    wind_dirs_maxs = np.ones(gas_field.n_sites) * 90
    wind_dirs_maxs[50] = 270
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.474)) / (1.36 * np.sqrt(2))) for f
                                  in points])
    comp_survey = Dm.comp_survey.CompSurvey(
        time,
        dispatch_object=None,
        survey_interval=None,
        survey_speed=100,
        labor=100,
        site_queue=[],
        detection_probability_points=[1, 2],
        detection_probabilities=[0, 1],
        ophrs={'begin': 8, 'end': 17}
    )
    tech = Dm.site_survey.SiteSurvey(
        time,
        survey_interval=50,
        sites_per_day=100,
        ophrs={'begin': 8, 'end': 17},
        site_cost=100,
        dispatch_object=comp_survey,
        op_envelope={
            'wind speed': {'class': 1, 'min': 1, 'max': 10},
            'wind direction': {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
        },
        detection_probability_points=points,
        detection_probabilities=probs
    )
    tech.site_queue = list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int))
    np.random.seed(0)
    sites_surveyed = tech.sites_surveyed(gas_field, time)
    if sites_surveyed != [50]:
        raise ValueError("sites_surveyed is not returning sites correctly")


def test_ldar_program():
    gas_field = basic_gas_field()
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    probs[0] = 0
    ogi = Dm.comp_survey.CompSurvey(
        time,
        site_queue=list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)),
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    ogi_no_survey = Dm.comp_survey.CompSurvey(
        time,
        site_queue=[],
        survey_interval=None,
        survey_speed=100,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    points = [0.5, 0.6]
    probs = [0, 1]
    probs[0] = 0
    plane_survey = Dm.site_survey.SiteSurvey(
        time,
        survey_interval=50,
        sites_per_day=200,
        site_cost=100,
        dispatch_object=ogi_no_survey,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    # test __init__
    ogi_survey = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), {'ogi': ogi}
    )
    em = ogi_survey.emissions.get_current_emissions(time)
    if np.sum(em.flux) != 100:
        raise ValueError("Unexpected emission rate in LDAR program initialization")

    # test end_emissions
    ogi_survey.emissions.emissions.loc[0, 'end_time'] = 0
    # test action
    ogi_survey.action(time, gas_field)
    em = ogi_survey.emissions.get_current_emissions(time)
    if np.sum(em.flux) != 84:
        raise ValueError("Unexpected emission rate after LDAR program action")
    # test combined program
    tech_dict = {
        'plane': plane_survey,
        'ogi': ogi_no_survey
    }
    tiered_survey = Dm.ldar_program.LDARProgram(
        gas_field, tech_dict
    )
    # test action
    np.random.seed(0)
    tiered_survey.action(time, gas_field)
    em = tiered_survey.emissions.get_current_emissions(time)
    if np.sum(em.flux) != 86:
        raise ValueError("Unexpected emission rate after LDAR program action with tiered survey")
    time.current_time = 1
    tiered_survey.calc_rep_costs(time)
    if tiered_survey.repair_cost.time_value[0][0] != 0 or len(tiered_survey.repair_cost.time_value) != 14 or \
            tiered_survey.repair_cost.time_value[-1][-1] != 2:
        raise ValueError("tiered_survey.calc_rep_costs is not calculating the correct repair costs")


def test_scenario_run():
    gas_field = basic_gas_field()
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=2)
    rep = Dm.repair.Repair(repair_delay=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    ogi = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        site_queue=[]
    )
    ogi_no_survey = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=None,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        site_queue=[]
    )
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.474)) / (1.36 * np.sqrt(2))) for f
                                  in points])
    plane_survey = Dm.site_survey.SiteSurvey(
        timeobj,
        survey_interval=50,
        sites_per_day=200,
        site_cost=100,
        dispatch_object=ogi_no_survey,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    ogi_survey = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), {'ogi': ogi}
    )
    tech_dict = {
        'plane': plane_survey,
        'ogi': ogi_no_survey
    }
    tiered_survey = Dm.ldar_program.LDARProgram(
        gas_field, tech_dict
    )
    scenario = sc.Scenario(time=timeobj, gas_field=gas_field, ldar_program_dict={'tiered': tiered_survey,
                                                                                 'ogi': ogi_survey})
    scenario.run(dir_out='ResultsTemp', display_status=False, save_method='pickle')
    with open('ResultsTemp/realization0.p', 'rb') as f:
        res = pickle.load(f)
    if res.ldar_program_dict['tiered'].emissions_timeseries[-1] >= \
            res.ldar_program_dict['ogi'].emissions_timeseries[-1]:
        raise ValueError("Scenario.run is not returning emission reductions as expected")
    if res.ldar_program_dict['ogi'].emissions_timeseries[-1] >= res.ldar_program_dict['Null'].emissions_timeseries[-1]:
        raise ValueError("Scenario.run is not returning emission reductions as expected")

    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


def test_check_op_envelope():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    rep = Dm.repair.Repair(repair_delay=0)
    wind_dirs_mins = np.zeros(gas_field.n_sites)
    wind_dirs_maxs = np.ones(gas_field.n_sites) * 90
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        op_envelope={
            'wind speed': {'class': 1, 'min': 1, 'max': 10},
            'wind direction': {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
        },
        site_queue=[],
        detection_probability_points=points,
        detection_probabilities=probs
    )
    op_env = tech.check_op_envelope(gas_field, time, 0)
    if op_env != 'site fail':
        raise ValueError("check_op_envelope is not returning 'site fail' as expected")
    if len(tech.op_env_site_fails.get_vals(0, 1)) != 1 or tech.op_env_site_fails.get_vals(0, 1)[0] != 1:
        raise ValueError("DetectionMethod.op_env_site_fails is not updated correctly")
    wind_dirs_mins = np.zeros([gas_field.n_sites, 2])
    wind_dirs_mins[:, 1] = 145
    wind_dirs_maxs = np.ones([gas_field.n_sites, 2]) * 90
    wind_dirs_maxs[:, 1] += 145
    tech.op_envelope['wind direction'] = {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
    op_env = tech.check_op_envelope(gas_field, time, 0)
    if op_env != 'site pass':
        raise ValueError("check_op_envelope is not passing as expected")
    tech.op_envelope['wind speed']['max'] = [2]
    op_env = tech.check_op_envelope(gas_field, time, 0)
    if op_env != 'field fail':
        raise ValueError("check_op_envelope is not retruning 'field fail' as expected")
    if tech.op_env_field_fails.get_sum_val(0, 1) != 1:
        raise ValueError("DetectionMethod.op_env_field_fails is not updated correctly")


def test_get_current_conditions():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    rep = Dm.repair.Repair(repair_delay=0)
    prob_points, detect_probs = ex_prob_detect_arrays()
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean', 'wind speed': 'mean'},
        detection_probability_points=prob_points,
        detection_probabilities=detect_probs,
        site_queue=[]
    )
    # pd.DataFrame initializes a new DataFrame to avoid chained indexing warnings
    emissions = pd.DataFrame(gas_field.emissions.get_current_emissions(time))
    n_em = 100
    emissions.loc[:, 'flux'] = np.linspace(0.1, 10, n_em)
    em_indexes = np.linspace(0, n_em - 1, n_em, dtype=int)
    ret = tech.get_current_conditions(time, gas_field, emissions, em_indexes)
    if np.any(ret[:, 0] != emissions.flux):
        raise ValueError("get_current_conditions not returning the correct values")
    if np.any(ret[:, 1] != np.mean(gas_field.met['wind speed'][8:17])):
        raise ValueError("get_current_conditions not returning the correct values")


def test_empirical_interpolator():
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    rep = Dm.repair.Repair(repair_delay=0)
    prob_points, detect_probs = ex_prob_detect_arrays()
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean', 'wind speed': 'mean'},
        detection_probability_points=prob_points,
        detection_probabilities=detect_probs,
        site_queue=[]
    )
    probs = tech.empirical_interpolator(tech.detection_probability_points, tech.detection_probabilities,
                                        np.array([[0.01, 1], [0.05, 1]]))
    if np.abs(probs[0] - 0.24509709) > 1e-5 or np.abs(probs[1] - 0.25784611) > 1e-5:
        raise ValueError("empirical_interpolator is not returning the correct probabilities")
    probs = tech.empirical_interpolator(tech.detection_probability_points, tech.detection_probabilities,
                                        np.array([0.03, 1]))
    if not min(tech.detection_probabilities[:2]) <= probs[0] <= max(tech.detection_probabilities[:2]):
        raise ValueError("empirical_interpolator is not interpolating correctly")
    probs = tech.empirical_interpolator(tech.detection_probability_points, tech.detection_probabilities,
                                        np.array([0.01, 1.5]))
    if not tech.detection_probabilities[0] >= probs[0] >= tech.detection_probabilities[6]:
        raise ValueError("empirical_interpolator is not interpolating correctly")


def test_choose_sites():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    rep = Dm.repair.Repair(repair_delay=0)
    # wind_dirs_mins = np.zeros(gas_field.n_sites)
    # wind_dirs_maxs = np.ones(gas_field.n_sites) * 90
    tech = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        op_envelope={
            # 'wind speed': {'class': 1, 'min': 1, 'max': 10},
            # 'wind direction': {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
        },
        site_queue=[],
        detection_probability_points=[1, 2],
        detection_probabilities=[0, 1]
    )
    tech.site_queue = list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int))
    siteinds = tech.choose_sites(gas_field, time, 10)
    if not (siteinds == np.linspace(0, 9, 10, dtype=int)).all():
        raise ValueError("choose_sites() is not selecting the correcting sites")

    tech.site_queue = []
    siteinds = tech.choose_sites(gas_field, time, 10)
    if siteinds:
        raise ValueError("choose_sites() fails for empty site_queue queue")


def test_site_monitor():
    gas_field = basic_gas_field()
    gas_field.met_data_path = 'TMY-DataExample.csv'
    time = sc.Time(delta_t=1, end_time=10, current_time=0)
    gas_field.met_data_maker()
    rep = Dm.repair.Repair(repair_delay=0)
    ogi = Dm.comp_survey.CompSurvey(
        time,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        op_envelope={
            # 'wind speed': {'class': 1, 'min': 1, 'max': 10},
            # 'wind direction': {'class': 2, 'min': wind_dirs_mins, 'max': wind_dirs_maxs}
        },
        site_queue=[],
        detection_probability_points=[1, 2],
        detection_probabilities=[0, 1]
    )
    cm = Dm.site_monitor.SiteMonitor(
        time,
        time_to_detect_points=np.array([0.99, 1.0, 1.01]),
        time_to_detect_days=[np.infty, 1, 0],
        detection_variables={'flux': 'mean'},
        dispatch_object=ogi,
        capital=1000,
        ophrs={'begin': 0, 'end': 24},
        site_queue=list(np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int))
    )
    site_inds = list(range(0, 10))
    emissions = copy.copy(gas_field.emissions.get_current_emissions(time))
    detect = cm.detect_prob_curve(time, gas_field, site_inds, emissions)
    must_detect = [0, 9]
    must_not_detect = [2, 7, 8]
    for md in must_detect:
        if md not in detect:
            raise ValueError("site_monitor.detect_prob_curve not flagging the correct sites")
    for mnd in must_not_detect:
        if mnd in detect:
            raise ValueError("site_monitor.detect_prob_curve flagging sites that it should not")
    ttd_list = []
    ttd = 0
    time.delta_t = 0.1
    for ind in range(1000):
        ttd += time.delta_t
        detect = cm.detect_prob_curve(time, gas_field, site_inds, emissions)
        if 1 in detect:
            ttd_list.append(ttd)
            ttd = 0
    if np.abs(np.mean(ttd_list) - 1) > 0.5:
        raise ValueError("Mean time to detection deviates from expected value by >5 sigma in site_monitor test")
    if cm.deployment_cost.get_sum_val(0, 1) != 1000:
        raise ValueError("SiteMonitor deployment cost is not calculated correctly.")
    cm.detect(time, gas_field, emissions)


test_repair()
test_comp_survey()
test_check_time()
test_site_survey()
test_ldar_program()
test_scenario_run()
test_check_op_envelope()
test_sitedetect_sites_surveyed()
test_comp_survey_emitters_surveyed()
test_get_current_conditions()
test_empirical_interpolator()
test_choose_sites()
test_site_monitor()


print("Successfully completed LDAR tests.")
