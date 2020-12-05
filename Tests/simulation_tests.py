import numpy as np
import copy
import feast
import os
import pickle
from feast import DetectionModules as Dm
from feast import EmissionSimModules as Esm
from feast.ResultsProcessing import results_analysis_functions as raf
from Tests.test_helper import basic_gas_field
import pandas as pd


def test_results_analysis():
    gf = basic_gas_field()
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=2)
    rep = Dm.repair.Repair(repair_delay=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    for ind in range(3):
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
        ogi_survey = Dm.ldar_program.LDARProgram(
            copy.deepcopy(gf), {'ogi': ogi}
        )
        sc = Esm.simulation_classes.Scenario(time=timeobj, gas_field=gf, ldar_program_dict={'ogi': ogi_survey})
        # todo: test 'json' save method
        # sc2 = copy.deepcopy(sc)
        sc.run(display_status=False, dir_out='ResultsTemp', save_method='pickle')
    null_npv, emissions, techs = raf.results_analysis('ResultsTemp', 0.08, 2e-4)
    if len(null_npv.keys()) != 4:
        raise ValueError("results analysis function returning the wrong number of keys")
    if null_npv['Finding'].shape != (1, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or LDAR programs in "
                         + "null npv")
    if emissions.shape != (2, 2, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or "
                         + "LDAR programs in emissions")
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')
    # sc2.run(display_status=False, dir_out='ResultsTemp', save_method='json')
    # res = pd.read_json(os.path.join('ResultsTemp', 'realization0.json'))


def test_npv_calculator():
    file_out = 'temp_emissions.p'
    rep_file_out = 'temp_rep_costs.p'
    em_array = np.ones(100)
    emissions = feast.input_data_classes.LeakData()
    leak_data = {'All': em_array}
    well_counts = {'All': 100}  # Number of wells in the study
    comp_counts = {'All': 100}  # Assumed components per well
    emissions.define_data(leak_data=leak_data, well_counts=well_counts, comp_counts=comp_counts)
    with open(file_out, 'wb') as f:
        pickle.dump(emissions, f)
    repair_out = feast.input_data_classes.RepairData()
    repair_out.define_data(repair_costs=np.ones(5) * 2)
    with open(rep_file_out, 'wb') as f:
        pickle.dump(repair_out, f)

    comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
        name='Fugitive emitters',
        emission_data_path=file_out,
        emission_per_comp=0.1,
        repair_cost_path=rep_file_out,
        emission_production_rate=0
    )
    n_sites = 100
    site_dict = {}
    basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
        name='basic pad',
        comp_dict={
            'Fugitive ': {'number': 100, 'parameters': comp_fug},
        }
    )
    site_dict['basic pad'] = {'number': n_sites, 'parameters': basicpad}
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=2)
    n_em = 1000
    initial_leaks = feast.EmissionSimModules.emission_class_functions.Emission(
        flux=np.ones(n_em), site_index=np.random.randint(0, n_sites, n_em),
        comp_index=np.random.randint(0, 100, n_em), end_time=np.infty, repair_cost=np.ones(n_em) * 2
    )
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj,
    )
    gas_field.emissions.emissions = gas_field.emissions.emissions[gas_field.emissions.emissions.start_time > 0]
    gas_field.emissions.emissions.index = list(np.linspace(n_em, n_em + len(gas_field.emissions.emissions.flux) - 1,
                                               len(gas_field.emissions.emissions.flux), dtype=int))
    gas_field.emissions.extend(initial_leaks)
    rep = Dm.repair.Repair(repair_delay=0)
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    ogi = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=50,
        survey_speed=10000,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=rep,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        site_queue=[]
    )
    ogi_survey = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), {'ogi': ogi}
    )
    sc = Esm.simulation_classes.Scenario(time=timeobj, gas_field=gas_field, ldar_program_dict={'ogi': ogi_survey})
    sc.run(dir_out='ResultsTemp', display_status=False, save_method='pickle')
    npv = feast.ResultsProcessing.results_analysis_functions.npv_calculator('ResultsTemp/realization0.p', 0.08, 2e-4)
    if npv['Repair'] != 2000:
        raise ValueError("npv_calculator not returning the expected repair cost")
    if np.abs(npv['Gas'] - 34560) > 100:  # 34560 = gas_value * 1000 g/s * 2 days * 3600 * 24 sec/day--no discount rate.
        raise ValueError("npv_calculator not returning the expected value of gas saved")
    os.remove(file_out)
    os.remove(rep_file_out)
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


def test_ResultsAggregate():
    res = Esm.result_classes.ResultAggregate(units='g/s')
    if res.units != 'g/s':
        raise ValueError("ResultsAggregate is not initializing self.units correctly")
    if res.time_value:
        raise ValueError("ResultsAggregate is not initializing self.time or self.value corerctly")
    res.append_entry([35, 6.2])
    if res.time_value[-1][0] != 35:
        raise ValueError("ResultsAggregate is not appending times correctly")
    if res.time_value[-1][1] != 6.2:
        raise ValueError("ResultsAggregate is not appending values correctly")
    v = res.get_vals(35)
    if v[0] != 6.2:
        raise ValueError("ResultsAggregate.get_value is not returning the correct entry")
    v = res.get_vals(36)
    if len(v) != 0:
        raise ValueError("ResultsAggregate.get_value is not returning an empty array if the requested time has no "
                         "value associated with it")


def test_ResultsDiscrete():
    res = Esm.result_classes.ResultDiscrete(units='g/s')
    if res.units != 'g/s':
        raise ValueError("ResultsDiscrete is not initializing self.units correctly")
    if res.time_value:
        raise ValueError("ResultsDiscrete is not initializing self.time or self.value corerctly")
    res.append_entry([35, 6.2])
    res.append_entry([37, 3.4])
    if res.time_value[0][0] != 35:
        raise ValueError("ResultsDiscrete is not appending times correctly")
    if res.time_value[-1][1] != 3.4:
        raise ValueError("ResultsDiscrete is not appending values correctly")
    if res.get_cumulative_vals(0, 38)[1][0] != 6.2 or res.get_cumulative_vals(0, 38)[1][1] != 9.6:
        raise ValueError("ResultsDiscrete is not computing cumulative vals correctly")
    if res.get_sum_val(0, 38) != 9.6:
        raise ValueError("ResultsDiscrete.get_sum_val is not computing sums correctly")


def test_ResultsContinuous():
    res = Esm.result_classes.ResultContinuous(units='g/s')
    if res.units != 'g/s':
        raise ValueError("ResultsDiscrete is not initializing self.units correctly")
    res.append_entry([0, 6.2])
    res.append_entry([37, 3.4])
    if res.get_time_integrated(0, 30) != 30 * 6.2:
        raise ValueError('ResultsContinuous.get_time_integrated is not integrating correctly')
    if res.get_time_integrated(0, 40) != 6.2 * 37 + (40 - 37) * 3.4:
        raise ValueError('ResultsContinuous.get_time_integrated is not integrating corretly')


test_results_analysis()
test_npv_calculator()
test_ResultsAggregate()
test_ResultsDiscrete()
test_ResultsContinuous()

print("Successfully completed simulation tests.")
