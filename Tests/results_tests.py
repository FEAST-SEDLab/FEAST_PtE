import numpy as np
import copy
import feast
import os
import pickle
from feast import DetectionModules as Dm
from feast.ResultsProcessing import results_analysis_functions as raf
from Tests.test_helper import basic_gas_field


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
            timeobj, copy.deepcopy(gf), {'ogi': ogi}
        )
        feast.field_simulation.field_simulation(gf, time=timeobj,
                                                ldar_program_dict={'ogi': ogi_survey},
                                                dir_out='ResultsTemp', display_status=False)
    null_npv, emissions, costs, techs = raf.results_analysis('ResultsTemp')
    if len(null_npv.keys()) != 4:
        raise ValueError("results analysis function returning the wrong number of keys")
    if null_npv['Finding'].shape != (1, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or LDAR programs in "
                         + "null npv")
    if emissions.shape != (2, 3, 3):
        raise ValueError("results analysis function returning the wrong number of realizations or "
                         + "LDAR programs in emissions")
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


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
    initial_leaks = feast.EmissionSimModules.emission_class_functions.Emission(
        flux=np.ones(1000), site_index=np.random.randint(0, n_sites, 1000),
        comp_index=np.random.randint(0, 100, 1000), endtime=np.infty, repair_cost=np.ones(1000) * 2
    )
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj,
        initial_emissions=initial_leaks
    )
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
        timeobj, copy.deepcopy(gas_field), {'ogi': ogi}
    )
    econ_set = feast.EmissionSimModules.simulation_classes.FinanceSettings(gas_price=2e-4, discount_rate=0.08)
    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field, dir_out='ResultsTemp',
        display_status=False, ldar_program_dict={'ogi': ogi_survey},
        econ_set=econ_set
    )
    npv = feast.ResultsProcessing.results_analysis_functions.npv_calculator('ResultsTemp/realization0.p')
    if npv['Repair'] != 2000:
        raise ValueError("npv_calculator not returning the expected repair cost")
    if np.abs(npv['Gas'] - 34560) > 100:  # 34560 = gas_value * 1000 g/s * 2 days * 3600 * 24 sec/day--no discount rate.
        raise ValueError("npv_calculator not returning the expected value of gas saved")
    os.remove(file_out)
    os.remove(rep_file_out)
    for f in os.listdir('ResultsTemp'):
        os.remove(os.path.join('ResultsTemp', f))
    os.rmdir('ResultsTemp')


test_results_analysis()
test_npv_calculator()


print("Successfully completed results tests.")
