import feast
import numpy as np
import copy

comp_fug = feast.EmissionSimModules.infrastructure_classes.Component(
    name='Fugitive emitters',
    emission_data_path='ExampleData/DataObjectInstances/production_emissions.p',
    emission_per_comp=0.00231,  # Fraction of components expected to be emitting at the beginning of the simulation.
    # emission_production_rate is set to about 5 new emissions per well per year
    emission_production_rate=5.4 / 650 / 365,  # number of new emissions per component per day
    repair_cost_path='ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
    base_reparable=True
)

wells_per_site = 2
comps_per_well = 650
basic_site = feast.EmissionSimModules.infrastructure_classes.Site(
    name='basic site',
    comp_dict={
        # ------ The number of components is proportional to the number of wells, ind2
        'Fugitive ': {'number': comps_per_well * wells_per_site,
                     'parameters': copy.deepcopy(comp_fug)},
    },
)


time = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=365)

site_dict = {
    'basic site': {
        'number': 100, # The number of basic_site instances to simulate in the gas_field
        'parameters': basic_site
    }
    # If additional site types are needed, they should be added here as separate dict entries
}
gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
    sites=site_dict,
    time=time,
    met_data_path='ExampleData/TMY-DataExample.csv'
)

rep3 = feast.DetectionModules.repair.Repair(repair_delay=3)


def make_ogi(dispatch_obj):
    points = np.array([0.5, 1, 2, 3, 7]) * 0.01157 # g/s
    probs = np.array([0, 0.25, 0.5, 0.75, 1])
    ogi = feast.DetectionModules.comp_survey.CompSurvey(
        time,
        survey_interval=180, # days
        survey_speed=500, # comps/hr
        ophrs={'begin': 8, 'end': 17},
        labor=400, # $/hr
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,  # g/s
        detection_probabilities=probs,
        dispatch_object=copy.deepcopy(dispatch_obj),
        site_queue=[],
    )
    return ogi


# An OGI-like example that does not automatically survey components but can be dispatched
def make_ogi_no_survey(dispatch_obj):
    points = np.array([0.5, 1, 2, 3, 7]) * 0.01157 # g/s
    probs = np.array([0, 0.25, 0.5, 0.75, 1])
    ogi_no_survey = feast.DetectionModules.comp_survey.CompSurvey(
        time,
        survey_interval=None, # days
        survey_speed=500, # comps/hr
        ophrs={'begin': 8, 'end': 17},
        labor=400, # $/hr
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,  # g/s
        detection_probabilities=probs,
        dispatch_object=copy.deepcopy(dispatch_obj),
        site_queue=[],
    )
    return ogi_no_survey


def make_plane_survey(dispatch_obj):
    points = np.array([25, 50, 100, 200, 400]) * 0.01157 # g/s
    probs = np.array([0, 0.25, 0.5, 0.75, 1])
    plane_survey = feast.DetectionModules.site_survey.SiteSurvey(
        time,
        survey_interval=180, # days
        sites_per_day=200,
        site_cost=100, # $/site
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        dispatch_object=copy.deepcopy(dispatch_obj),
        site_queue=[],
        ophrs={'begin': 8, 'end': 17}
    )
    return plane_survey


def make_cont_monitor(dispatch_obj):
    points = [[0.5, 1],
              [1.0, 1],
              [1.1, 1],
              [0.5, 5],
              [1.0, 5],
              [1.1, 5],
              [0.5, 5.1],
              [1.0, 5.1],
              [1.1, 5.1]]
    time_to_detect_days = [np.infty, 1, 0, np.infty, 5, 0, np.infty, np.infty, np.infty]
    cont_monitor = feast.DetectionModules.site_monitor.SiteMonitor(
        time,
        time_to_detect_points=points,
        time_to_detect_days=time_to_detect_days,
        detection_variables={'flux': 'mean', 'wind speed': 'mean'},
        site_queue=list(range(gas_field.n_sites)),
        dispatch_object=copy.deepcopy(dispatch_obj),
        ophrs={'begin': 8, 'end': 17}
    )
    return cont_monitor


def make_iteration(ind):
    print("Currently evaluating iteration number {:0.0f}".format(ind))
    time = feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=365)
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=time,
        met_data_path='ExampleData/TMY-DataExample.csv'
    )
    ogi = make_ogi(rep3)
    ogi_no_survey = make_ogi_no_survey(rep3)
    plane_survey = make_plane_survey(ogi_no_survey)
    # LDAR programs
    ogi_survey = feast.DetectionModules.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), {'ogi': ogi},
    )
    # tiered survey
    tech_dict = {
        'plane': plane_survey,
        'ogi': plane_survey.dispatch_object
    }
    plane_ogi_survey = feast.DetectionModules.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), tech_dict,
    )
    ldar_dict = {
        'ogi': ogi_survey,
        'plane': plane_ogi_survey
    }
    scenario = feast.EmissionSimModules.simulation_classes.Scenario(time=time,
                                                                gas_field=gas_field,
                                                                ldar_program_dict=ldar_dict)
    scenario.run(dir_out='TutorialResults-MC', display_status=False, save_method='pickle')