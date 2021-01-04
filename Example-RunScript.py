"""
This script provides an example of a custom FEAST simulation. The script defines a set of four custom components,
and combines those components to create a variety of sites. The sites are then combined into a gas field that is used in
the simulation. A component level survey, a tiered site level survey and a tiered site monitor program are simulated.
"""
import numpy as np
import copy
import feast.EmissionSimModules.infrastructure_classes
import feast.EmissionSimModules.simulation_classes as sc
from feast.EmissionSimModules.infrastructure_classes import Component
import feast.DetectionModules as Dm
import pickle
import time

# The random seed below can be un-commented to generate reproducible realizations.
# np.random.seed(0)
n_montecarlo = 5
a = time.time()


def define_emitters():
    """
    Defines all emitters to be used in the simulation using the Component class
    :return comp_fug: source of fugitive emissions
    :return misc_vent: source of miscelaneous vents
    :return plunger: source of unloading emissions due to wells with plungers
    :return noplunger: source of unloading emission due to wells without plungers
    """
    # Generates reparable fugitive emissions
    comp_fug = Component(
        name='Fugitive emitters',
        emission_data_path='ExampleData/DataObjectInstances/production_emissions.p',
        emission_per_comp=0.00231,  # Fraction of components expected to be emitting at the beginning of the simulation.
        emission_production_rate=5.4 / 650 / 365,  # number of new emissions per component per day
        repair_cost_path='ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        base_reparable=True
    )

    # Generates miscelaneous emissions that cannot be repaired
    # (many tank emissions and pneumatic controller emissions can be modeled in this way).
    misc_vent = Component(
        name='misc vents',
        emission_data_path='ExampleData/DataObjectInstances/production_emissions.p',
        emission_per_comp=0.00231,  # Fraction of components expected to be emitting at the beginning of the simulation.
        emission_production_rate=5.4 / 650 / 365,  # number of new emissions per component per day
        repair_cost_path='ExampleData/DataObjectInstances/fernandez_leak_repair_costs_2006.p',
        base_reparable=False,
    )

    # Simulates unloading emissions from a well with a plunger lift
    plunger = feast.EmissionSimModules.infrastructure_classes.Component(
        name='plunger',
        emission_production_rate=0,  # Expected rate of production for fugitive emissions
        episodic_emission_sizes=[50.36],  # gram-per-sec
        episodic_emission_duration=2059 / 3600 / 24,  # Average length of emission (days)
        episodic_emission_per_day=0.0165  # Number of events per day per plunger well
    )

    # Simulates unloading emissions from a well with out a plunger lift
    noplunger = feast.EmissionSimModules.infrastructure_classes.Component(
        name='no plunger',
        emission_production_rate=0,
        episodic_emission_sizes=[76.98],  # gram-per-sec
        episodic_emission_duration=4973.5 / 3600 / 24,  # Average length of emission (days)
        episodic_emission_per_day=0.0002111  # Number of events per day for unloading wells without a plunger
    )
    return comp_fug, misc_vent, plunger, noplunger


def define_sites(comp_fug, misc_vent, plunger, noplunger):
    """
    Defines all sites to be used in the simulation. Sites consist of a collection of Component objects
    :param comp_fug: fugitive emitters
    :param misc_vent: miscelaneous vents
    :param plunger: unloading emissions from wells with plungers
    :param noplunger: unloading emissions from wells without plungers
    :return site_dict: dict of all sites to be simulated
    """
    # Define the number of wells at each site
    # The well per site distribution was built using data from the Colorado Oil and Gas Conservation Commission
    counts = [78, 8, 4, 3, 2, 1, 1, 1, 1, 1]
    n_wells = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    # ------Each iteration of this for loop generates one realization of the simulation
    site_dict = {}

    # Assign components to sites
    for n_wells_in_site in n_wells:
        basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
            name='basic pad',
            comp_dict={
                # ------ The number of components is proportional to the number of wells, ind2
                'Fugitive ' + str(n_wells_in_site): {'number': 300 * n_wells_in_site,
                                                     'parameters': copy.copy(comp_fug)},
                'misc vent ' + str(n_wells_in_site): {'number': 350 * n_wells_in_site,
                                                      'parameters': copy.copy(misc_vent)}
            },
        )
        site_dict['basic pad ' + str(n_wells_in_site)] = {'number': counts[n_wells_in_site - 1],
                                                          'parameters': basicpad}
    plungersite = feast.EmissionSimModules.infrastructure_classes.Site(
        name='plunger well',
        comp_dict={
            'Fugitive_plunger': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'plunger': {'number': 2, 'parameters': copy.copy(plunger)},
            'misc vent plunger': {'number': 600, 'parameters': copy.copy(misc_vent)}
        },
    )
    unloadnp = feast.EmissionSimModules.infrastructure_classes.Site(
        name='unload no plunger',
        comp_dict={
            'Fugitive_np': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'no plunger': {'number': 2, 'parameters': copy.copy(noplunger)},
            'misc vent np': {'number': 600, 'parameters': copy.copy(misc_vent)}
        },
    )
    site_dict['plunger'] = {'number': 6, 'parameters': plungersite}
    site_dict['unload no plunger'] = {'number': 1, 'parameters': unloadnp}
    return site_dict


def define_time_settings():
    """
    Sets the time resolution and duration of the simulation
    :return: a Time object containing simulation settings
    """
    return feast.EmissionSimModules.simulation_classes.Time(delta_t=1, end_time=365)


def define_gas_field(timeobj, site_dict):
    """
    Creates a GasField object to be used in the simulation
    :param timeobj: A time object containing time settings for the simulation
    :param site_dict: A dict of sites to be included in the gas field
    :return gas_field: A GasField object to be used in the simulation
    """
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    gas_field.met_data_path = 'ExampleData/TMY-DataExample.csv'
    gas_field.met_data_maker()
    return gas_field


def define_detection_methods(timeobj):
    """
    Define detection methods to be used in LDAR programs
    :param timeobj: A time object for simulation settings
    :return ogi: A component survey method representing OGI with periodic surveys
    :return ogi_no_survey: A component survey method representing OGI deployed by a site-level detection method
    :return plane: A site survey method representing a plane based detection program with periodic surveys
    :return cont_monitor: A site monitor method representing continuous monitors deployed at a site
    :return rep0: A repair method with 0 delay between detection and repair
    :return rep7: A repair method with a delay of 7 days between detection and repair
    """
    points = np.logspace(-3, 1, 100)  # emission rates
    rep0 = Dm.repair.Repair(repair_delay=0)
    rep7 = Dm.repair.Repair(repair_delay=7)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    probs[0] = 0
    ogi = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=180,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,  # g/s
        detection_probabilities=probs,
        dispatch_object=rep0,
        site_queue=[],
    )
    ogi_no_survey = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=None,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        dispatch_object=copy.copy(rep0),
        site_queue=[],
    )
    points = np.logspace(-3, 1, 100)
    # 0.474
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(1.5)) / (1.36 * np.sqrt(2))) for f
                                  in points])
    probs[0] = 0
    plane_survey = Dm.site_survey.SiteSurvey(
        timeobj,
        survey_interval=180,
        sites_per_day=200,
        site_cost=100,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs,
        dispatch_object=ogi_no_survey,
        site_queue=[],
        ophrs={'begin': 8, 'end': 17}
    )
    cont_monitor = Dm.site_monitor.SiteMonitor(
        timeobj,
        time_to_detect_points=[[0.5, 1], [1.0, 1], [1.1, 1], [0.5, 5], [1.0, 5], [1.1, 5],
                               [0.5, 5.1], [1.0, 5.1], [1.1, 5.1]],
        time_to_detect_days=[np.infty, 1, 0, np.infty, 5, 0, np.infty, np.infty, np.infty],
        detection_variables={'flux': 'mean', 'wind speed': 'mean'},
        site_queue=list(range(gas_field.n_sites)),
        dispatch_object=copy.deepcopy(rep0),
        ophrs={'begin': 8, 'end': 17}
    )
    return ogi, ogi_no_survey, plane_survey, cont_monitor, rep0, rep7


def define_ldar_programs(gas_field, ogi, ogi_no_survey, plane_survey, cont_monitor, rep0, rep7):
    """
    Define LDAR programs using the detection and repair methods defined previously
    :param gas_field: Emission simulation settings
    :param ogi:  component survey method representing OGI with periodic surveys
    :param ogi_no_survey: A component survey method representing OGI deployed by a site-level detection method
    :param plane_survey: A site survey method representing a plane based detection program with periodic surveys
    :param cont_monitor: A site monitor method representing continuous monitors deployed at a site
    :param rep0: A repair method with 0 delay between detection and repair
    :param rep7: A repair method with a delay of 7 days between detection and repair
    :return ldar_dict: A dict of LDAR programs to be simulated
    """
    # Add dispatch methods and site specific conditions to detection methods
    # Good practice to use copies so that LDAR programs do not interfere with eachother in the simulation
    ogi.dispatch_object = copy.deepcopy(rep0)
    ogi_no_survey.dispatch_object = copy.deepcopy(rep0)
    plane_ogi = copy.deepcopy(ogi_no_survey)
    plane_survey.dispatch_object = plane_ogi
    cm_ogi = copy.deepcopy(ogi_no_survey)
    cont_monitor.dispatch_object = cm_ogi
    cont_monitor.site_queue = np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int)
    cont_monitor.op_envelope = {
        'wind direction': {'class': 2,
                           'min': [[45, 225]] * gas_field.n_sites,
                           'max': [[135, 315]] * gas_field.n_sites}
    }
    # Define LDAR programs
    ogi_survey = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), {'ogi': ogi},
    )
    # tiered survey
    tech_dict = {
        'plane': plane_survey,
        'ogi': plane_ogi
    }
    plane_ogi_survey = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), tech_dict,
    )

    # continuous monitor
    tech_dict = {
        'cm': cont_monitor,
        'ogi': cm_ogi
    }
    cm_ogi = Dm.ldar_program.LDARProgram(
        copy.deepcopy(gas_field), tech_dict,
    )

    # All programs
    ldar_dict = {
        'cm': cm_ogi,
        'ogi': ogi_survey,
        'plane': plane_ogi_survey
    }
    return ldar_dict


for ind in range(n_montecarlo):
    print('Iteration number: {:0.0f}'.format(ind))
    comp_fug, misc_vent, plunger, noplunger = define_emitters()
    site_dict = define_sites(comp_fug, misc_vent, plunger, noplunger)
    timeobj = define_time_settings()
    gas_field = define_gas_field(timeobj, site_dict)
    ogi, ogi_no_survey, plane_survey, cont_monitor, rep0, rep7 = define_detection_methods(timeobj)
    ldar_dict = define_ldar_programs(gas_field, ogi, ogi_no_survey, plane_survey, cont_monitor, rep0, rep7)
    scenario = sc.Scenario(time=timeobj, gas_field=gas_field, ldar_program_dict=ldar_dict)
    scenario.run(dir_out='ExampleRunScriptResults', display_status=True, save_method='pickle')

b = time.time()
print("run time {:0.2f} seconds".format(b - a))
