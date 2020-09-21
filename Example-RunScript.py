"""
This script provides an example of a custom FEAST simulation. The script defines a set of four custom components,
and combines those components to create a variety of sites. The sites are then combined into a gas field that is used in
the simulation. Three different survey-based LDAR programs are also defined and used in the simulation.
"""
import numpy as np
import copy
import feast.EmissionSimModules.infrastructure_classes
from feast.EmissionSimModules.infrastructure_classes import Component
import feast.DetectionModules as Dm

# The random seed below can be un-commented to generate reproducible realizations.
# np.random.seed(0)

# Defining component types in the gas field
# Generates reparable fugitive emissions
comp_fug = Component(
    name='Fugitive emitters',
    emission_data_path='production_emissions.p',
    emission_per_comp=0.00231,  # Fraction of components expected to be emitting at the beginning of the simulation.
    emission_production_rate=5.4 / 650 / 365  # number of new emissions per component per day
)

# Generates miscelaneous emissions that cannot be repaired (many tank emissions and pneumatic controller emissions can
# be modeled in this way).
misc_vent = Component(
    name='misc vents',
    emission_data_path='production_emissions.p',
    emission_per_comp=0.00231,  # Fraction of components expected to be emitting at the beginning of the simulation.
    emission_production_rate=5.4 / 650 / 365,  # number of new emissions per component per day
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


# Define the number of wells at each site
# Assign the number of wells per site
n_sites = 93  # This is the number of "basic pads" to be included in the simulation
# The well per site distribution was built using data from the Colorado Oil and Gas Conservation Commission
well_dist = np.cumsum([0.7842, 0.0762, 0.0382, 0.0280, 0.0150, 0.0114, 0.0082, 0.0111, 0.0047, 0.0045, 0.0035, 0.0038,
                      0.0027, 0.0016, 0.0018, 0.0026, 0.0010, 0.0011, 0.0006])
n_wells_odds = np.random.uniform(0, 1, n_sites)
n_wells_samp = np.zeros(n_sites, dtype=int)
for ind in range(n_sites):
    n_wells_samp[ind] = np.min(np.where(n_wells_odds[ind] < well_dist)[0]) + 1

counts, n_wells_bin = np.histogram(n_wells_samp,
                                   bins=np.linspace(0.5, np.max(n_wells_samp) + 0.5, np.max(n_wells_samp) + 1))


# ------Each iteration of this for loop generates one realization of the simulation
for ind in range(2):
    print('Iteration number: {:0.0f}'.format(ind))
    site_dict = {}
    # Assign components to sites
    for n_wells_in_site in range(1, np.max(n_wells_samp)):
        if counts[n_wells_in_site - 1] == 0:
            continue
        cond = (n_wells_samp == n_wells_in_site)
        basicpad = feast.EmissionSimModules.infrastructure_classes.Site(
            name='basic pad',
            comp_dict={
                # ------ The number of components is proportional to the number of wells, ind2
                'Fugitive ' + str(n_wells_in_site): {'number': 300 * n_wells_in_site,
                                                     'parameters': copy.copy(comp_fug)},
                'misc vent ' + str(n_wells_in_site): {'number': 350 * n_wells_in_site,
                                                      'parameters': copy.copy(misc_vent)}
            }
        )

        site_dict['basic pad ' + str(n_wells_in_site)] = {'number': counts[n_wells_in_site - 1], 'parameters': basicpad}
    plungersite = feast.EmissionSimModules.infrastructure_classes.Site(
        name='plunger well',
        comp_dict={
            'Fugitive_plunger': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'plunger': {'number': 2, 'parameters': copy.copy(plunger)},
            'misc vent plunge': {'number': 600, 'parameters': copy.copy(misc_vent)}
        }
    )
    unloadnp = feast.EmissionSimModules.infrastructure_classes.Site(
        name='Site 4',
        comp_dict={
            'Fugitive_4': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'no plunger': {'number': 2, 'parameters': copy.copy(noplunger)},
            'misc vent np': {'number': 600, 'parameters': copy.copy(misc_vent)}
        }
    )
    site_dict['plunger'] = {'number': 6, 'parameters': plungersite}
    site_dict['unload no plunger'] = {'number': 1, 'parameters': unloadnp}

    # Define time and gas field parameters
    timeobj = feast.EmissionSimModules.simulation_classes.Time(delta_t=1 / 24, end_time=3 * 365)
    gas_field = feast.EmissionSimModules.infrastructure_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    gas_field.met_data_path = 'TMY-DataExample.csv'
    gas_field.met_data_maker()

    # --- Define LDAR programs ----
    rep0 = Dm.repair.Repair(repair_delay=0)
    rep7 = Dm.repair.Repair(repair_delay=7)
    points = np.logspace(-3, 1, 100)  # emission rates
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.02)) / (0.8 * np.sqrt(2))) for f
                                  in points])
    probs[0] = 0
    ogi = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=50,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=copy.copy(rep0),
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,  # g/s
        detection_probabilities=probs
    )
    ogi_no_survey = Dm.comp_survey.CompSurvey(
        timeobj,
        survey_interval=None,
        survey_speed=150,
        ophrs={'begin': 8, 'end': 17},
        labor=100,
        dispatch_object=copy.copy(rep0),
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    points = np.logspace(-3, 1, 100)
    probs = 0.5 + 0.5 * np.array([np.math.erf((np.log(f) - np.log(0.474)) / (1.36 * np.sqrt(2))) for f
                                  in points])
    probs[0] = 0
    plane_survey = Dm.site_survey.SiteSurvey(
        timeobj,
        survey_interval=50,
        sites_per_day=200,
        site_cost=100,
        mu=0.1,
        dispatch_object=ogi_no_survey,
        detection_variables={'flux': 'mean'},
        detection_probability_points=points,
        detection_probabilities=probs
    )
    cm_ogi = copy.deepcopy(ogi_no_survey)
    cont_monitor = Dm.site_monitor.SiteMonitor(
        timeobj,
        time_to_detect_points=[[19, 1], [20, 1], [21, 1], [19, 5], [20, 5], [21, 5], [19, 5.1], [20, 5.1], [21, 5.1]],
        time_to_detect_days=[np.infty, 1, 0, np.infty, 5, 0, np.infty, np.infty, np.infty],
        detection_variables={'flux': 'mean', 'wind speed': 'mean'},
        dispatch_object=cm_ogi,
        site_queue=np.linspace(0, gas_field.n_sites - 1, gas_field.n_sites, dtype=int),
        op_envelope={
            'wind direction': {'class': 2, 'min': [[45, 225]]*gas_field.n_sites, 'max': [[135, 315]]*gas_field.n_sites}
        },
    )
    # Define LDAR programs
    ogi_survey = Dm.ldar_program.LDARProgram(
        timeobj, copy.deepcopy(gas_field), {'ogi': ogi}
    )
    # tiered survey
    tech_dict = {
        'plane': plane_survey,
        'ogi': ogi_no_survey
    }
    plane_ogi_survey = Dm.ldar_program.LDARProgram(
        timeobj, copy.deepcopy(gas_field), tech_dict
    )

    # continuous monitor
    tech_dict = {
        'cm': cont_monitor,
        'ogi': cm_ogi
    }
    cm_ogi = Dm.ldar_program.LDARProgram(
        timeobj, copy.deepcopy(gas_field), tech_dict
    )
    ldar_dict = {
        'cm': cm_ogi,
        'ogi': ogi_survey,
        'plane': plane_ogi_survey
    }
    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field,
        ldar_program_dict=ldar_dict,
        dir_out='ExampleRunScriptResults', display_status=False
    )
