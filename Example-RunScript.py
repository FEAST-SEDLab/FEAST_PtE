"""
This script provides an example of a custom FEAST simulation. The script defines a set of four custom components,
and combines those components to create a variety of sites. The sites are then combined into a gas field that is used in
the simulation. Three different survey-based LDAR programs are also defined and used in the simulation.
"""
import numpy as np
import copy
import feast
import feast.EmissionSimModules.infrastructure_classes
from feast.EmissionSimModules.infrastructure_classes import Component

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
for ind in range(30):
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

    # Define LDAR programs
    tech = {
        'ogi': {'survey_interval': 180,  # Days between surveys
                'mu': 0.0018,  # Emission rate with 50% probability of detection (g/s)
                'lam': 2.23,  # Width parameter for the probability of detection curve
                'survey_speed': 400  # Components per day
                },
    }
    tech_dict = {}
    for tech, params in tech.items():
        tech_dict[tech] = feast.DetectionModules.tech_detect.TechDetect(
            timeobj, gas_field,
            **params
        )
    tiered_tech = {
        'Plane': {
            'survey_interval': 180,  # Days between surveys
            'mu': 0.474,  # Emission rate with 50% probability of detection (g/s)
            'lam': 3.88,  # Width parameter for the probability of detection curve
            'sites_per_day': 222,  # Defines the speed of the primary survey (sites per day)
            'lam2': 2.23,  # Width parameter for the secondary survey probability of detection curve.
            'mu2': 0.002,  # Emission size with a 50% probability of detection for the secondary survey (g/s)
            'secondary_comps_hr': 400,  # Defines the speed of the secondary survey (components per hour)
            'site_cost': 100,
            'labor': 100
        }
    }

    for tech, params in tiered_tech.items():
        tech_dict[tech] = feast.DetectionModules.tiered_detect.TieredDetect(
            timeobj, gas_field,
            **params
        )

    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field,
        tech_dict=tech_dict, dir_out='ExampleRunScriptResults', display_status=False
    )
