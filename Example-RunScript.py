import numpy as np
import copy
import feast
from feast.GeneralClassesFunctions.simulation_classes import Component

# np.random.seed(0)


# Generates reparable fugitive emissions
comp_fug = Component(
    name='Fugitive emitters',
    leak_data_path='production-emissions.p',
    leaks_per_comp=0.00231,
    leak_production_rate=5.4 / 650 / 365
)

# Generates miscelaneous emissions that cannot be repaired
misc_vent = Component(
    name='misc vents',
    leak_data_path='production-emissions.p',
    leaks_per_comp=0.00231,
    leak_production_rate=5.4 / 650 / 365,
    base_reparable=False,
)


plunger = feast.GeneralClassesFunctions.simulation_classes.Component(
    name='plunger',
    leak_production_rate=0,
    episodic_emission_sizes=[50.36],  # gram-per-sec
    episodic_emission_duration=2059 / 3600 / 24,  # see Zaimes Figure S5
    episodic_emission_per_day=0.0165  # totals['n_plunger_events'] / totals['n_plunger_wells'] / 365
)

noplunger = feast.GeneralClassesFunctions.simulation_classes.Component(
    name='no plunger',
    leak_production_rate=0,
    episodic_emission_sizes=[76.98],  # [np_methane / np_events / 4973.5 * 1e6], # gram-per-sec
    episodic_emission_duration=4973.5 / 3600 / 24,  # see Zaimes Figure S3
    episodic_emission_per_day=0.0002111  # np_events / np_wells / 365
)


# Assign the number of wells per site
n_sites = 93
# The well per site distribution was built using data from the Colorado Oil and Gas Conservation Commission
well_dist = np.cumsum([0.7842, 0.0762, 0.0382, 0.0280, 0.0150, 0.0114, 0.0082, 0.0111, 0.0047, 0.0045, 0.0035, 0.0038,
                      0.0027, 0.0016, 0.0018, 0.0026, 0.0010, 0.0011, 0.0006])
n_wells_odds = np.random.uniform(0, 1, n_sites)
n_wells_samp = np.zeros(n_sites)
for ind in range(n_sites):
    n_wells_samp[ind] = np.min(np.where(n_wells_odds[ind] < well_dist)[0]) + 1

counts, n_wells_bin = np.histogram(n_wells_samp,
                                   bins=np.linspace(0.5, np.max(n_wells_samp) + 0.5, np.max(n_wells_samp) + 1))
    

for ind in range(1):
    print('Iteration number: {:0.0f}'.format(ind))
    site_dict = {}
    # Assign a number of components to each site based on the number of wells
    for ind2 in range(1, np.max(n_wells_samp)):
        if counts[ind2 - 1] == 0:
            continue
        cond = (n_wells_samp == ind2)
        basicpad = feast.GeneralClassesFunctions.simulation_classes.Site(
            name='basic pad',
            comp_dict={
                'Fugitive ' + str(ind2): {'number': 300 * ind2, 'parameters': copy.copy(comp_fug)},
                'misc vent ' + str(ind2): {'number': 350 * ind2, 'parameters': copy.copy(misc_vent)}
            }
        )

        site_dict['basic pad' + str(ind2)] = {'number': counts[ind2 - 1], 'parameters': basicpad}
    plungersite = feast.GeneralClassesFunctions.simulation_classes.Site(
        name='plunger well',
        comp_dict={
            'Fugitive_plunger': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'plunger': {'number': 2, 'parameters': copy.copy(plunger)},
            'misc vent plunge': {'number': 600, 'parameters': copy.copy(misc_vent)}
        }
    )
    unloadnp = feast.GeneralClassesFunctions.simulation_classes.Site(
        name='Site 4',
        comp_dict={
            'Fugitive_4': {'number': 700, 'parameters': copy.copy(comp_fug)},
            'no plunger': {'number': 2, 'parameters': copy.copy(noplunger)},
            'misc vent np': {'number': 600, 'parameters': copy.copy(misc_vent)}
        }
    )
    site_dict['plunger'] = {'number': 6, 'parameters': plungersite}
    site_dict['unload no plunger'] = {'number': 1, 'parameters': unloadnp}
    timeobj = feast.GeneralClassesFunctions.simulation_classes.Time(delta_t=1/24, end_time=3*365)
    gas_field = feast.GeneralClassesFunctions.simulation_classes.GasField(
        sites=site_dict,
        time=timeobj
    )
    tech = {
        'ogi': {'survey_interval': 180, 'mu': 0.0018,
                'lam': 2.23, 'survey_speed': 400},
        'drone': {'survey_interval': 180, 'mu': 0.0057,
                  'lam': 4.52, 'survey_speed': 983}
    }
    tech_dict = {}
    for tech, params in tech.items():
        tech_dict[tech] = feast.DetectionModules.tech_detect.TechDetect(
            timeobj, gas_field,
            **params
        )
    tiered_tech = {
        'Plane': {
            'survey_interval': 180, 'mu': 0.474,
            'lam': 3.88, 'sites_per_day': 222,
            'lam2': 2.23, 'mu2': 0.002,
            'secondary_comps_hr': 400}
    }
  
    for tech, params in tiered_tech.items():
        tech_dict[tech] = feast.DetectionModules.tiered_detect.TieredDetect(
            timeobj, gas_field,
            **params
        )

    feast.field_simulation.field_simulation(
        time=timeobj, gas_field=gas_field,
        tech_dict=tech_dict, dir_out='Results', display_status=True
    )
