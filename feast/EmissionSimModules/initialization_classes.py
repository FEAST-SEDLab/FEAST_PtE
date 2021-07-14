import pathlib
import copy
import math

import pandas as pd
import numpy as np

import feast.EmissionSimModules.infrastructure_classes
from .infrastructure_classes import Component


class SiteInitializer:
    """
    Class to read model formulation data from spreadsheets and generate a dictionary of sites
    """
    def __init__(self, gfield_file, feast_data_path, returncomps=False, returnequip=False, returnsites=False):
        """
        :param gfield_file: spreadsheet (.xlsx) file containing gas-field, site-type, equipment-type, and component-type
        parameters.
        :feast_data_path: path to location of feast data (i.e. repair cost data and emissions distributions)
        :param returncomps: option to return a dictionary of components
        :param returnequip: option to return a dictionary of equipment
        :param returnsites: option to return a dictionary of sites
        """
        self.gfield_params = pd.read_excel(gfield_file, sheet_name='Gas_Field')
        self.site_params = pd.read_excel(gfield_file, sheet_name='Site')
        self.equip_params = pd.read_excel(gfield_file, sheet_name='Equipment')
        self.comp_params = pd.read_excel(gfield_file, sheet_name='Components')
        self.feast_data_path = feast_data_path
        self.returncomps = returncomps
        self.returnequip = returnequip
        self.returnsites = returnsites
        self.components = {}
        self.equipment = {}
        self.sitedict = {}

    def component_maker(self):
        """
        Method to read component parameters from the spreadsheet and instantiate component objects for components that
        will be used during the simulation. Equipment used is derived from the "Equipment" sheet in the xlsx workbook.
        """
        # Identify Components used during the simulation
        site_types_used = self.gfield_params.site_type.unique()
        equipment_types_used = self.site_params.loc[self.site_params.Site_type.isin(site_types_used),
                                                    'Equipment_type'].unique()
        component_types_used = self.equip_params.loc[self.equip_params.Equipment_Type.isin(equipment_types_used),
                                                     'Component_type'].unique()

        # Instantiate the component dictionary
        self.components = {comp: None for comp in component_types_used}

        for comp in self.components.keys():
            # LOAD COMPONENT PARAMETERS FROM FILE

            # Emissions distribution ---------------------------------- #
            emiss_dist = self.comp_params.loc[self.comp_params.Component_type == comp, 'emissions_distribution'].unique()[0]
            if isinstance(emiss_dist, str):
                if emiss_dist == 'None':
                    emiss_dist_file = None
                else:
                    emiss_dist_file = list(pathlib.Path(f"{self.feast_data_path}").glob(f"{emiss_dist}.p"))[0]
            elif math.isnan(emiss_dist):
                emiss_dist_file = None
            else:
                raise Exception("Unknown value passed for emissions distribution. Check component parameters sheet.")

            # Repair cost data ---------------------------------------------- #
            rep_cost = self.comp_params.loc[self.comp_params.Component_type == comp, 'Repair_costs'].unique()[0]
            if isinstance(rep_cost, str):
                if rep_cost == 'None':
                    rep_cost_file = None
                else:
                    rep_cost_file = list(pathlib.Path(f"{self.feast_data_path}").glob(f"{rep_cost}.p"))[0]
            elif math.isnan(rep_cost):
                rep_cost_file = None
            else:
                raise Exception("Unknown value passed for repair costs. Check component parameters sheet.")

            # Reparable status (yes/no - True/False) ------------------------------- #
            rep_cat = self.comp_params.loc[self.comp_params.Component_type == comp, 'reparable'].unique()[0]
            if rep_cat == 'yes':
                rep_status = True
            elif rep_cat == 'no':
                rep_status = False
            else:
                raise Exception(
                    'Repair category must be either yes or no. Reconfigure Parameter sheet and re-run simulation.'
                )

            # Emission Production Rate ------------------------------------------- #
            emission_prod_rate = self.comp_params.loc[self.comp_params.Component_type == comp,
                                                      'emission_production_rate'].unique()[0]
            if math.isnan(emission_prod_rate):
                emission_prod_rate = 0

            # Emission per component ---------------------------------------------- #
            emission_per_comp = self.comp_params.loc[self.comp_params.Component_type == comp,
                                                     'emission_per_component'].unique()[0]
            if math.isnan(emission_per_comp):
                emission_per_comp = None

            # Episodic emission sizes ------------------------------------------- #
            episodic_sizes = self.comp_params.loc[self.comp_params.Component_type == comp,
                                                  'episodic_emission_sizes'].unique()[0]
            if isinstance(episodic_sizes, str):
                episodic_sizes = episodic_sizes[1:]
                es = episodic_sizes.split(',')
                es = [float(i) for i in es]
            elif np.isnan(episodic_sizes):
                es = [0]
            elif isinstance(episodic_sizes, float):
                es = [episodic_sizes]
            else:
                raise Exception("""
                    Unsupported value entered for Episodic Emission Sizes. Check the component parameters sheet.
                    Values should be an integer, comma separated integers (no spaces, formatted as string), or the cell 
                    should be left blank.
                """)

            # Episodic emissions per day ---------------------------------------- #
            episodic_emiss_per_day = self.comp_params.loc[self.comp_params.Component_type == comp,
                                                          'episodic_emission_per_day'].unique()[0]
            if math.isnan(episodic_emiss_per_day):
                episodic_emiss_per_day = 0

            # Episodic emission duration ---------------------------------------- #
            episodic_emiss_duration = self.comp_params.loc[self.comp_params.Component_type == comp,
                                                           'episodic_emission_duration'].unique()[0]
            if math.isnan(episodic_emiss_duration):
                episodic_emiss_duration = 0

            # Vent sizes -------------------------------------------------------- #
            vent_sizes = self.comp_params.loc[self.comp_params.Component_type == comp, 'vent_sizes'].unique()[0]
            if isinstance(vent_sizes, str):
                vent_sizes = vent_sizes[1:]
                vs = vent_sizes.split(',')
                vs = [float(i) for i in vs]
            elif np.isnan(vent_sizes):
                vs = [0]
            elif isinstance(vent_sizes, float):
                vs = vent_sizes
            else:
                raise Exception("""
                    Unsupported value entered for Episodic Emission Sizes. Check the component parameters sheet.
                    Values should be an integer, comma separated integers (no spaces, formatted as string), or the cell 
                    should be left blank.
                                """)

            # Vent period ------------------------------------------------------- #
            vent_period = self.comp_params.loc[self.comp_params.Component_type == comp, 'vent_period'].unique()[0]
            if math.isnan(vent_period):
                vent_period = np.infty

            # Vent Starts ------------------------------------------------------- #
            vent_starts = self.comp_params.loc[self.comp_params.Component_type == comp, 'vent_starts'].unique()[0]
            if np.isnan(vent_starts):
                vent_starts = np.array([])
            else:
                vent_starts = str(vent_starts)
                vent_starts = vent_starts.split(',')
                vent_starts = [float(i) for i in vent_starts]

            # Vent duration ----------------------------------------------------- #
            vent_duration = self.comp_params.loc[self.comp_params.Component_type == comp, 'vent_duration'].unique()[0]
            if math.isnan(vent_duration):
                vent_duration = 0

            # INSTANTIATE THE COMPONENT OBJECT
            locals()[comp] = Component(
                name=comp,
                repair_cost_path=rep_cost_file,
                base_reparable=rep_status,
                emission_data_path=emiss_dist_file,
                emission_production_rate=emission_prod_rate,
                emission_per_comp=emission_per_comp,
                episodic_emission_sizes=es,
                episodic_emission_per_day=episodic_emiss_per_day,
                episodic_emission_duration=episodic_emiss_duration,
                vent_sizes=vs,
                vent_period=vent_period,
                vent_starts=vent_starts,
                vent_duration=vent_duration
            )
            self.components[comp] = locals()[comp]

        if self.returncomps:
            return self.components

    def equipment_maker(self, component_dict=None):

        # Generate components or pull-in components from dictionary
        if len(self.components) == 0:
            if component_dict is not None:
                self.components = component_dict
            else:
                self.component_maker()

        # Get equipment types used in the simulation
        site_types_used = self.gfield_params.site_type.unique()
        equipment_types_used = self.site_params.loc[self.site_params.Site_type.isin(site_types_used),
                                                    'Equipment_type'].unique()

        # Instantiate the equipment dictionary
        self.equipment = {equip: {} for equip in equipment_types_used}

        # Populate the equipment dictionary
        for equip in equipment_types_used:
            equip_components = self.equip_params.loc[
                self.equip_params.Equipment_Type == equip, 'Component_type'
            ].values
            for comp in equip_components:
                name = f"{equip}__{comp}"
                self.equipment[equip].update({name: {'number': None, 'parameters': None}})
                comp_ct = self.equip_params.loc[(self.equip_params.Equipment_Type == equip) &
                                                (self.equip_params.Component_type == comp), 'Component_Count'].values[0]
                self.equipment[equip][name]['number'] = comp_ct
                self.equipment[equip][name]['parameters'] = copy.deepcopy(
                    self.components[comp]
                )

        if self.returnequip:
            return self.equipment

    def site_maker(self, component_dict=None, equipment_dict=None):

        if len(self.components) == 0:
            if component_dict is not None:
                self.components = component_dict
            else:
                self.component_maker()
        if len(self.equipment) == 0:
            if equipment_dict is not None:
                self.equipment = equipment_dict
            else:
                self.equipment_maker()

        # Get site types used in the simulation
        site_types_used = self.gfield_params.site_type.unique()

        for site_type in site_types_used:
            site_ct = self.gfield_params.loc[self.gfield_params.site_type == site_type, 'site_count'].values[0]
            for site_num in range(site_ct):
                site_equip = self.site_params.loc[self.site_params.Site_type == site_type, 'Equipment_type'].values
                site_name = self.gfield_params.loc[self.gfield_params.site_type == site_type, 'Site_name'].values[0]
                self.sitedict.update({f"{site_name}_{site_num + 1}": {'number': 1, 'parameters': None}})

                equip_dict = {}
                for equip in site_equip:
                    equip_ct = int(self.site_params.loc[(self.site_params.Site_type == site_type) &
                                                        (self.site_params.Equipment_type == equip),
                                                        'Equipment_Count'].values[0])
                    comp_types = self.equipment[equip]

                    for comp in comp_types:
                        comp_types[comp]['number'] = comp_types[comp]['number'] * equip_ct
                        equip_dict.update({f"{site_name}_{site_num + 1}__{comp}": comp_types[comp]})
                        # for i in range(equip_ct):
                        #     equip_dict.update({f"{site_name}_{comp}_{i}": self.equipment[equip][comp]})
                self.sitedict[f"{site_name}_{site_num + 1}"]['parameters'] =\
                    feast.EmissionSimModules.infrastructure_classes.Site(
                    name=f"{site_name}_{site_num + 1}",
                    comp_dict=equip_dict
                )
        if self.returnsites:
            return self.sitedict
