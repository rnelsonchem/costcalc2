'''
Chemical Reaction Cost Cacluation Routines.
Adapted from the Excel spreadsheets prepared by Saeed Ahmad, PhD.
(C) Ryan Nelson
'''
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gspread

from oauth2client.client import GoogleCredentials
from google.colab import auth

from google.colab import files

# Set up some plotting stuff for the notebooks
plt.style.use('ggplot')
plt.rc('figure', dpi=150)

class GenericCost(object):
    # Not meant to be instantiated directly
    def rxn_cost(self, prod, amp=1.0):
        '''A recursive function for calculating the cost of an arbitrary
        reaction route.
        
        This is the workhorse function of the whole process.
        
        prod : str
            The name of the reaction to cost. This should also be the name of
            the final product for that reaction.

        amp : float 
            The number is an amplifier for masses of materials, for
            example. You shouldn't need to change this.  
        '''
        data = self.fulldata.loc[prod]

        # Kg of nonsolvent materials used per equivalent
        amount_kg = data['Equiv']*data['MW']#/(data['Density'])

        # Amount of solvent
        # First figure out which materials are solvents 
        mask = ~data['Volumes'].isna()
        # Which material are the volumes relative to? What is the kg?
        amt_rel = amount_kg[data.loc[mask, 'Relative']].values
        # Calculate the mass of solvent. Take into account the solvent recycyling 
        amount_kg[mask] = data.loc[mask, 'Volumes']*data.loc[mask, 'Density']*(1 - 
                data.loc[mask, 'Sol Recyc'])*amt_rel
        
        # Set the kg amounts in the large data table
        self.fulldata.loc[prod, 'kg/kg rxn'] = \
                (amount_kg/amount_kg[prod]).values

        # Calculate unknown costs. Looks for any empty values in the "Cost" column
        #unknown_cost = ~data['Cost calc'].isna()
        unknown_cost = data['Cost'].isna()
        # This is recursive. The final cost per kg of product will be
        # amplified by each subsequent step, which is where the new_amp
        # calculation comes into play
        for cpd, row in data.loc[unknown_cost].iterrows():
            # Don't do this for the product of the current reaction
            if cpd == prod: 
                continue
            # The amounts needed will be amplified by the appropriate kg ratio.
            # Set that ratio
            new_amp = data.loc[cpd, 'kg/kg rxn']
            # Run the cost calculation for the unknown compound
            cst = self.rxn_cost(cpd, amp*new_amp)
            # Set the calculated cost in the larger data table
            self.fulldata.loc[(prod, cpd), 'Cost'] = cst

        # Normalize the cost for each materials to 1 kg of reaction product
        self.fulldata.loc[prod, 'RM cost/kg rxn'] = \
                (data['kg/kg rxn']*data['Cost']).values
        # The product cost will be the sum of all the reactant/solvent costs
        self.fulldata.loc[(prod, prod), 'RM cost/kg rxn'] = \
                data['RM cost/kg rxn'].sum()
        
        # Calculate % costs for individual rxn
        # = (RM cost/kg rxn)/(RM cost/kg rxn for the rxn product)
        p_rm_cost = data['RM cost/kg rxn']*100/data.loc[prod, 'RM cost/kg rxn']
        self.fulldata.loc[prod, '% RM cost/kg rxn'] = p_rm_cost.values
        # Remove the % cost for the rxn product
        self.fulldata.loc[(prod, prod), '% RM cost/kg rxn'] = np.nan
        
        # These are the costs for ultimate product
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn cost
        # However, for feeder reactions this will get amplified by each step   
        self.fulldata.loc[prod, 'RM cost/kg prod'] = \
                self.fulldata.loc[prod, 'RM cost/kg rxn'].values*amp
        
        # Set the "Cost" to the calculated value
        self.fulldata.loc[(prod, prod), 'Cost'] = data.loc[prod, 'RM cost/kg rxn']
        
        # Return the calculated product cost, which is required for the recurisive 
        # nature of the algorithm
        return data.loc[prod, 'RM cost/kg rxn']

    
    def rxn_data_setup(self, ):
        '''Setup the full data set for the upcoming cost calculations. 
        
        For example, columns are added for caculated values. Some values are
        set to null.  Etc. This also cleans things up for a rerun of the cost
        analysis.  
        '''
        fulldata = pd.merge(self.materials, self.rxns, on='Compound', how='right')
        if fulldata['MW'].isna().any():
            # raise ValueError('You are missing a material from a rxn.')
            print('There is a mismatch between the reaction and materials file.')
            print('Material missing from materials sheet.')

        fulldata.set_index(['Prod', 'Compound'], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not copies
        fulldata = fulldata.sort_index()
        
        # Set the costs to NaN for materials that will have costs calculated 
        cost_recalc_mask = ~fulldata['Cost calc'].isna()
        fulldata.loc[cost_recalc_mask, 'Cost'] = np.nan
        
        # Create a bunch of columns that will be populated later during cost
        # calculation
        #fulldata['Amount_kg'] = np.nan
        fulldata['kg/kg rxn'] = np.nan
        #fulldata['RM cost/eqivs'] = np.nan
        fulldata['RM cost/kg rxn'] = np.nan
        fulldata['% RM cost/kg rxn'] = np.nan
        fulldata['RM cost/kg prod'] = np.nan
        fulldata['% RM cost/kg prod'] = np.nan
        
        self.fulldata = fulldata


    def rxn_data_post(self,):
        '''Calculate the % costs of the raw materials for a route.
        '''
        prod = self.final_prod
        # Calculate % overall costs relative to the prod
        prod_cost = self.fulldata.loc[(prod, prod), 'RM cost/kg rxn'] 
        self.fulldata['% RM cost/kg prod'] = \
                self.fulldata['RM cost/kg prod']*100/prod_cost
        
        # Remove the cost and %s for cost-calculated materials
        # This is necessary so that this column adds up to 100%
        mask = ~self.fulldata['Cost calc'].isna()
        self.fulldata.loc[mask, '% RM cost/kg prod'] = np.nan
        # This filters some of the costs which are simply the sum of raw materials
        # from eariler rxns. The sum of this column will now be equal to the cost
        # of the final product.
        self.fulldata.loc[mask, 'RM cost/kg prod'] = np.nan
     

class ColabCost(GenericCost):
    def __init__(self, materials_sheet_key, rxn_sheet_key, final_prod,
            materials_worksheet=0, rxn_worksheet=0, alt_mat_key=None,
            alt_mat_sheet=0):
        # Authenticate the Colab environment 
        auth.authenticate_user()
        self._gc = gspread.authorize(GoogleCredentials.get_application_default())
        
        # Set up the reaction DataFrame
        self._rxn_sheet_key = rxn_sheet_key
        self._rxn_worksheet = rxn_worksheet
        self.final_prod = final_prod
        self._rxn_read()
        
        # Create the Materials DataFrame from a main sheet and an optional
        # alternate sheet.
        self._material_sheet_key = materials_sheet_key
        self._materials_worksheet = materials_worksheet
        self._alt_mat_key = alt_mat_key
        self._alt_mat_sheet = alt_mat_sheet
        self._materials_build()                

        # Combine the reaction/materials sheets, add the new columns
        self.rxn_data_setup()
        # Run the costing/post processing
        self.cost = self.rxn_cost(final_prod)
        self.rxn_data_post()
        # Make a copy of the costing data for plotting
        self._fdcopy = self.fulldata.copy()
        
    def _materials_build(self, ):
        '''This function combines the materials DataFrames.'''
        materials = self._materials_read(self._materials_key,
                self._materials_sheet)

        # If an alternative materials key is given, combine that materials
        # sheet with the main one
        if self._alt_mat_key:
            alt_mats = self._materials_read(self._alt_mat_key,
                    self._alt_mat_sheet)
            # Concatenate the sheets. Reset the index so that it is
            # consecutively numbered
            materials = pd.concat([materials, alt_mats])\
                    .reset_index(drop=True)
            # If the two materials sheets have any of the same materials, it
            # could cause some problems. Throw an error if this is the case
            if materials.duplicated('Compound').any():
                raise ValueError('Duplicated materials will cause errors')

        # Set the final materials sheet
        self.materials = materials

    def _materials_read(self, mat_key, wsheet):
        '''Read a materials Google sheet.'''
        # Grab the Google sheet handle, pull down all values and make a DataFrame
        mat_sh = self._gc.open_by_key(mat_key)
        ws = mat_sh.get_worksheet(wsheet)
        vals = ws.get_all_values()
        mats = pd.DataFrame(data=vals[1:], columns=vals[0])
        
        # Convert empty cells to NaN
        mask = (mats == '')
        mats[mask] = np.nan
        # Drop empty rows
        mats.dropna(how='all', inplace=True)
        
        # Convert numeric/date columns
        mats['MW'] = pd.to_numeric(mats['MW'])
        mats['Density'] = pd.to_numeric(mats['Density'])
        mats['Cost'] = pd.to_numeric(mats['Cost'])
        mats['Date'] = pd.to_datetime(mats['Date'])
        
        return mats
        
    def _rxn_read(self, ):
        '''Read a Google Sheet of rxn info and merge with materials.'''
        rxn_sh = self._gc.open_by_key(self._rxn_sheet_key)
        ws = rxn_sh.get_worksheet(self._rxn_worksheet)
        vals = ws.get_all_values()
        rxns = pd.DataFrame(data=vals[1:], columns=vals[0])
        
        # Convert empty cells to NaN
        mask = rxns == ''
        rxns[mask] = np.nan
        # Drop empty rows
        rxns.dropna(how='all', inplace=True)
        
        rxns['Equiv'] = pd.to_numeric(rxns['Equiv'])
        rxns['Volumes'] = pd.to_numeric(rxns['Volumes'])
        rxns['Sol Recyc'] = pd.to_numeric(rxns['Sol Recyc'])
        
        self.rxns = rxns
        
    def excel_download(self, fname):
        self.fulldata.to_excel(fname)
        # There seems to be a bit of a lag before you can download
        # the file, this delay might fix some of the errors this causes
        time.sleep(2)
        files.download(fname)

    def value_mod(self, cpd, vals, scan_type='Cost', step=None):
        '''Manually set/scan the cost/equiv of a given material.'''
        # If a single value was given, convert to a list
        # Set this flag to undo the list at the end of the function
        val_list = True
        if isinstance(vals, (float, int)):
            vals = [vals,]
            val_list = False
        
        all_costs = []
        for val in vals:
            # You'll need to reset some of the columns, luckily we have a copy
            # of the main DataFrame
            self.rxn_data_setup()
           
            # Set the value. The first one sets all values w/ the compound
            # name. The second one sets only a the values for a specific
            # reaction.
            if not step:
                self.fulldata.loc[(slice(None), cpd), scan_type] = val
            else: 
                self.fulldata.loc[(step, cpd), scan_type] = val
                
            # Cost this modified version, append the cost to the list
            int_cost = self.rxn_cost(self.final_prod)
            all_costs.append(int_cost)
        
        # Make a modified DataFrame copy for later
        self.moddata = self.fulldata.copy()
        # Reset the original costing DataFrame
        self.fulldata = self._fdcopy.copy()

        # When a single value was used, return just that one value. Otherwise,
        # a list will be returned
        if val_list == False:
            all_costs = all_costs[0]
        
        return all_costs
