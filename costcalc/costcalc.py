'''
Chemical Reaction Cost Cacluation Routines.
Adapted from the Excel spreadsheets prepared by Saeed Ahmad, PhD.
(C) Ryan Nelson
'''
import time
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gspread


# Set up some plotting stuff for the notebooks
plt.style.use('ggplot')
plt.rc('figure', dpi=150)

# For future reference -- Can't use this in conjunction with "fillna" becuase
# the DF doesn't display correctly 
# Set Pandas precision
#pd.set_option('precision', 2)

# Setting the display to print function
# This gets changed for Jupyter Notebook sessions
disp = print

class GenericCost(object):
    '''Base costing class. 

    Not meant to be instantiated directly, but defines some common functions
    for any subclasses.

    Parameters
    ----------
    final_prod : str
        Defines the final product name for costing calculations.
    '''
    def __init__(self, final_prod):
        # This is a common attribute for all subclasses.
        self.final_prod = final_prod        

    def rxn_data_setup(self, ):
        '''Setup the full data set for the upcoming cost calculations. 
        
        This method merges the materials and reaction sheets into a combined
        DataFrame called `fulldata`. This function also serves as a "reset" of
        sorts after a costing calculation, if you want to start over with a
        different costing model.
        '''
        mat_drops = ['Notes', 'CAS Num']
        rxn_drops = ['Notes',]
        fulldata = pd.merge(self.materials.drop(mat_drops, axis=1), 
                            self.rxns.drop(rxn_drops, axis=1), 
                            on='Compound', how='right')
        
        # Check for a missing material -- Everything should have a MW
        if fulldata['MW'].isna().any():
            # raise ValueError('You are missing a material from a rxn.')
            print('There is a mismatch between the reaction and materials file.')
            print('Material might be missing from materials sheet.')
            
        # Check for a missing cost, which is not being calculated
        # This is a tricky error because the costing will run just fine with 
        # NaNs. Check for this by looking for NaN in both Cost and Calc columns
        if (fulldata['Cost calc'].isna() & fulldata['Cost'].isna()).any():
            print('You are missing a material cost that is not being calculated.')

        fulldata.set_index(['Prod', 'Compound'], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not copies
        fulldata = fulldata.sort_index()
        
        # Save the full data set
        self.fulldata = fulldata

        # Add a modified variable placeholder. This will store modified
        # values for later processing
        self._mod_vals = []
        # Add the empty columns
        self._column_clear()

    def _column_clear(self, ):
        '''Clear out calculated values.

        This method will reset all the calculated values in the `fulldata`
        DataFrame. This is perhaps not strictly necessary, but it should help
        to avoid unwanted errors in recalculations due to old data still being
        present. This is not the same as a reset, though, because manually
        modified values with `value_mod` method will be re-modified. 
        '''
        # Set the costs to NaN for materials that will have costs calculated 
        cost_recalc_mask = ~self.fulldata['Cost calc'].isna()
        self.fulldata.loc[cost_recalc_mask, 'Cost'] = np.nan

        # Modify stored mod variables
        for mod in self._mod_vals:
            self._set_val(*mod)

        # Create or clear a bunch of columns that will be populated during 
        # cost calculation. 
        empty_cols = ['kg/kg rxn', 'RM cost/kg rxn', '% RM cost/kg rxn',
                'kg/kg prod', 'RM cost/kg prod', '% RM cost/kg prod',
                ]
        for col in empty_cols:
            self.fulldata[col] = np.nan

    def value_mod(self, cpd, val, val_type='Cost', step=None):
        '''Manually set a value for a given material.

        Parameters
        ----------
        cpd : str
            This the compound name for which the value will be modified.

        val : int, float
            This is the modified value for the parameter.

        val_type : str, optional (Default = 'Cost')
            This is the column name for the parameter that you'll be changing.
            This must be for a non-calculated column, such as 'Cost', 'Equiv',
            'OPEX', etc.

        step : None, str, optional (Default = None)
            The name of the reaction step for which this value will be
            changed. If this is `None` (default), then all the values for the
            given compound (`cpd`) will be set to the same value. This is
            mostly important for something like `val_type`='Equiv'. Clearly,
            you would only want to change the number of equivalents for a
            specific reaction. If this parameter is left as `None`, the
            equivalents for a given compound in all reactions will be set to
            the same value.
        
        Note
        ----
        This method will *NOT* recalculate the cost; this must be done as a
        separate step.
        '''
        # Store the values
        self._mod_vals.append( (cpd, val, val_type, step) )
        # This will clear out the old calculation data and set the modified
        # value. Keeps folks from getting confused b/c calculated values are
        # unchanged.
        self._column_clear()
            
    def _set_val(self, cpd, val, val_type, step):
        '''Set a modified value in the `fulldata` DataFrame
        '''
        # The first one sets all values w/ the compound name. The second one
        # sets only a value for a specific reaction.
        if not step:
            cells = (slice(None), cpd)
        else: 
            cells = (step, cpd)
            
        self.fulldata.loc[cells, val_type] = val
        # The "Cost calc" flag must be set to np.nan when setting a cost. 
        # This is necessary for % RM cost calcs, e.g.
        if val_type == 'Cost':
            self.fulldata.loc[cells, 'Cost calc'] = np.nan

    def value_scan(self, cpd, vals, val_type='Cost', step=None):
        '''Scan a range of values for a given material.
        
        Parameters
        ----------
        See `value_mod` method description, except for the following.

        vals : container of int/float values, int, float
            This is the container of values for which to scan through. If you
            want, this can be a single value, although the `value_mod` method
            may be more appropriate for that. 

        Returns
        -------
        list of floats
            This is the costs associate with each value in the input
            container. 

        Notes
        -----
        Although this method recalculates the cost for every value, it does
        not modify the original `fulldata` attribute. However, the `cost`
        attribute will reflect the final cost calculation.
        '''
        # If a single value was given, convert to a list
        # Set this flag to undo the list at the end of the function
        val_list = True
        if isinstance(vals, (float, int)):
            vals = [vals,]
            val_list = False
       
        # I need a copy of the full data set in order to reset for each
        # iteration. Otherwise, I was noticing some issues.
        fd_copy = self.fulldata.copy()
        all_costs = []
        for val in vals:
            self.value_mod(cpd, val, val_type, step)
            self.calc_cost()
            all_costs.append(self.cost)
            # Reset the full data set 
            self.fulldata = fd_copy.copy()
            # Pop out the mod value, otherwise this list will get really long
            self._mod_vals.pop()

        # Reset the final cost
        self.cost = self.fulldata.loc[(self.final_prod, self.final_prod), 
                                  'RM cost/kg rxn']
        
        # When a single value was used, return just that one value. Otherwise,
        # a list will be returned
        if val_list == False:
            all_costs = all_costs[0]
        
        return all_costs

    def calc_cost(self, ):
        '''Calculate the cost of the route. 
        
        This function combines the column clearing method, the reaction
        costing algorithm, and a post processing function.
        '''
        # Save a time stamp so it can be displayed later
        self._now = datetime.now().strftime('%Y-%m-%d')
        # Prep the DataFrame
        self._column_clear()
        # Run the costing and set the cost attribute
        self.cost = self._rxn_cost(self.final_prod)
        # Post process the DataFrame
        self._rxn_data_post()
        
    def _rxn_cost(self, prod, amp=1.0):
        '''The recursive cost calculating function. 
        
        This is the workhorse function of the whole process, but is not meant
        to be called on its own. Typically, you'll want to call `calc_cost` to
        get the costing information.
        
        Parameters
        ----------
        prod : str
            The name of the reaction to cost. This should also be the name of
            the final product for that reaction.

        amp : float, optional (default = 1.0)
            This number is an amplifier that increases some of the values,
            e.g. masses of materials, based on how much material is being used
            for the overall final product.   
        '''
        # Select out the reaction of interest from the full data set. Saves
        # some typing.
        data = self.fulldata.loc[prod]

        # Kg of nonsolvent materials used per equivalent
        amount_kg = data['Equiv']*data['MW']#/(data['Density'])

        # Amount of solvent
        # First figure out which materials are solvents 
        mask = ~data['Volumes'].isna()
        # Which material are the volumes relative to? What is the kg?
        amt_rel = amount_kg[data.loc[mask, 'Relative']].values
        # Calculate the mass of solvent. Take into account the solvent 
        # recycyling 
        # kg sol = Volume*Density*(1-Recycle)*(kg SM)
        amount_kg[mask] = data.loc[mask, 'Volumes']*data.loc[mask, 'Density']\
            *(1 - data.loc[mask, 'Sol Recyc'])*amt_rel
        
        # Set the kg/rxn amounts in the large data table. This is normalized
        # to make the product kg = 1
        self.fulldata.loc[prod, 'kg/kg rxn'] = \
                (amount_kg/amount_kg[prod]).values

        # Calculate unknown costs. Looks for any empty values in the "Cost" 
        # column. Don't use the 'Cost calc' column directly, because some of
        # the costs may have been manually set using the `value_mod` method
        # unknown_cost = ~data['Cost calc'].isna()
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
            cst = self._rxn_cost(cpd, amp*new_amp)
            # Set the calculated cost in the larger data table
            self.fulldata.loc[(prod, cpd), 'Cost'] = cst

        # Calculate the cost for each material in the reaction
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
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn 
        # cost. However, for feeder reactions this will get amplified by 
        # each step   
        self.fulldata.loc[prod, 'RM cost/kg prod'] = \
                self.fulldata.loc[prod, 'RM cost/kg rxn'].values*amp
        
        # Set the "Cost" to the calculated value
        self.fulldata.loc[(prod, prod), 'Cost'] = \
                data.loc[prod, 'RM cost/kg rxn']

        # This sets the number of kg of each material per kilogram of product
        # This is done by multiplying the per reaction value by the amplifier
        # This isn't necessary for costing, but we can use it for PMI
        self.fulldata.loc[prod, 'kg/kg prod'] = \
                self.fulldata.loc[prod, 'kg/kg rxn'].values*amp
        
        # Return the calculated product cost, which is required for the 
        # recurisive nature of the algorithm. In addition, an optional OPEX
        # may be added to take into acount production costs of the cpd
        if np.isnan(data.loc[prod, 'OPEX']):
            return data.loc[prod, 'RM cost/kg rxn']
        else:
            return data.loc[prod, 'RM cost/kg rxn'] + data.loc[prod, 'OPEX']
    
    def _rxn_data_post(self,):
        '''Calculate some values after the final cost of the reaction is
        determined. 
        '''
        prod = self.final_prod
        
        # If an OPEX for the final reaction is given, add that to the cost
        # of the final product
        opex = self.fulldata.loc[(prod, prod), 'OPEX']
        if not np.isnan(opex):
            self.fulldata.loc[(prod, prod), 'Cost'] = self.cost
                
        # Calculate % overall costs relative to the prod
        self.fulldata['% RM cost/kg prod'] = \
                self.fulldata['RM cost/kg prod']*100/self.cost
        
        # Filter out certain values to simplify full data set
        # Remove the cost and %s for cost-calculated materials
        # This is necessary so that this column adds up to 100% (w/o OPEX)
        mask = ~self.fulldata['Cost calc'].isna()
        self.fulldata.loc[mask, '% RM cost/kg prod'] = np.nan
        # This filters some of the costs which are simply the sum of raw materials
        # from eariler rxns. The sum of this column will now be equal to the cost
        # of the final product.
        self.fulldata.loc[mask, 'RM cost/kg prod'] = np.nan
        # This filters out the kg/kg prod values that were calculated, so that
        # the sum of this column is the PMI
        self.fulldata.loc[mask, 'kg/kg prod'] = np.nan
        # But we are making 1 kg of final product so that needs to be reset
        self.fulldata.loc[(prod, prod), 'kg/kg prod'] = 1.

        # PMI Calculations
        # Need to append this prefix for sorting purposes
        # There will be a funny column in this DF with these values...
        self._pre = 'zzzz'
        
        # First of all, calculate the PMI for each reaction individually
        gb = self.fulldata[['kg/kg rxn']].groupby('Prod')
        rxn_pmi = gb.sum().reset_index()
        rxn_pmi['Compound'] = self._pre + rxn_pmi['Prod'] + ' PMI'
        
        # The full route PMI is not the sum of the above, but is the sum of
        # the 'kg/kg prod' column. We need to make this into a DataFrame to
        # merge with the per reaction values above
        df_vals = {'kg/kg prod': [self.fulldata['kg/kg prod'].sum()], 
                   'Prod': [self.final_prod],
                   'Compound': [self._pre*2 + 'Full Route PMI']
                   }
        full_pmi = pd.DataFrame(df_vals)

        # Merge the per-reaction and full PMI
        self.pmi = pd.concat([rxn_pmi, full_pmi], 
                             sort=False).set_index('Prod')

    def results(self, style='compact', decimals=2, fill='-'):
        '''Print all the results of the costing calculation.
        '''
        # Print the time the calculation was run
        print('As of', self._now, '--')
        
        # Print a string about the final cost of the product
        if decimals:
            dec_str = ':.{:d}f'.format(decimals)
        else:
            dec_str = ':f'
        cost_str = 'The final cost of {} is ${' + dec_str + '}/kg.'
        print(cost_str.format(self.final_prod, self.cost))
            
        # For compact display, these are the most important columns
        comp_col = ['Cost', 'Equiv', 'Volumes', 'Sol Recyc', 'OPEX', 
                    'kg/kg rxn', 'RM cost/kg rxn', '% RM cost/kg rxn',
                    'kg/kg prod', 'RM cost/kg prod', '% RM cost/kg prod']
        
        fd = self._df_combine()
        
        # Display the DataFrames for different permutations of kwargs
        # The fillna removes NaN from the displayed tables in Notebook, but
        # may goof up printing
        if decimals:
            if style == 'full':
                disp(fd.round(decimals).fillna(fill))
            elif style == 'compact':
                disp(fd[comp_col].round(decimals).fillna(fill))
        else:
            if style == 'full':
                disp(fd.fillna(fill))
            elif style == 'compact':
                disp(fd[comp_col].fillna(fill))

    def _df_combine(self, ):
        '''Combine the DataFrames for saving/exporting.'''
        # Copy the original DFs, and remove indexing
        fd = self.fulldata.reset_index()
        pmi = self.pmi.reset_index()
        
        # Combine the DFs. Set the index and then sort. Undo the multiindex
        # So compound names can be fixed
        concated = pd.concat([fd, pmi], sort=False)\
                    .set_index(['Prod', 'Compound']).sort_index()\
                    .reset_index()
        
        # Fix the compound names
        concated['Compound'] = concated['Compound'].str.replace(self._pre, '*')
        
        # Reset the index and return the DF
        return concated.set_index(['Prod', 'Compound'])                
            

class ColabCost(GenericCost):
    '''Costing class designed for the Colab Python environment.

    Parameters
    ----------
    materials_key : str
        The Google Sheet key to a materials list. This value can be found
        in the URL of the sheet.

    rxn_key : str
        The Google Sheet key to the reaction list. This value can be found
        in the URL of the sheet.

    final_prod : str
        The name of the overall final product of this route. This should
        be the same as a product from the reaction sheet.

    materials_sheet : int, str, optional (default = 0)
        The sheet to pull out of the materials Google Sheet. The default
        (0) is the first sheet in the spreadsheet. You could use a
        different number or a name if you want a different sheet from the
        spreadsheet.

    rxn_sheet : int, str, optional (default = 0)
        See `materials_sheet` description, except this is for the reaction
        Google Sheet.

    alt_mat_key : str, optional (default = None)
        A Google Sheet key for an optional, secondary materials sheet.
        This is useful if you have separate master and user materials
        sheets, for example.

    alt_mat_sheet : int, str, optional (default = 0)
        The sheet number/name for the secondary materials sheet. See
        `materials_sheet` description. 

    Attributes
    ----------
    final_prod : str
        The name of the overall final product of this route.
        
    rxns : DataFrame
        A DataFrame describing the original reactions in the given route. This 
        is idential to the reactions Google Sheet, and is not changed by any 
        of the helper functions in this class.
        
    materials : DataFrame
        A DataFrame describing all of the known materials. This will contain
        all of the materials from both of the given materials Google Sheets.
        
    cost : Numpy Float64
        The final cost of the described route. This value *will* include OPEX
        for the final reaction, if that value is given. 
        
    fulldata : DataFrame
        A DataFrame containing all the costing related values for the given 
        route. If the cost for the route has been calculated and an OPEX was 
        given, the product "Cost" values *will* include this additional OPEX 
        cost. The "RM cost/kg rxn" value will be the cost without the OPEX.

    Notes
    -----
    If there is a missing material or reaction, you'll get a printed
    error. Materials that are marked as being cost calculated will have
    their costs deleted, so they will need reactions defined in order to
    reset their costs.
    '''
    def __init__(self, materials_key, rxn_key, final_prod, materials_sheet=0,
            rxn_sheet=0, alt_mat_key=None, alt_mat_sheet=0):

        # Do some imports that are only possible in the Colab environment
        # This should prevent these from running in a non-Colab environment
        from oauth2client.client import GoogleCredentials
        from google.colab import auth
        from google.colab import files
        from IPython.display import display as disp
        # These will have to be made global
        global GoogleCredentials
        global auth
        global files
        global disp

        # Fix the final product and setup a mod variable
        super(ColabCost, self).__init__(final_prod)

        # Authenticate the Colab environment 
        auth.authenticate_user()
        self._gc = gspread.authorize(GoogleCredentials.get_application_default())
        
        # Set up the reaction DataFrame
        self._rxn_key = rxn_key
        self._rxn_sheet = rxn_sheet
        self._rxn_read()
        
        # Create the Materials DataFrame from a main sheet and an optional
        # alternate sheet.
        self._materials_key = materials_key
        self._materials_sheet = materials_sheet
        self._alt_mat_key = alt_mat_key
        self._alt_mat_sheet = alt_mat_sheet
        self._materials_build()                

        # Combine the reaction/materials sheets, add the new columns
        self.rxn_data_setup()
        
    def _materials_build(self, ):
        '''This function creates and combines the main and an optional
        alternate materials sheets.'''
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
        '''Read a Google sheet that defines the materials used in costing.
        '''
        mats = self._get_gsheet_vals(mat_key, wsheet)
        
        # Convert numeric/date columns
        mats['MW'] = pd.to_numeric(mats['MW'])
        mats['Density'] = pd.to_numeric(mats['Density'])
        mats['Cost'] = pd.to_numeric(mats['Cost'])
        #mats['Date'] = pd.to_datetime(mats['Date'])
        
        return mats
        
    def _rxn_read(self, ):
        '''Read a Google Sheet of reaction info.
        '''
        rxns = self._get_gsheet_vals(self._rxn_key,
                                     self._rxn_sheet)
        
        rxns['Equiv'] = pd.to_numeric(rxns['Equiv'])
        rxns['Volumes'] = pd.to_numeric(rxns['Volumes'])
        rxns['Sol Recyc'] = pd.to_numeric(rxns['Sol Recyc'])
        rxns['OPEX'] = pd.to_numeric(rxns['OPEX'])
        
        self.rxns = rxns
        
    def _get_gsheet_vals(self, key, sheet):
        '''General code for getting Google Sheet values and returning a 
        DataFrame.
        '''
        # Grab the Google sheet handle, pull down all values and make a 
        # DataFrame
        gsh = self._gc.open_by_key(key)
        # Differentiate between string and integer worksheets 
        if isinstance(sheet, str):
            ws = gsh.worksheet(sheet)
        else:
            ws = gsh.get_worksheet(sheet)
        vals = ws.get_all_values()
        val_df = pd.DataFrame(data=vals[1:], columns=vals[0])
        
        # Convert empty cells to NaN
        mask = (val_df == '')
        val_df[mask] = np.nan
        # Drop empty rows
        val_df.dropna(how='all', inplace=True)
        
        return val_df

    def excel_download(self, fname, decimals=2):
        '''Download the costing DataFrame as an Excel file.

        Note
        ----
        In some cases, this function will throw an error. In that case, try
        running this again in order to get it to work. 
        '''
        # Can set some keyword arguments here
        kwargs = {}
        # If decimals is given, set that value to the rounding for float
        # formatting in the output
        if decimals:
            kwargs['float_format'] = '%.{:d}f'.format(decimals)
            
        fd = self._df_combine()
        
        # Create the excel file
        with pd.ExcelWriter(fname) as writer:
            fd.to_excel(writer, sheet_name='Route Cost', **kwargs)
            
        # There seems to be a bit of a lag before you can download
        # the file, this delay might fix some of the errors this causes
        time.sleep(2)
        files.download(fname)
        

