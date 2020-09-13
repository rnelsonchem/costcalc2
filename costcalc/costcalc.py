'''
Chemical Reaction Cost Cacluation Routines.
Adapted from the Excel spreadsheets prepared by Saeed Ahmad, PhD.
(C) Ryan Nelson
'''
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Setting the display to print function
# This gets changed for Jupyter Notebook/IPython sessions
try:
    from IPython.display import display as disp
    from IPython.display import Javascript
except:
    disp = print

# Set up some plotting stuff for the notebooks
plt.style.use('ggplot')
plt.rc('figure', dpi=150)

# For future reference -- Can't use this in conjunction with "fillna" becuase
# the DF doesn't display correctly 
# Set Pandas precision
#pd.set_option('precision', 2)

# Set pandas to display lots of DataFrame rows so things don't get cut out
pd.options.display.max_rows = 1000


class ExcelCost(object):
    '''Costing class designed for local Excel/csv spreadsheets.

    This can also act as the base class for other subclasses.

    Parameters
    ----------
    materials_file : str
        The name/path of the materials list file. Can be xlsx or csv.

    rxn_file : str
        The name/path of the file defining the reactions. can be xlsx or csv.

    alt_mat_file : str, optional (default = None)
        The name/path of an optional, secondary materials sheet.  sheet. This
        is useful if you have separate master and user materials sheets, for
        example. Can be an xlsx or csv.

    final_prod : str
        Defines the final product name for costing calculations. This should
        be the same name as in the material/reaction sheet.

    materials_sheet : int, str, optional (default = 0)
        The sheet to pull out of the materials spreadsheet. The default
        (0) is the first sheet in the spreadsheet. You could use a
        different number or a name if you want a different sheet from the
        spreadsheet.

    rxn_sheet : int, str, optional (default = 0)
        See `materials_sheet` description, except this is for the reaction
        Google Sheet.

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

    pmi : DataFrame
        A DataGrame containing the PMI for each reaction and the overall
        route. There will be an extra column with a bunch of weird names. This
        is necessary for sorting and can be ignored.

    Notes
    -----
    If there is a missing material or reaction, you'll get a printed
    error. Materials that are marked as being cost calculated will have
    their costs deleted, so they will need reactions defined in order to
    reset their costs.
    '''
    def __init__(self, materials_file, rxn_file, final_prod,
            materials_sheet=0, rxn_sheet=0, alt_mat_file=None,
            alt_mat_sheet=0):
        # We need to store this for costing
        self.final_prod = final_prod        

        # Set up the reaction DataFrame
        self._rxn_file = rxn_file
        self._rxn_sheet = rxn_sheet
        self._rxn_read()

        # Create the Materials DataFrame from a main sheet and an optional
        # alternate sheet.
        self._materials_file = materials_file
        self._materials_sheet = materials_sheet
        self._alt_mat_file = alt_mat_file
        self._alt_mat_sheet = alt_mat_sheet
        self._materials_build()                

        # Combine the reaction/materials sheets, add the new columns
        self.rxn_data_setup()

        # Look for common errors in the input.
        self._sanity_check()

    def _rxn_read(self, ):
        '''Read an Excel sheet that defines the reactions.
        '''
        self.rxns = self._excel_csv_reader(self._rxn_file, self._rxn_sheet,
                            dtypes={'Step':str, 'Cost calc':str})

    def _excel_csv_reader(self, fname, fsheet, dtypes=None):
        '''A simple Excel/CSV reader function for both reaction and materials
        files.
        '''
        # Read the file, drop NaN-only and commented rows.
        if fname[-4:].lower() == 'xlsx':
            df = pd.read_excel(fname, fsheet, dtype=dtypes, comment='#')\
                            .dropna(how='all')
        elif fname[-3:].lower() == 'csv':
            df = pd.read_csv(fname, dtype=dtypes, comment='#')\
                            .dropna(how='all')
        return df
        
    def _materials_build(self, ):
        '''Read and combine the main and an optional alternate materials
        sheets.
        '''
        materials = self._materials_read(self._materials_file,
                self._materials_sheet)

        # If an alternative materials key is given, combine that materials
        # sheet with the main one
        if self._alt_mat_file:
            alt_mats = self._materials_read(self._alt_mat_file,
                    self._alt_mat_sheet)
            # Concatenate the sheets. Reset the index so that it is
            # consecutively numbered
            materials = pd.concat([materials, alt_mats], sort=False)\
                    .reset_index(drop=True)

        # Set the final materials sheet
        self.materials = materials

    def _materials_read(self, mat_file, wsheet):
        '''Read an Excel file sheet that defines the materials used in
        costing.
        '''
        mats = self._excel_csv_reader(mat_file, wsheet,)
        return mats

    def rxn_data_setup(self, ):
        '''Setup the full data set for the upcoming cost calculations. 
        
        This method merges the materials and reaction sheets into a combined
        DataFrame called `fulldata`. It also serves as a "reset" of sorts
        after a costing calculation, if you want to start over with a
        different costing model.
        '''
        # Merge the materials and reaction DataFrames. A few columns are
        # dropped, which are not necessary for calculations. The merge happens
        # on the rxns DataFrame ('right'), which means that missing materials
        # will be fairly obvious (no MW, e.g.).
        mat_keeps = ['Compound', 'MW', 'Density', 'Cost']
        rxn_keeps = ['Step', 'Compound', 'Equiv', 'Volumes', 'Relative', 
                     'Sol Recyc', 'Cost calc', 'OPEX',]
        fulldata = pd.merge(self.materials[mat_keeps], self.rxns[rxn_keeps],
                            on='Compound', how='right')

        # Find the step number for the final product. 
        fp_mask = fulldata.Compound == self.final_prod
        self._fp_idx = fulldata.loc[fp_mask, 'Cost calc'].iloc[0]

        # Set MultiIndex
        fulldata.set_index(['Step', 'Compound'], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not
        # copies
        fulldata = fulldata.sort_index()
        
        # Save the full data set
        self.fulldata = fulldata

        # Add a modified variable placeholder. This will store modified
        # values for later processing
        self._mod_vals = []
        # Add the empty columns
        self._column_clear()

    def _sanity_check(self, ):
        '''Run some sanity checks on the DataFrames to catch common errors.

        Some of the errors are tricky, and the error output is unhelpful.
        These are some checks that will throw some "sane" errors if things are
        not correct.
        ''' 
        # Check for a missing material -- Everything should have a MW
        mw_mask = self.fulldata['MW'].isna()
        if mw_mask.any():
            print('You are missing a MW!!!')
            print('May be a mismatch between the reaction and materials file.')
            print('Material might be missing from materials sheet.')
            print('Check these materials.')
            disp(self.fulldata.loc[mw_mask, 'MW'])
            raise ValueError('Yikes! Read the note above.')
            
        # Check for a missing cost, which is not being calculated
        # This is a tricky error because the costing will run just fine with 
        # NaNs. Check for this by looking for NaN in both Cost and Calc columns
        cost_mask = (self.fulldata['Cost calc'].isna() & \
                    self.fulldata['Cost'].isna())
        if cost_mask.any():
            print('You are missing a necessary material cost!!')
            print('You may need to indicate a Step in the "Cost calc" column.')
            print('Check these columns.')
            disp(self.fulldata.loc[cost_mask, ['Cost', 'Cost calc']])
            raise ValueError('Yikes! Read the note above.')
        
        # Check for duplicated materials. This will probably be a big issue
        # with two materials sheets.
        dup_cpds = self.materials['Compound'].duplicated()
        if dup_cpds.any():
            print('You have a duplicate material!!!')
            print("These compounds are duplicated in your materials sheet.")
            disp(self.materials.loc[dup_cpds, 'Compound'])
            raise ValueError('Yikes! Read the note above.')

        # Check for duplicated materials in a single reaction.
        # When you select a single value from a reaction, you'll get a series
        # and not a float, e.g.
        dup_rxn = self.fulldata.loc[(self._fp_idx, self.final_prod), 'MW']
        if isinstance(dup_rxn, pd.Series):
            print('You have a duplicated material in a single reaction.')
            print('Check these lines:')
            gb = self.fulldata.groupby(['Prod', 'Compound'])
            for prod, group in gb:
                if group.shape[0] > 1:
                    disp(prod)
            raise ValueError('Yikes! Read the note above.')
            
        # Check that all the solvent information is given
        sol_mask = ~self.fulldata['Volumes'].isna() 
        sol_chk = self.fulldata.loc[sol_mask, 'Relative'].isna() | \
                self.fulldata.loc[sol_mask, 'Density'].isna() | \
                self.fulldata.loc[sol_mask, 'Sol Recyc'].isna()
        # If anything is missing, print a note
        if sol_chk.any():
            sol_cols = ['Density', 'Volumes', 'Relative', 'Sol Recyc']
            print('You are missing some solvent information.')
            print('Check the following lines.')
            sol_mask2 = sol_mask & sol_chk
            disp(self.fulldata.loc[sol_mask2, sol_cols])
            raise ValueError('Yikes! Read the note above.')
        
        # Check to make sure that the "Relative" compound for a solvent
        # is acutally contained in the step
        sol_mask = ~self.fulldata['Relative'].isna()
        for cpd in self.fulldata[sol_mask].itertuples():
            new_idx = (cpd.Index[0], cpd.Relative)
            if new_idx not in self.fulldata.index:
                print('One of your "Relative" compounds is not correct.')
                print(f'"{cpd.Relative}" is not in Step {cpd.Index[0]}.')
                raise ValueError('Yikes! Read the note above.')
        
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

        step : None, int, optional (Default = None)
            The reaction step number for which this value will be changed. If
            this is `None` (default), then all the values for the given
            compound (`cpd`) will be set to the same value. This is mostly
            important for something like `val_type`='Equiv'. Clearly, you
            would only want to change the number of equivalents for a specific
            reaction. If this parameter is left as `None`, the equivalents for
            a given compound in all reactions will be set to the same value.
        
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
            # The step needs to be a string, not in int as might be expected
            cells = (str(step), cpd)
        
        # Does this position actually exist???
        try:
            self.fulldata.loc[cells, val_type]
        except KeyError:
            # If not, remove the values from the list and throw an error
            self._mod_vals.pop()
            if not step:
                err = '"' + cpd + '"'
            else:
                err = 'Step "' + str(step) + '", "' + cpd + '"'
            print('Oops! ' + err + " doesn't exist. Check") 
            print("your reactions to find the correct Step/Compound.")
            raise 
            
        self.fulldata.loc[cells, val_type] = val
        # The "Cost calc" flag must be set to np.nan when setting a cost. 
        # This is necessary for % RM cost calcs, e.g.
        if val_type == 'Cost':
            self.fulldata.loc[cells, 'Cost calc'] = np.nan

    def value_scan(self, cpd, start, stop, npts, val_type='Cost', step=None):
        '''Scan a range of values for a given material.
        
        Parameters
        ----------
        See `value_mod` method description, except for the following.

        start : int, float
            The starting value for the scan.

        stop : int, float
            The ending value for the scan.

        npts : int
            The numbers of points to calculate the costs between `start` and
            `stop`

        Returns
        -------
        Pandas DataFrame
            This DataFrame has 'Values' and 'Costs' columns for the input
            values and output costs, repectively. 

        Notes
        -----
        Although this method recalculates the cost for every value, it does
        not modify the original `fulldata` or `cost` attributes. 
        '''
        # Create the values array
        vals = np.linspace(start, stop, npts)
       
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
        self.cost = self.fulldata.loc[(self._fp_idx, self.final_prod), 
                                  'RM cost/kg rxn']
        
        return pd.DataFrame({'Values':vals, 'Costs':all_costs})

    def plot_scan(self, cpd, start, stop, npts, val_type='Cost', step=None, 
                legend=None):
        '''Plot a range of values.

        Parameters
        ----------
        See `value_scan` method description, except for the following.

        legend : bool, str
            This will add a legend to the plot. If you use the value `True`,
            the compound name will be added to the legend. Otherwise, you can
            pass a custom string if you want that to be in the legend instead.

        '''
        # Calculate all the costs for the given values
        costs = self.value_scan(cpd, start, stop, npts, val_type=val_type,
                            step=step)
        
        # If legend is selected, make sure the label is set properly
        if legend == True:
            label = cpd
        else:
            label = legend

        # Plot the values
        plt.plot(costs['Values'], costs['Costs'], 'o', label=label)
        # Generate the legend, if requested
        if legend:
            plt.legend()


    def swap(self, cpd_old, cpd_new, step=None):
        '''Swap one compound for another in the route.

        Parameters
        ----------
        cpd_old : str
            This is the name of the compound that you want to remove from the
            route costing.

        cpd_new : str or Cost class
            This can either be the name of the new compound to add into the
            route or a Cost class instance, such as ExcelCost. If it is just a
            name, the material properties will be pulled out of the materials
            database. If the Cost class instance is used, the the final
            product information will be pulled from there. (In that case,
            then, the `cost_calc` method must have been run on that instance
            to define the final cost of that material.

        step : int or None, optional (default = None)
            This is the step number for which to do the swap. The default
            (None) is to swap all instances of the particular compound.

        Note
        ----
            If you swap out a compound that was marked to have its cost
            calculated, this method will switch off this feature, which is
            most likely what was intended. 
        '''
        # If the new compound is a string, look in the materials database
        if isinstance(cpd_new, str):
            mat_mask = self.materials.Compound == cpd_new
            # Throw an error if it is not there
            if not mat_mask.any():
                raise ValueError("Oops! Your new compound isn't in the"
                                " materials list.")
            # There should only be one entry, so we'll select that one
            mat_vals = self.materials[mat_mask].iloc[0]
            # Set some values that will be used for updating
            cpd_name = cpd_new
            mw = mat_vals['MW']
            density = mat_vals['Density']
            cost = mat_vals['Cost']
        
        # If the new compound is a costing instance...
        elif isinstance(cpd_new, (ExcelCost, ColabCost)):
            # The value selection is a little different
            cpd_name = cpd_new.final_prod
            cpd_loc = (cpd_new._fp_idx, cpd_new.final_prod)
            mw = cpd_new.fulldata.loc[cpd_loc, 'MW']
            density = cpd_new.fulldata.loc[cpd_loc, 'Density']
            cost = cpd_new.fulldata.loc[cpd_loc, 'Cost']
        
        # Else you've added the wrong kind of new compound
        else:
            raise ValueError("Oops!! Your new compound is not "
                            "of the correct type.")

        # Process the fulldata array
        # First, remove the MultiIndex
        fd_rst = self.fulldata.reset_index()
        # Check for the requested compound
        cpd_mask = fd_rst['Compound'] == cpd_old
        # If a particular step is chosen, select only that step
        if step:
            cpd_mask = cpd_mask & (fd_rst['Step'] == str(step))

        if not cpd_mask.any():
            raise ValueError("Oops! '" + cpd_old + "' isn't in your "
                            "current route.")
        # Swap out the compound names
        fd_rst.loc[cpd_mask, 'Compound'] = cpd_name
        # Reset index and fulldata attribute
        self.fulldata = fd_rst.set_index(['Step', 'Compound'])

        # Set the values using `value_mod`
        self.value_mod(cpd_name, mw, val_type='MW', step=step)
        self.value_mod(cpd_name, density, val_type='Density', step=step)
        self.value_mod(cpd_name, cost, step=step)

    def calc_cost(self, ):
        '''Calculate the cost of the route. 
        '''
        # Save a time stamp so it can be displayed later
        self._now = pd.Timestamp.now('US/Eastern').strftime('%Y-%m-%d %H:%M')
        # Prep the DataFrame
        self._column_clear()
        # Run the costing and set the cost attribute
        self.cost = self._rxn_cost(self.final_prod, self._fp_idx)
        # Post process the DataFrame
        self._rxn_data_post()
        
    def _rxn_cost(self, prod, step, amp=1.0):
        '''The recursive cost calculating function. 
        
        This is the workhorse function of the whole process, but is not meant
        to be called on its own. Typically, you'll want to call `calc_cost` to
        get the costing information.
        
        Parameters
        ----------
        prod : str
            The name of the reaction to cost. This should also be the name of
            the final product for that reaction.

        step : str
            The step number for which to calculate the cost. Even though the
            step numbers are numbers, this needs to be given as a string.

        amp : float, optional (default = 1.0)
            This number is an amplifier that increases some of the values,
            e.g. masses of materials, based on how much material is being used
            for the overall final product.   
        '''
        # Select out the reaction of interest from the full data set. Saves
        # some typing. (Make this a copy to enusre it isn't given as a view.)
        data = self.fulldata.loc[step].copy()

        # Kg of nonsolvent materials used per equivalent
        data['kg/kg rxn'] = data['Equiv']*data['MW']

        # Amount of solvent
        # First figure out which materials are solvents 
        mask = ~data['Volumes'].isna()
        sols = data[mask].copy()
        if sols.size != 0:
            # Find the material that the volumes are relative to (taking only
            # the first one... This might not be great.
            cpd_rel = sols['Relative'][0]
            # What is the kg of relative cpd?
            amt_rel = data.loc[cpd_rel, 'kg/kg rxn']
            # Calculate the mass of solvent. Take into account the solvent
            # recycyling 
            # kg sol = Volume*Density*(1-Recycle)*(kg SM)
            sols['kg/kg rxn'] = sols['Volumes']*sols['Density']*\
                    (1 - sols['Sol Recyc'])*amt_rel
            data[mask] = sols

        # Normalize the kg of reaction
        data['kg/kg rxn'] /= data.loc[prod, 'kg/kg rxn']

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
            # The new step is needed as well
            new_stp = data.loc[cpd, 'Cost calc']
            # Run the cost calculation for the unknown compound
            cst = self._rxn_cost(cpd, new_stp, amp*new_amp)
            # Set the calculated cost
            data.loc[cpd, 'Cost'] = cst

        # Calculate the cost for each material in the reaction
        data['RM cost/kg rxn'] = data['kg/kg rxn']*data['Cost']
        # The product cost will be the sum of all the reactant/solvent costs
        data.loc[prod, 'RM cost/kg rxn'] = data['RM cost/kg rxn'].sum()
        # Set the "Cost" to the calculated value
        self.fulldata.loc[(step, prod), 'Cost'] = \
                data.loc[prod, 'RM cost/kg rxn']

        # Calculate % costs for individual rxn
        # = (RM cost/kg rxn)/(RM cost/kg rxn for the rxn product)
        p_rm_cost = data['RM cost/kg rxn']*100/data.loc[prod, 'RM cost/kg rxn']
        data['% RM cost/kg rxn'] = p_rm_cost.values
        # Remove the % cost for the rxn product
        data.loc[prod, '% RM cost/kg rxn'] = np.nan

        # These are the costs for ultimate product
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn 
        # cost. However, for feeder reactions this will get amplified by 
        # each step   
        data['RM cost/kg prod'] = data['RM cost/kg rxn']*amp
        
        # This sets the number of kg of each material per kilogram of product
        # This is done by multiplying the per reaction value by the amplifier
        # This isn't necessary for costing, but we can use it for PMI
        data['kg/kg prod'] = data['kg/kg rxn']*amp
        
        # Set the values in the big DataFrame with this slice. This goofy call
        # is necessary to make sure the data are set correctly. 
        self.fulldata.loc[(step, data.index), data.columns] = data.values
        
        # Return the calculated product cost, which is required for the 
        # recursive nature of the algorithm. In addition, an optional OPEX
        # may be added to take into acount production costs of the cpd
        if np.isnan(data.loc[prod, 'OPEX']):
            return data.loc[prod, 'RM cost/kg rxn']
        else:
            return data.loc[prod, 'RM cost/kg rxn'] + data.loc[prod, 'OPEX']
    
    def _rxn_data_post(self,):
        '''Calculate some values after the final cost of the reaction is
        determined. 

        This includes the final "% RM cost/kg prod", setting the final cost
        with an optional OPEX, and all PMI calculations. Some values are
        filtered out as well to make the column sums sensible.
        '''
        prod = self.final_prod
        step = self._fp_idx
        
        # If an OPEX for the final reaction is given, add that to the cost
        # of the final product
        opex = self.fulldata.loc[(step, prod), 'OPEX']
        if not np.isnan(opex):
            self.fulldata.loc[(step, prod), 'Cost'] = self.cost
                
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
        self.fulldata.loc[(step, prod), 'kg/kg prod'] = 1.

        # PMI Calculations
        # Need to append this prefix for sorting purposes
        # There will be a funny column in this DF with these values...
        self._pre = 'zzzz'
        
        # First of all, calculate the PMI for each reaction individually
        gb = self.fulldata[['kg/kg rxn']].groupby('Step')
        rxn_pmi = gb.sum().reset_index()
        rxn_pmi['Compound'] = self._pre + 'Step ' + rxn_pmi['Step'] + ' PMI'
        
        # The full route PMI is not the sum of the above, but is the sum of
        # the 'kg/kg prod' column. We need to make this into a DataFrame to
        # merge with the per reaction values above
        df_vals = {'kg/kg prod': [self.fulldata['kg/kg prod'].sum()], 
                   'Step': [self._fp_idx],
                   'Compound': [self._pre*2 + 'Full Route PMI']
                   }
        full_pmi = pd.DataFrame(df_vals)

        # Merge the per-reaction and full PMI
        self.pmi = pd.concat([rxn_pmi, full_pmi], 
                             sort=False).set_index('Step')

    def results(self, style='compact', decimals=2, fill='-'):
        '''Print the results of the costing calculation.

        Parameters
        ----------
        style : str, optional (Default = 'compact')
            This sets the style of the displayed costing DataFrame.
            `'compact'` prints a DataFrame that has been truncated slightly.
            `'full'` prints the entire DataFrame.

        decimals : int or None, optional (Default = 2)
            How many decimal places to show in the table. Set this to `None`
            if you want full precision.

        fill : str or None, optional ('-')
            Fill NaN values with this string. This makes the table a little
            easier to read. Set this to `None` if you want to see the table
            with the typical NaN labels.
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
        comp_col = ['Cost', 'Equiv',] 
        if self.fulldata.Volumes.any():
            comp_col.extend(['Volumes', 'Sol Recyc',])
        if self.fulldata.OPEX.any():
            comp_col.append('OPEX') 
        comp_col.extend(['kg/kg rxn', 'RM cost/kg rxn', '% RM cost/kg rxn',
                    'kg/kg prod', 'RM cost/kg prod', '% RM cost/kg prod'])
        
        # Combine the fulldata and pmi DataFrames
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
        '''Combine the fulldata and pmi DataFrames for saving/exporting.
        '''
        # Copy the original DFs, and remove indexing
        fd = self.fulldata.reset_index()
        pmi = self.pmi.reset_index()
        
        # Combine the DFs. Set the index and then sort. Undo the multiindex
        # So compound names can be fixed
        concated = pd.concat([fd, pmi], sort=False)\
                    .set_index(['Step', 'Compound']).sort_index()\
                    .reset_index()
        
        # Fix the compound names
        concated['Compound'] = concated['Compound'].str.replace(self._pre, '*')
        
        # Reset the index and return the DF
        return concated.set_index(['Step', 'Compound'])                

    def sensitivity(self, col='Equiv', frac=0.1, decimals=2):
        '''Do a sensitivity analysis for the equivalents of reagents.

        Parameters
        ----------
        col : str, optional (Default = 'Equiv')
            Which column from the `fulldata` DataFrame should be used for the
            sensitivity analysis. 

        frac : float, optional (Default = 0.1)
            Fractional percentage to increase/decrease the values by before
            recosting. The default is 0.1, which is the same as +/- 10%. 

        decimals : int or None, optional (Default = 2)
            How many decimal places to display. If `None`, the full precision
            DataFrame will be displayed.
        '''
        # Make a new DF for sensitivity analysis
        # Make values that are a certain percent above and below the current
        # numbers
        sens = self.fulldata[[col]].dropna()
        sens['Val low'] = sens[col]*(1 - frac)
        sens['Val high'] = sens[col]*(1 + frac)
        
        # Re-run the costing under the current conditions, which resets the
        # cost and fulldata variables. 
        self.calc_cost()
        # Make copies of these values so they don't change
        cost_save = self.cost
        fd_save = self.fulldata.copy()

        # Loop through the values and calculate the cost if these values
        # increase or decease by the percent given
        for step_cpd, vals in sens.iterrows():
            step, cpd = step_cpd
            # Low values
            self.value_mod(cpd, vals['Val low'], val_type=col, step=step)
            self.calc_cost()
            cost_low = self.cost
            cost_low_per = (cost_low*100./cost_save) - 100
            # Just to be safe - Drop the modified value
            self._mod_vals.pop()

            # High values
            self.value_mod(cpd, vals['Val high'], val_type=col, step=step)
            self.calc_cost()
            cost_high = self.cost
            cost_high_per = (cost_high*100./cost_save) - 100
            self._mod_vals.pop()

            # Set the values in the sensitivity DataFrame
            sens.loc[(step, cpd), 'Cost low'] = cost_low
            sens.loc[(step, cpd), 'Cost high'] = cost_high
            sens.loc[(step, cpd), '% low'] = cost_low_per
            sens.loc[(step, cpd), '% high'] = cost_high_per

            # Reset the original values.  This is necessary so the next step
            # in the loop starts from the original position.
            self.cost = cost_save
            self.fulldata = fd_save.copy()
        
        if decimals:
            return sens.round(decimals)
        else:
            return sens

    def excel_save(self, fname, decimals=2):
        '''Save the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str
            The name you want to give to the Excel file.

        decimals : str or None, optional (default = 2)
            The number of decimal places to display in the Excel sheet. If
            `None`, then the full precision will be saved. 

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
        
        # Create the excel file. Can only save with the date and not the time
        with pd.ExcelWriter(fname) as writer:
            fd.to_excel(writer, 
                        sheet_name='As of ' + self._now.split()[0], 
                        **kwargs)
            

class ColabCost(ExcelCost):
    '''Costing class designed for the Colab Python environment.

    Parameters
    ----------
    materials_key : str
        The Google Sheet key to a materials list. This value can be found
        in the URL of the sheet.

    rxn_key : str
        The Google Sheet key to the reaction list. This value can be found
        in the URL of the sheet.

    alt_mat_key : str, optional (default = None)
        A Google Sheet key for an optional, secondary materials sheet.
        This is useful if you have separate master and user materials
        sheets, for example.

    '''
    def __init__(self, materials_key, rxn_key, final_prod, materials_sheet=0,
            rxn_sheet=0, alt_mat_key=None, alt_mat_sheet=0):
        # Do some imports that are only possible in the Colab environment
        # This should prevent these from running in a non-Colab environment
        from oauth2client.client import GoogleCredentials
        from google.colab import auth
        from google.colab import files
        import gspread
        # These will have to be made global
        global GoogleCredentials
        global auth
        global files
        global gspread

        # Authenticate the Colab environment 
        auth.authenticate_user()
        self._gc = gspread.authorize(GoogleCredentials.get_application_default())
        
        # Fix the final product and setup a mod variable
        super(ColabCost, self).__init__(materials_key, rxn_key, final_prod,
                materials_sheet, rxn_sheet, alt_mat_key, alt_mat_sheet)
        
    def _materials_read(self, mat_key, wsheet):
        '''Read a Google sheet that defines the materials used in costing.

        Parameters
        ----------
        mat_key : str
            The unique key for a Google spreadsheet that defines the
            materials. 

        wsheet : str or int
            The specific sheet to extract from the Google spreadsheet. 
        '''
        mats = self._get_gsheet_vals(mat_key, wsheet)

        # Convert numeric/date columns. Everything is read from a Google sheet
        # as strings
        num_cols = ['MW', 'Density', 'Cost'] 
        for nc in num_cols:
            mats[nc] = pd.to_numeric(mats[nc])
        #mats['Date'] = pd.to_datetime(mats['Date'])

        return mats
        
    def _rxn_read(self, ):
        '''Read a Google Sheet of reaction info.
        '''
        rxns = self._get_gsheet_vals(self._rxn_file,
                                     self._rxn_sheet)

        # Set some rxns columns to numeric values. Everything is read from a
        # Google sheet as strings
        num_cols = ['Equiv', 'Volumes', 'Sol Recyc', 'OPEX']
        for nc in num_cols:
            rxns[nc] = pd.to_numeric(rxns[nc])
        
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

        # Drop lines with comment markers
        mask = ~val_df.iloc[:, 0].str.startswith('#')
        val_df = val_df[mask]
        
        return val_df

    def results(self, style='compact', decimals=2, fill='-'):
        '''Print the results of the costing calculation.

        Parameters
        ----------
        style : str, optional (Default = 'compact')
            This sets the style of the displayed costing DataFrame.
            `'compact'` prints a DataFrame that has been truncated slightly.
            `'full'` prints the entire DataFrame.

        decimals : int or None, optional (Default = 2)
            How many decimal places to show in the table. Set this to `None`
            if you want full precision.

        fill : str or None, optional ('-')
            Fill NaN values with this string. This makes the table a little
            easier to read. Set this to `None` if you want to see the table
            with the typical NaN labels.
        '''
        # This makes the max Colab output window very large, so that
        # DataFrames are not put into separate scroll windows, which is very
        # annoying for users. Unfortunately, I need to copy/paste the doc
        # string for results method, though...
        # See: https://github.com/googlecolab/colabtools/issues/541
        iframe_h = 'google.colab.output.setIframeHeight(0, true, {\n'\
                   '  maxHeight: 5000,\n'\
                   '})'
        disp(Javascript(iframe_h)) 
        # Call results method from ExcelCost
        super(ColabCost, self).results(style=style, decimals=decimals,
                                       fill=fill)


    def excel_save(self, fname, decimals=2):
        '''Download the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str
            The name you want to give to the Excel file.

        decimals : str or None, optional (default = 2)
            The number of decimal places to display in the Excel sheet. If
            `None`, then the full precision will be saved. 

        Note
        ----
        In some cases, this function will throw an error. In that case, try
        running this again in order to get it to work. 
        '''
        super(ColabCost, self).excel_save(fname, decimals)

        # There seems to be a bit of a lag before you can download
        # the file, this delay might fix some of the errors this causes
        time.sleep(2)
        files.download(fname)
        

