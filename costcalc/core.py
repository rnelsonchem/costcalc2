'''
Chemical Reaction Cost Calculation Routine.
Adapted from the Excel spreadsheets prepared by Saeed Ahmad, PhD.
(C) Ryan Nelson
'''
import time

import numpy as np
import pandas as pd

# Column dictionary
# This dictionary will be used for creating dynamic excel sheets, for which
# I'll need to know the correct column labels
ecols = {'Step':'A', 'Compound':'B', 'MW':'C', 'Density':'D', 'Cost':'E',
        'Equiv':'F', 'Volumes':'G', 'Relative':'H', 'Sol Recyc':'I',
        'Cost Calc':'J', 'OPEX':'K', 'kg/kg rxn':'L', 'RM cost/kg rxn':'M',
        '% RM cost/kg rxn':'N', 'kg/kg prod':'O', 'RM cost/kg prod':'P',
        '% RM cost/kg prod':'Q',
        }

class CoreCost(object):
    def __init__(self, rxns, materials, final_prod):
        # We need to store the input values/DataFrames
        self.rxns = rxns
        self.materials = materials
        self.final_prod = final_prod        

        # Combine the reaction/materials sheets, add the new columns
        self.rxn_data_setup()

        # Look for common errors in the input.
        self._sanity_check()

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
        mat_keeps = ['Compound', 'MW', 'Density', 'Cost', 'Notes']

        rxn_keeps = ['Step', 'Compound', 'Equiv', 'Volumes', 'Relative',
                    'Sol Recyc', 'Cost calc', 'OPEX', 'Notes']
        # Check if an "Amount" column is present. This is used to calculate
        # equivalents. If this is present, add this column to our "keep" list
        amt = False
        if 'Amount' in self.rxns.columns:
            amt = True 
            rxn_keeps = rxn_keeps[:2] + ['Amount'] + rxn_keeps[2:]

        fulldata = pd.merge(self.materials[mat_keeps], self.rxns[rxn_keeps],
                            on='Compound', how='right')\
                            .rename({'Notes_x': 'Material Notes',
                                'Notes_y': 'Reaction Notes'}, axis=1)

        # Find the step number for the final product. 
        fp_mask = fulldata.Compound == self.final_prod
        self._fp_idx = fulldata.loc[fp_mask, 'Cost calc'].iloc[0]

        # If the "Amount" column is present, then you will need calculate the
        # "Equiv" based on the Amount given
        if amt:
            # Group everything by "Step" column
            grp = fulldata.groupby('Step')
            for idx, sb_grp in grp:
                # If there are no amounts for a given Step, then skip
                if ~sb_grp['Amount'].any():
                    continue
                # Mask out only the values that have amounts. This is
                # important in case a mixture of mass and equiv is used.
                # Also check to make sure this isn't a solvent
                mask = ~sb_grp['Amount'].isna() & sb_grp['Volumes'].isna()
                # Calculate the Equivalents
                sb_grp.loc[mask, 'Equiv'] = sb_grp.loc[mask, 'MW']*\
                                        sb_grp.loc[mask, 'Amount']
                # Assume the first entry is the limiting reagent. Normalize
                # the rest of the equivalents to that value
                sb_grp.loc[mask, 'Equiv'] = sb_grp.loc[mask, 'Equiv']/\
                                       (sb_grp.iloc[0]['Equiv'])
                # Set the values in the full DataFrame
                mask = fulldata['Step'] == idx
                fulldata.loc[mask, 'Equiv'] = sb_grp['Equiv']
            # Remove the Amount column
            fulldata.drop('Amount', axis=1, inplace=True)

        # Add an indexing column, which will be used for sorting purposes
        fulldata['idx_col'] = np.arange(fulldata.shape[0])
        # Set MultiIndex
        fulldata.set_index(['idx_col', 'Step', 'Compound'], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not
        # copies
        fulldata = fulldata.sort_index()
        # Remove sorting index, drop column
        fulldata = fulldata.reset_index('idx_col').drop('idx_col', axis=1)
        
        # Save the full data set
        self.fulldata = fulldata

        # Add a modified variable placeholder. This will store modified
        # values for later processing
        self._mod_vals = []
        # Add the empty columns
        self._column_clear()

    def _column_clear(self, excel=False):
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
        # For dynamic Excel
        if excel:
            excel_cols = ['Cost dyn', 'kg/kg rxn dyn', 'RM cost/kg rxn dyn', 
                        '% RM cost/kg rxn dyn', 'kg/kg prod dyn', 
                        'RM cost/kg prod dyn', '% RM cost/kg prod dyn',
                        ]
            for col in excel_cols:
                self.fulldata[col] = ''

    def _sanity_check(self, ):
        '''Run some sanity checks on the DataFrames to catch common errors.

        Some of the errors are tricky, and the error output is unhelpful.
        These are some checks that will throw some "sane" errors if things are
        not correct.
        ''' 
        # Check for a missing material -- Everything should have a MW
        # If not, this may suggest that there is a line that should not be
        # empty
        mw_mask = self.fulldata['MW'].isna()
        if mw_mask.any():
            err_line = 'Missing MW or bad line in input file! \n'\
                'You are missing one or more MW values or else there is \n'\
                'a cell entry in a "blank" line. Make sure that all empty\n'\
                'reaction/materials lines are completely blank, i.e. \n'\
                'entry free (even spaces!).'
            print(err_line)
            disp(self.fulldata[mw_mask])
            raise ValueError(err_line)

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
        
    def calc_cost(self, excel=False):
        '''Calculate the cost of the route. 
        '''
        # Save a time stamp so it can be displayed later
        self._now = pd.Timestamp.now('US/Eastern').strftime('%Y-%m-%d %H:%M')
        # Prep the DataFrame
        self._column_clear(excel=excel)
        # Set a column of unique row labels
        nrows = self.fulldata.shape[0]
        nrowdig = len(str(nrows))
        if excel:
            self.fulldata['rnum'] = [f'r{i:0{nrowdig:d}d}' for i in range(nrows)]
        # Run the costing and set the cost attribute
        self.cost = self._rxn_cost(self.final_prod, self._fp_idx, excel=excel)
        # Post process the DataFrame
        self._rxn_data_post(excel=excel)
        
    def _rxn_cost(self, prod, step, amp=1.0, eamp='', excel=False):
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
        # And for Excel
        if excel:
            data['kg/kg rxn dyn'] = '=' + ecols['Equiv'] + data['rnum'] + '*'\
                    + ecols['MW'] + data['rnum']
        # For Excel, normalize the data here. This will get overwritten for
        # solvents. This is a little weird for the product because it will be
        # the same in numerator and denominator, but that will keep things
        # fully interactive.
        if excel:
            data['kg/kg rxn dyn'] += '/(' + ecols['Equiv'] +\
                    data.loc[prod,'rnum'] + '*' + ecols['MW'] +\
                    data.loc[prod, 'rnum'] + ')'

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
            if excel:
                amt_rel_e = data.loc[cpd_rel, 'rnum']
            # Calculate the mass of solvent. Take into account the solvent
            # recycyling 
            # kg sol = Volume*Density*(1-Recycle)*(kg SM)
            sols['kg/kg rxn'] = sols['Volumes']*sols['Density']*\
                    (1 - sols['Sol Recyc'])*amt_rel
            # And for Excel
            if excel:
                sols['kg/kg rxn dyn'] = '=' + ecols['Volumes'] +\
                        sols['rnum'] + '*' + ecols['Density'] +\
                        sols['rnum'] + '*' + '(1 - ' + ecols['Sol Recyc']\
                        + sols['rnum'] + ')*' + ecols['kg/kg rxn'] +\
                        amt_rel_e
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
            # And for dynamic excel output
            if excel:
                new_amp_enum = data.loc[cpd, 'rnum']
                new_amp_eden = data.loc[prod, 'rnum']
                kg_col = ecols['kg/kg rxn']
            # This will be '*(C#/D#)'. This is different for Excel because we
            # can't normalize the kg/kg rxn column
                new_eamp = f'*({kg_col}{new_amp_enum}/{kg_col}{new_amp_eden})'

            # The new step is needed as well
            new_stp = data.loc[cpd, 'Cost calc']
            # Run the cost calculation for the unknown compound
            if excel:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp, 
                                 eamp + new_eamp, excel=excel)
            else:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp)
            # Set the calculated cost
            data.loc[cpd, 'Cost'] = cst
            # And for Excel -- This cost will need to get swapped out later.
            # Also need to check if an OPEX is necessary
            if excel:
                if np.isnan(self.fulldata.loc[(new_stp, cpd), 'OPEX']):
                    data.loc[cpd, 'Cost dyn'] = '=' + ecols['Cost'] +\
                            self.fulldata.loc[(new_stp, cpd), 'rnum']
                else:
                    data.loc[cpd, 'Cost dyn'] = '=' + ecols['Cost'] +\
                            self.fulldata.loc[(new_stp, cpd), 'rnum'] + '+'\
                            + ecols['OPEX'] +\
                            self.fulldata.loc[(new_stp, cpd), 'rnum']

        # Calculate the cost for each material in the reaction
        data['RM cost/kg rxn'] = data['kg/kg rxn']*data['Cost']
        # The product cost will be the sum of all the reactant/solvent costs
        data.loc[prod, 'RM cost/kg rxn'] = data['RM cost/kg rxn'].sum()
        # Set the "Cost" to the calculated value
        data.loc[prod, 'Cost'] = data.loc[prod, 'RM cost/kg rxn']
        # And for Excel
        if excel:
            data['RM cost/kg rxn dyn'] = '=' + ecols['kg/kg rxn'] +\
                    data['rnum'] + '*' + ecols['Cost'] + data['rnum']
            # And for Excel, first we need the rows that are not the product
            mask = data.index != prod
            # Then we need to make a comma-separated list of these rows
            cells = [f'{ecols["RM cost/kg rxn"]}{r}' for r in\
                    data.loc[mask, 'rnum']]
            rs = ','.join(cells)
            # Combine them together into a sum
            data.loc[prod, 'RM cost/kg rxn dyn'] = '=SUM(' + rs + ')'
            # Need to divide by the total number of kgs
            data.loc[prod, 'Cost dyn'] = '=' + ecols['RM cost/kg rxn'] +\
                    data.loc[prod, 'rnum'] + '/' + ecols['kg/kg rxn'] +\
                    data.loc[prod, 'rnum']

        # Calculate % costs for individual rxn
        # = (RM cost/kg rxn)/(RM cost/kg rxn for the rxn product)
        p_rm_cost = data['RM cost/kg rxn']*100/data.loc[prod, 'RM cost/kg rxn']
        data['% RM cost/kg rxn'] = p_rm_cost.values
        # And for Excel
        if excel:
            data['% RM cost/kg rxn dyn'] = '=' + ecols['RM cost/kg rxn'] +\
                    data['rnum'] + '*100/' + ecols['RM cost/kg rxn'] +\
                    data.loc[prod, 'rnum']
        # Remove the % cost for the rxn product
        data.loc[prod, '% RM cost/kg rxn'] = np.nan
        # And for Excel
        if excel:
            data.loc[prod, '% RM cost/kg rxn dyn'] = ''

        # These are the costs for ultimate product
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn 
        # cost. However, for feeder reactions this will get amplified by 
        # each step   
        data['RM cost/kg prod'] = data['RM cost/kg rxn']*amp
        # And for Excel
        if excel:
            data['RM cost/kg prod dyn'] = '=' + ecols['RM cost/kg rxn'] +\
                    data['rnum'] + eamp
        
        # This sets the number of kg of each material per kilogram of product
        # This is done by multiplying the per reaction value by the amplifier
        # This isn't necessary for costing, but we can use it for PMI
        data['kg/kg prod'] = data['kg/kg rxn']*amp
        # And for Excel
        if excel:
            data['kg/kg prod dyn'] = '=' + ecols['kg/kg rxn'] + data['rnum']\
                    + eamp
        
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
    
    def _rxn_data_post(self, excel=False):
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
            # And for Excel
            if excel:
                self.fulldata.loc[(step, prod), 'Cost dyn'] += '+' +\
                        ecols['OPEX'] +\
                        self.fulldata.loc[(step, prod), 'rnum']
                
        # Calculate % overall costs relative to the prod
        self.fulldata['% RM cost/kg prod'] = \
                self.fulldata['RM cost/kg prod']*100/self.cost
        # And for Excel
        if excel:
            self.fulldata['% RM cost/kg prod dyn'] = '=' +\
                    ecols['RM cost/kg prod'] + self.fulldata['rnum'] +\
                    '*100/' + ecols['RM cost/kg rxn'] +\
                    self.fulldata.loc[(step, prod), 'rnum']
        
        # Filter out certain values to simplify full data set
        # Remove the cost and %s for cost-calculated materials
        # This is necessary so that this column adds up to 100% (w/o OPEX)
        mask = ~self.fulldata['Cost calc'].isna()
        self.fulldata.loc[mask, '% RM cost/kg prod'] = np.nan
        if excel:
            self.fulldata.loc[mask, '% RM cost/kg prod dyn'] = ''
        # This filters some of the costs which are simply the sum of raw
        # materials from eariler rxns. The sum of this column will now be
        # equal to the cost of the final product.
        self.fulldata.loc[mask, 'RM cost/kg prod'] = np.nan
        if excel:
            self.fulldata.loc[mask, 'RM cost/kg prod dyn'] = ''
        # This filters out the kg/kg prod values that were calculated, so that
        # the sum of this column is the PMI
        self.fulldata.loc[mask, 'kg/kg prod'] = np.nan
        if excel:
            self.fulldata.loc[mask, 'kg/kg prod dyn'] = ''
        # But we are making 1 kg of final product so that needs to be reset
        self.fulldata.loc[(step, prod), 'kg/kg prod'] = 1.
        if excel:
            self.fulldata.loc[(step, prod), 'kg/kg prod dyn'] = '=1.'

        # PMI Calculations
        # Adding a prefix for display purposes
        # There will be a funny column in this DF with these values...
        self._pre = '*'
        
        # First of all, calculate the PMI for each reaction individually
        gb = self.fulldata.groupby('Step')
        if excel:
            rxn_pmi = gb.agg({'kg/kg rxn':'sum', 'rnum':self._excel_pmi})\
                            .rename({'rnum': 'kg/kg rxn dyn'}, axis=1)\
                            .reset_index()
        else:
            rxn_pmi = gb.agg({'kg/kg rxn':'sum'})\
                            .reset_index()
        rxn_pmi['Compound'] = self._pre + 'Step ' + rxn_pmi['Step'] + ' PMI'
        
        # The full route PMI is not the sum of the above, but is the sum of
        # the 'kg/kg prod' column. We need to make this into a DataFrame to
        # merge with the per reaction values above.
        df_vals = {'kg/kg prod': [self.fulldata['kg/kg prod'].sum()], 
                   'Step': [self._fp_idx],
                   'Compound': [self._pre*2 + 'Full Route PMI'],
                   }
        if excel:
            mask = ~self.fulldata['kg/kg prod'].isna()
            all_cells = [f'{ecols["kg/kg prod"]}{i}' for i in\
                         self.fulldata.loc[mask, 'rnum']]
            cell_sum = '=SUM(' + ','.join(all_cells) + ')'
            df_vals['kg/kg prod dyn'] = [cell_sum]
        full_pmi = pd.DataFrame(df_vals)

        # Merge the per-reaction and full PMI
        self.pmi = pd.concat([rxn_pmi, full_pmi], 
                             sort=False).set_index('Step')

    def _excel_pmi(self, col):
        # A function for creating the dynamic Excel cells for PMI
        cells = [f'{ecols["kg/kg rxn"]}{i}' for i in col]
        return "=SUM(" + ','.join(cells) + ")"

    def results(self, style='compact'):
        '''Preps an output DataFrame for the results method.
        '''
        # For compact display, these are the most important columns
        comp_col = ['Cost', 'Equiv',] 
        if self.fulldata.Volumes.any():
            comp_col.extend(['Volumes', 'Sol Recyc',])
        if self.fulldata.OPEX.any():
            comp_col.append('OPEX') 
        comp_col.extend(['kg/kg rxn', 'RM cost/kg rxn', '% RM cost/kg rxn',
                    'kg/kg prod', 'RM cost/kg prod', '% RM cost/kg prod'])

        # For full display, don't show the Excel columns
        ecol_mask = ~self.fulldata.columns.str.contains('dyn|rnum')
        no_ecol = self.fulldata.columns[ecol_mask]
        
        # Combine the fulldata and pmi DataFrames
        fd = self._df_combine()

        # Display the DataFrames for different permutations of kwargs
        # The fillna removes NaN from the displayed tables in Notebook, but
        # may goof up printing
        if style == 'full':
            fd = fd[no_ecol]
        elif style == 'compact':
            fd = fd[comp_col]

        return fd

    def _df_combine(self, ):
        '''Combine the fulldata and pmi DataFrames for saving/exporting.
        '''
        # Copy the original DFs, and remove indexing
        fd = self.fulldata.reset_index()
        pmi = self.pmi.reset_index()
        
        # Concats the fulldata and pmi based on step. The groupby function
        # will sort by "Step" as it processes things. This will need to be
        # modified if this behavior is undesired, but as is, it is not trivial
        # to make this fix. The apply function doubles up the "Step" column.
        # So one needs to be removed.
        gb = fd.groupby('Step')
        concated = gb.apply(self.__fd_pmi_concat, pmi)
        concated = concated.drop('Step', axis=1)\
                            .reset_index().drop('level_1', axis=1)
        
        # Reset the index and return the DF. No sorting is necessary here.
        return concated.set_index(['Step', 'Compound'])                

    def __fd_pmi_concat(self, df, pmi):
        '''Concats the fulldata and pmi DataFrames.

        This is only used in a DataFrame.apply call in the _df_combine
        method.'''
        # The step number
        step = df['Step'].iloc[0]
        # Mask out the correct PMI values
        pmi_mask = pmi['Step'] == step
        pmi_small = pmi[pmi_mask]
        # Combine the route step info with the PMI. In this way, the pmi will
        # be at the end.
        return pd.concat([df, pmi_small], ignore_index=True)

    def excel(self, ):
        '''Prepare a DataFrame for Excel export.

        See `excel_save` method for the docstring info.
        '''
        # Run the cost calculation again, but using the excel keyword
        self.calc_cost(excel=True)

        # Combine the fulldata and pmi DataFrames
        fd = self._df_combine()
        # Move the Notes columns to the end of the combined DataFrame
        mat_note = fd.pop('Material Notes')
        fd.insert(fd.shape[1], 'Material Notes', mat_note)
        rxn_note = fd.pop('Reaction Notes')
        fd.insert(fd.shape[1], 'Reaction Notes', rxn_note)

        # Convert the unique row ID with an Excel sheet row number
        nrows = fd.shape[0]
        for r, n in zip(fd['rnum'], np.arange(2, nrows+2)):
            if isinstance(r, float):
                continue
            for col in fd:
                if 'dyn' in col:
                    fd[col] = fd[col].str.replace(r, str(n))
                    
        # Set dynamic cost for the "Cost" of calculated products
        mask = fd['Cost dyn'] != ''
        fd.loc[mask, 'Cost'] = fd.loc[mask, 'Cost dyn']
        
        # Drop the non-dynamic columns
        d_cols = ['Cost dyn', 'kg/kg rxn', 'RM cost/kg rxn', '% RM cost/kg rxn', 
                  'kg/kg prod', 'RM cost/kg prod', '% RM cost/kg prod',
                  'rnum']
        fd = fd.drop(d_cols, axis=1)
        
        # Rename columns to remove "dyn" suffix
        col_str = fd.columns.str.replace(' dyn', '')
        new_cols = {c:cn for (c, cn) in zip(fd.columns, col_str)}
        fd = fd.rename(new_cols, axis=1)
        
        # Rerun the cost calculation without the excel stuff to get rid of all
        # the other columns
        self.calc_cost(excel=False)

        return fd


