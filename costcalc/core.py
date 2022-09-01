'''
Chemical Reaction Cost Calculation and Excel Export Algorithms.
Adapted from the Excel spreadsheets prepared by Saeed Ahmad, PhD.
(C) Ryan Nelson
'''
import time

import numpy as np
import pandas as pd

# Setting the display to print function
# This gets changed for Jupyter Notebook/IPython sessions, so that DataFrames
# are displayed in a fancier format
try:
    from IPython.display import display as disp
    from IPython.display import Javascript
except:
    disp = print

# Column name mapper: This is setting variable names for many of the important
# DataFrame columns. 
# Reactions Table Columns
rxn_stp = 'Step' # The step labels
rxn_cpd = 'Compound' # The compound names column, same for materials table
rxn_eq = 'Equiv' # Molar equivalents of reagent
rxn_ms = 'Mass' # Mass of material used, typically kg, just be consistent
rxn_vol = 'Volumes' # Volumes of solvent L/kg of limiting reagent (usually)
rxn_rel = 'Relative' # Limiting reagent name for volumes->mass calculation
rxn_rcy = 'Sol Recyc' # Fractional amount of solvent that is recycled
rxn_cst = 'Cost step' # Step where cost is calculated
rxn_opx = 'OPEX' # Operational expenditures in $/kg, only for rxn product
notes = 'Notes' # The notes column name, same for materials table
# Materials Table Columns
mat_mw = 'MW' # Molecular weight
mat_den = 'Density' # Density, only used for solvents
mat_cst = 'Cost' # The estimated material price/cost ($/kg)
# Calculated columns
rxn_kg = 'kg/kg rxn' # kg of materials used per kg reaction/step product
rxn_rmc = 'RM cost/kg rxn' # Cost of a material to make 1 kg rxn product
rxn_rmp = '% RM cost/kg rxn' # % material cost per 1 kg rxn product
prd_kg = 'kg/kg prod' # kg of material used per kg of route product
prd_rmc = 'RM cost/kg prod' # Cost of a material to make 1 kg route prod
prd_rmp = '% RM cost/kg prod' # % material cost per 1 kg route product
# Excel temporary column names, these are renamed on Excel export
# These are the same as the calculated columns except for the suffix
suffix = ' dyn'
rnum = 'rnum' # Row number column for indexing Excel cells
dyn_cst = mat_cst + suffix # Dynamic product cost
dyn_rkg = rxn_kg + suffix # dynamic kg/kg rxn
dyn_rrmc = rxn_rmc + suffix # dynamic RM cost/kg rxn
dyn_rrmp = rxn_rmp + suffix # dynamic %RM cost/kg rxn
dyn_pkg = prd_kg + suffix # dynamic kg/kg prod
dyn_prmc = prd_rmc + suffix # dynamic RM cost/kg prod
dyn_prmp = prd_rmp + suffix # dynamic %RM cost/kg prod


# Excel Column dictionary
# This dictionary will be used for creating dynamic excel sheets, for which
# I'll need to know the DataFrame<->Excel column name mapping
ecols = {rxn_stp:'A', rxn_cpd:'B', mat_mw:'C', mat_den:'D', mat_cst:'E',
        rxn_eq:'F', rxn_vol:'G', rxn_rel:'H', rxn_rcy:'I', rxn_cst:'J',
        rxn_opx:'K', rxn_kg:'L', rxn_rmc:'M', rxn_rmp:'N', prd_kg:'O',
        prd_rmc:'P', prd_rmp:'Q', }


class CoreCost(object):
    def __init__(self, materials, rxns, final_prod, disp_err_df=False):
        # We need to store the input values/DataFrames
        self.rxns = rxns
        self.materials = materials
        self.final_prod = final_prod        
        self._disp_err_df = disp_err_df

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
        mat_keeps = [rxn_cpd, mat_mw, mat_den, mat_cst, notes]

        rxn_keeps = [rxn_stp, rxn_cpd, rxn_eq, rxn_vol, rxn_rel,
                    rxn_rcy, rxn_cst, rxn_opx, notes]
        # Check if an rxn_ms column is present. This is used to calculate
        # equivalents. If this is present, add this column to our "keep" list
        amt = False
        if rxn_ms in self.rxns.columns:
            amt = True 
            rxn_keeps = rxn_keeps[:2] + [rxn_ms] + rxn_keeps[2:]

        fulldata = pd.merge(self.materials[mat_keeps], self.rxns[rxn_keeps],
                            on=rxn_cpd, how='right')\
                            .rename({'Notes_x': 'Material Notes',
                                'Notes_y': 'Reaction Notes'}, axis=1)

        # Find the step number for the final product. 
        fp_mask = fulldata.Compound == self.final_prod
        self._fp_idx = fulldata.loc[fp_mask, rxn_cst].iloc[0]

        # If the rxn_ms column is present, then you will need calculate the
        # "Equiv" based on the Amount given
        if amt:
            # Group everything by "Step" column
            grp = fulldata.groupby(rxn_stp)
            for idx, sb_grp in grp:
                # If there are no amounts for a given Step, then skip
                if ~sb_grp[rxn_ms].any():
                    continue
                # Mask out only the values that have amounts. This is
                # important in case a mixture of mass and equiv is used.
                # Also check to make sure this isn't a solvent
                mask = ~sb_grp[rxn_ms].isna() & sb_grp[rxn_vol].isna()
                # Calculate the Equivalents
                sb_grp.loc[mask, rxn_eq] = sb_grp.loc[mask, mat_mw]*\
                                        sb_grp.loc[mask, rxn_ms]
                # Assume the first entry is the limiting reagent. Normalize
                # the rest of the equivalents to that value
                sb_grp.loc[mask, rxn_eq] = sb_grp.loc[mask, rxn_eq]/\
                                       (sb_grp.iloc[0][rxn_eq])
                # Set the values in the full DataFrame
                mask = fulldata[rxn_stp] == idx
                fulldata.loc[mask, rxn_eq] = sb_grp[rxn_eq]
            # Remove the Amount column
            fulldata.drop(rxn_ms, axis=1, inplace=True)

        # Add an indexing column, which will be used for sorting purposes
        fulldata['idx_col'] = np.arange(fulldata.shape[0])
        # Set MultiIndex
        fulldata.set_index(['idx_col', rxn_stp, rxn_cpd], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not
        # copies
        fulldata = fulldata.sort_index()
        # Remove sorting index, drop column
        fulldata = fulldata.reset_index('idx_col').drop('idx_col', axis=1)
        
        # Save the full data set
        self.fulldata = fulldata

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
        cost_recalc_mask = ~self.fulldata[rxn_cst].isna()
        self.fulldata.loc[cost_recalc_mask, mat_cst] = np.nan

        # Create or clear a bunch of columns that will be populated during 
        # cost calculation. 
        empty_cols = [rxn_kg, rxn_rmc, rxn_rmp, prd_kg, prd_rmc, prd_rmp,]
        
        for col in empty_cols:
            self.fulldata[col] = np.nan
        # For dynamic Excel
        if excel:
            excel_cols = [dyn_cst, dyn_rkg, dyn_rrmc, dyn_rrmp, dyn_pkg,
                    dyn_prmc, dyn_prmp,]

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
        mw_mask = self.fulldata[mat_mw].isna()
        if mw_mask.any():
            err_line = 'Missing MW in input file! \n'\
                'You are missing one or more MW values.\n'\
                'Make sure all reaction compounds are defined in the'\
                ' materials table.'
            print(err_line)
            if self._disp_err_df:
                disp(self.fulldata[mw_mask])
            raise ValueError(err_line)

        # Check for duplicated materials. This will probably be a big issue
        # with two materials sheets.
        dup_cpds = self.materials[rxn_cpd].duplicated()
        if dup_cpds.any():
            err_line = 'You have a duplicate material!!!\n'\
                    'Check that you do not have two or more of the same\n'\
                    'compound name in your material sheet(s).'
            print(err_line)
            if self._disp_err_df:
                print("These compounds are the duplicated compounds.")
                disp(self.materials.loc[dup_cpds, rxn_cpd])
            raise ValueError(err_line)

        # Check for duplicated materials in a single reaction.
        # When you select a single value from a reaction, you'll get a series
        # and not a float, e.g.
        dup_rxn = self.fulldata.loc[(self._fp_idx, self.final_prod), mat_mw]
        if isinstance(dup_rxn, pd.Series):
            err_line = 'You have a duplicated material in a single '\
                    'reaction step.\nPlease combine multiple uses of'\
                    ' a material into one entry.'
            print(err_line)
            if self._disp_err_df:
                print('Check these lines:')
                gb = self.fulldata.groupby([rxn_stp, rxn_cpd])
                for prod, group in gb:
                    if group.shape[0] > 1:
                        print('Step: ' + prod[0] + '; Cpd: ' + prod[1])
            raise ValueError(err_line)
            
        # Check for a missing cost, which is not being calculated
        # This is a tricky error because the costing will run just fine with 
        # NaNs. Check for this by looking for NaN in both Cost and Calc columns
        cost_mask = (self.fulldata[rxn_cst].isna() & \
                    self.fulldata[mat_cst].isna())
        if cost_mask.any():
            err_line = 'You are missing a necessary material cost!!\n'\
                   f'You may need to indicate a Step in the "{rxn_cst}"'\
                   ' column.' 
            print(err_line)
            if self._disp_err_df:
                print('Check these columns.')
                disp(self.fulldata.loc[cost_mask, [mat_cst, rxn_cst]])
            raise ValueError(err_line)
        
        # Check that all the solvent information is given
        sol_mask = ~self.fulldata[rxn_vol].isna() 
        sol_chk = self.fulldata.loc[sol_mask, rxn_rel].isna() | \
                self.fulldata.loc[sol_mask, mat_den].isna() | \
                self.fulldata.loc[sol_mask, rxn_rcy].isna()
        # If anything is missing, print a note
        if sol_chk.any():
            err_line = 'You are missing some solvent information.'
            print(err_line)
            if self._disp_err_df:
                print('Check the following lines.')
                sol_cols = [mat_den, rxn_vol, rxn_rel, rxn_rcy]
                sol_mask2 = sol_mask & sol_chk
                disp(self.fulldata.loc[sol_mask2, sol_cols])
            raise ValueError(err_line)
        
        # Check to make sure that the "Relative" compound for a solvent
        # is acutally contained in the step
        msg = False
        sol_mask = ~self.fulldata[rxn_rel].isna()
        bad = []
        for cpd in self.fulldata[sol_mask].itertuples():
            new_idx = (cpd.Index[0], getattr(cpd, rxn_rel) )
            if new_idx not in self.fulldata.index:
                msg = True
                bad.append( cpd.Index )
        if msg:
            err_line = 'You have a bad solvent recycle entry.\n'\
                    f'One of your "{rxn_rel}" compounds is not correct.'
            print(err_line)
            if self._disp_err_df:
                print(f'Check the following Step/Solvents:')
                [print(f'Step {s}: Solvent {sol}') for s, sol in bad]
            raise ValueError(err_line)
        
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
            self.fulldata[rnum] = [f'r{i:0{nrowdig:d}d}' for i in range(nrows)]
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
        data[rxn_kg] = data[rxn_eq]*data[mat_mw]
        # And for Excel
        if excel:
            data[dyn_rkg] = '=' + ecols[rxn_eq] + data[rnum] + '*'\
                    + ecols[mat_mw] + data[rnum]
        # For Excel, normalize the data here. This will get overwritten for
        # solvents. This is a little weird for the product because it will be
        # the same in numerator and denominator, but that will keep things
        # fully interactive.
        if excel:
            data[dyn_rkg] += '/(' + ecols[rxn_eq] +\
                    data.loc[prod,rnum] + '*' + ecols[mat_mw] +\
                    data.loc[prod, rnum] + ')'

        # Amount of solvent
        # First figure out which materials are solvents 
        mask = ~data[rxn_vol].isna()
        sols = data[mask].copy()
        if sols.size != 0:
            # Find the material that the volumes are relative to (taking only
            # the first one... This might not be great.
            cpd_rel = sols[rxn_rel][0]
            # What is the kg of relative cpd?
            amt_rel = data.loc[cpd_rel, rxn_kg]
            if excel:
                amt_rel_e = data.loc[cpd_rel, rnum]
            # Calculate the mass of solvent. Take into account the solvent
            # recycyling 
            # kg sol = Volume*Density*(1-Recycle)*(kg SM)
            sols[rxn_kg] = sols[rxn_vol]*sols[mat_den]*\
                    (1 - sols[rxn_rcy])*amt_rel
            # And for Excel
            if excel:
                sols[dyn_rkg] = '=' + ecols[rxn_vol] +\
                        sols[rnum] + '*' + ecols[mat_den] +\
                        sols[rnum] + '*' + '(1 - ' + ecols[rxn_rcy]\
                        + sols[rnum] + ')*' + ecols[rxn_kg] +\
                        amt_rel_e
            data[mask] = sols

        # Normalize the kg of reaction
        data[rxn_kg] /= data.loc[prod, rxn_kg]

        # Calculate unknown costs. Looks for any empty values in the "Cost" 
        # column. Don't use the rxn_cst column directly, because some of
        # the costs may have been manually set using the `value_mod` method
        # unknown_cost = ~data[rxn_cst].isna()
        unknown_cost = data[mat_cst].isna()
        # This is recursive. The final cost per kg of product will be
        # amplified by each subsequent step, which is where the new_amp
        # calculation comes into play
        for cpd, row in data.loc[unknown_cost].iterrows():
            # Don't do this for the product of the current reaction
            if cpd == prod: 
                continue
            # The amounts needed will be amplified by the appropriate kg ratio.
            # Set that ratio
            new_amp = data.loc[cpd, rxn_kg]
            # And for dynamic excel output
            if excel:
                new_amp_enum = data.loc[cpd, rnum]
                new_amp_eden = data.loc[prod, rnum]
                kg_col = ecols[rxn_kg]
            # This will be '*(C#/D#)'. This is different for Excel because we
            # can't normalize the kg/kg rxn column
                new_eamp = f'*({kg_col}{new_amp_enum}/{kg_col}{new_amp_eden})'

            # The new step is needed as well
            new_stp = data.loc[cpd, rxn_cst]
            # Run the cost calculation for the unknown compound
            if excel:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp, 
                                 eamp + new_eamp, excel=excel)
            else:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp)
            # Set the calculated cost
            data.loc[cpd, mat_cst] = cst
            # And for Excel -- This cost will need to get swapped out later.
            # Also need to check if an OPEX is necessary
            if excel:
                if np.isnan(self.fulldata.loc[(new_stp, cpd), rxn_opx]):
                    data.loc[cpd, dyn_cst] = '=' + ecols[mat_cst] +\
                            self.fulldata.loc[(new_stp, cpd), rnum]
                else:
                    data.loc[cpd, dyn_cst] = '=' + ecols[mat_cst] +\
                            self.fulldata.loc[(new_stp, cpd), rnum] + '+'\
                            + ecols[rxn_opx] +\
                            self.fulldata.loc[(new_stp, cpd), rnum]

        # Calculate the cost for each material in the reaction
        data[rxn_rmc] = data[rxn_kg]*data[mat_cst]
        # The product cost will be the sum of all the reactant/solvent costs
        data.loc[prod, rxn_rmc] = data[rxn_rmc].sum()
        # Set the "Cost" to the calculated value
        data.loc[prod, mat_cst] = data.loc[prod, rxn_rmc]
        # And for Excel
        if excel:
            data[dyn_rrmc] = '=' + ecols[rxn_kg] + data[rnum] + '*' + \
                            ecols[mat_cst] + data[rnum]
            # And for Excel, first we need the rows that are not the product
            mask = data.index != prod
            # Then we need to make a comma-separated list of these rows
            cells = [f'{ecols[rxn_rmc]}{r}' for r in data.loc[mask, rnum]]
            rs = ','.join(cells)
            # Combine them together into a sum
            data.loc[prod, dyn_rrmc] = '=SUM(' + rs + ')'
            # Need to divide by the total number of kgs
            data.loc[prod, dyn_cst] = '=' + ecols[rxn_rmc] +\
                    data.loc[prod, rnum] + '/' + ecols[rxn_kg] +\
                    data.loc[prod, rnum]

        # Calculate % costs for individual rxn
        # = (RM cost/kg rxn)/(RM cost/kg rxn for the rxn product)
        p_rm_cost = data[rxn_rmc]*100/data.loc[prod, rxn_rmc]
        data[rxn_rmp] = p_rm_cost.values
        # And for Excel
        if excel:
            data[dyn_rrmp] = '=' + ecols[rxn_rmc] +\
                    data[rnum] + '*100/' + ecols[rxn_rmc] +\
                    data.loc[prod, rnum]
        # Remove the % cost for the rxn product
        data.loc[prod, rxn_rmp] = np.nan
        # And for Excel
        if excel:
            data.loc[prod, dyn_rrmp] = ''

        # These are the costs for ultimate product
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn 
        # cost. However, for feeder reactions this will get amplified by 
        # each step   
        data[prd_rmc] = data[rxn_rmc]*amp
        # And for Excel
        if excel:
            data[dyn_prmc] = '=' + ecols[rxn_rmc] + data[rnum] + eamp
        
        # This sets the number of kg of each material per kilogram of product
        # This is done by multiplying the per reaction value by the amplifier
        # This isn't necessary for costing, but we can use it for PMI
        data[prd_kg] = data[rxn_kg]*amp
        # And for Excel
        if excel:
            data[dyn_pkg] = '=' + ecols[rxn_kg] + data[rnum] + eamp
        
        # Set the values in the big DataFrame with this slice. This goofy call
        # is necessary to make sure the data are set correctly. 
        self.fulldata.loc[(step, data.index), data.columns] = data.values
        
        # Return the calculated product cost, which is required for the 
        # recursive nature of the algorithm. In addition, an optional OPEX
        # may be added to take into account production costs of the cpd
        if np.isnan(data.loc[prod, rxn_opx]):
            return data.loc[prod, rxn_rmc]
        else:
            return data.loc[prod, rxn_rmc] + data.loc[prod, rxn_opx]
    
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
        opex = self.fulldata.loc[(step, prod), rxn_opx]
        if not np.isnan(opex):
            self.fulldata.loc[(step, prod), mat_cst] = self.cost
            # And for Excel
            if excel:
                self.fulldata.loc[(step, prod), dyn_cst] += '+' + \
                        ecols[rxn_opx] + \
                        self.fulldata.loc[(step, prod), rnum]
                
        # Calculate % overall costs relative to the prod
        self.fulldata[prd_rmp] = self.fulldata[prd_rmc]*100/self.cost

        # And for Excel
        if excel:
            self.fulldata[dyn_prmp] = '=' +\
                    ecols[prd_rmc] + self.fulldata[rnum] +\
                    '*100/' + ecols[rxn_rmc] +\
                    self.fulldata.loc[(step, prod), rnum]
        
        # Filter out certain values to simplify full data set
        # Remove the cost and %s for cost-calculated materials
        # This is necessary so that this column adds up to 100% (w/o OPEX)
        mask = ~self.fulldata[rxn_cst].isna()
        self.fulldata.loc[mask, prd_rmp] = np.nan
        if excel:
            self.fulldata.loc[mask, dyn_prmp] = ''

        # This filters some of the costs which are simply the sum of raw
        # materials from earlier rxns. The sum of this column will now be
        # equal to the cost of the final product.
        self.fulldata.loc[mask, prd_rmc] = np.nan
        if excel:
            self.fulldata.loc[mask, dyn_prmc] = ''

        # This filters out the kg/kg prod values that were calculated, so that
        # the sum of this column is the PMI
        self.fulldata.loc[mask, prd_kg] = np.nan
        if excel:
            self.fulldata.loc[mask, dyn_pkg] = ''
        # But we are making 1 kg of final product so that needs to be reset
        self.fulldata.loc[(step, prod), prd_kg] = 1.
        if excel:
            self.fulldata.loc[(step, prod), dyn_pkg] = '=1.'

        # PMI Calculations
        # Adding a prefix for display purposes
        # There will be a funny column in this DF with these values...
        self._pre = '*'
        
        # First of all, calculate the PMI for each reaction individually
        gb = self.fulldata.groupby(rxn_stp)
        if excel:
            rxn_pmi = gb.agg({rxn_kg:'sum', rnum:self._excel_pmi})\
                            .rename({rnum: dyn_rkg}, axis=1)\
                            .reset_index()
        else:
            rxn_pmi = gb.agg({rxn_kg:'sum'}).reset_index()

        rxn_pmi[rxn_cpd] = self._pre + 'Step ' + rxn_pmi[rxn_stp] + ' PMI'
        
        # The full route PMI is not the sum of the above, but is the sum of
        # the prd_kg column. We need to make this into a DataFrame to
        # merge with the per reaction values above.
        df_vals = {prd_kg: [self.fulldata[prd_kg].sum()], 
                   rxn_stp: [self._fp_idx],
                   rxn_cpd: [self._pre*2 + 'Full Route PMI'],
                   }
        if excel:
            mask = ~self.fulldata[prd_kg].isna()
            all_cells = [f'{ecols[prd_kg]}{i}' for i in\
                         self.fulldata.loc[mask, rnum]]
            cell_sum = '=SUM(' + ','.join(all_cells) + ')'
            df_vals[dyn_pkg] = [cell_sum]
        full_pmi = pd.DataFrame(df_vals)

        # Merge the per-reaction and full PMI
        self.pmi = pd.concat([rxn_pmi, full_pmi], sort=False)\
                .set_index(rxn_stp)

    def _excel_pmi(self, col):
        # A function for creating the dynamic Excel cells for PMI
        cells = [f'{ecols[rxn_kg]}{i}' for i in col]
        return "=SUM(" + ','.join(cells) + ")"

    def results(self, style='compact'):
        '''Preps an output DataFrame for the results method.
        '''
        # For compact display, these are the most important columns
        comp_col = [mat_cst, rxn_eq,] 
        if self.fulldata.Volumes.any():
            comp_col.extend([rxn_vol, rxn_rcy,])
        if self.fulldata.OPEX.any():
            comp_col.append(rxn_opx) 
        comp_col.extend([rxn_kg, rxn_rmc, rxn_rmp, prd_kg, prd_rmc, prd_rmp])

        # For full display, don't show the Excel columns
        ecol_mask = ~self.fulldata.columns.str.contains(suffix+'|'+rnum)
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
        gb = fd.groupby(rxn_stp)
        concated = gb.apply(self.__fd_pmi_concat, pmi)
        concated = concated.drop(rxn_stp, axis=1)\
                            .reset_index()\
                            .drop('level_1', axis=1)
        
        # Reset the index and return the DF. No sorting is necessary here.
        return concated.set_index([rxn_stp, rxn_cpd])                

    def __fd_pmi_concat(self, df, pmi):
        '''Concats the fulldata and pmi DataFrames.

        This is only used in a DataFrame.apply call in the _df_combine
        method.'''
        # The step number
        step = df[rxn_stp].iloc[0]
        # Mask out the correct PMI values
        pmi_mask = pmi[rxn_stp] == step
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
        for r, n in zip(fd[rnum], np.arange(2, nrows+2)):
            if isinstance(r, float):
                continue
            for col in fd:
                if suffix in col:
                    fd[col] = fd[col].str.replace(r, str(n))
                    
        # Set dynamic cost for the "Cost" of calculated products
        mask = fd[dyn_cst] != ''
        fd.loc[mask, mat_cst] = fd.loc[mask, dyn_cst]
        
        # Drop the non-dynamic columns
        d_cols = [dyn_cst, rxn_kg, rxn_rmc, rxn_rmp, prd_kg, prd_rmc, prd_rmp,
                rnum]
        fd = fd.drop(d_cols, axis=1)
        
        # Rename columns to remove "dyn" suffix
        col_str = fd.columns.str.replace(suffix, '')
        new_cols = {c:cn for (c, cn) in zip(fd.columns, col_str)}
        fd = fd.rename(new_cols, axis=1)
        
        # Rerun the cost calculation without the excel stuff to get rid of all
        # the other columns
        self.calc_cost(excel=False)

        return fd


