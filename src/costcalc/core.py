import numpy as np
import pandas as pd

from .constants import *
from .exceptions import *

class CoreCost(object):
    def __init__(self, materials, rxns, final_prod, disp_err_df=False):
        '''Core class for calculating route costs and PMI.

        Parameters
        ----------
        materials : DataFrame
            A DataFrame containing  materials properties. At a minimum, this
            table must include the materials used in the route; however, it
            can contain an arbitrary number of additional materials as well,
            as long as the names are unique. This DataFrame must contain the
            following columns:

                * Compound : str, the compound names.
                * MW : float, the compound molecular weights.
                * Density : float, the density (g/mL) of the compound.
                Although this column must be present, density values are only
                necessary for reaction solvents.
                * Cost : float, the cost/price of the compound in $/kg. Can be
                left blank for unknown costs, like intermediates and final
                products
                * Notes : str, optional notes about the compound.

            As noted above, the `Cost` column can be left blank for materials
            that will have their costs calculated. However, placeholder prices
            can be put here, if desired. These prices will be ignored if the
            compounds RM cost is being calculated in the route.

        rxns : DataFrame
            A DataFrame containing reaction information. This DataFrame must
            contain the following columns:

                * Step : str, a unique identifier for each step such as "1" or
                "1a"
                * Compound: str, the compound names
                * Equiv : float, the equivalents used for each compound, for
                reaction products, this value should be the fractional yield
                (e.g. 75% yield is 0.75 equiv)
                * Mass : float, (Optional), the mass of material used. This is
                used to calculate equivalents. If used, the first compound
                for each reaction must be the limiting reagent.
                * Volumes : float, volume equivalents of solvent (L/kg)
                * Relative : str, the compound name to use as the reference
                for converting solvent volumes to kg. This must correspond to
                one of the Compound strs for the step
                * Sol Recyc : float, fractional percentage of solvent that can
                be recycled. E.g. 75% of solvent can be recycled = 0.75
                * Cost step : str, the "Step" where compound costs are
                calculated. For reaction products, this will be the current
                "Step". This is used to trace the reaction network.
                * OPEX : float, a $/kg-step charge that is (optionally) added
                to the final raw-material cost of a reaction product. Is only
                valid for reaction products.

        final_prod : str
            Name of the final product in the route.

        disp_err_df : Boolean (False)
            If a `CostError` exception is raised, this controls whether the
            offending section of the DataFrame will be printed along with the
            error. This is useful for debugging, but is typically set to
            False, so that potentially sensitive information is not printed to
            the error logs.

        Attributes
        ----------
        final_prod : str
            The name of the overall final product of this route.
            
        rxns : DataFrame
            A DataFrame describing the reactions defined for the route. 
            
        materials : DataFrame
            A DataFrame describing all of the known materials. This will contain
            all of the materials from both the `materials_file` and the
            `alt_mat_file`.
            
        cost : float
            The final cost of the described route. This value is only present
            after the `calc_cost` method is called. 
            
        fulldata : DataFrame
            A combined DataFrame containing all of the reaction, material, and
            costing related values for the given route. 

        pmi : DataFrame
            A DataGrame containing the PMI for each reaction and the overall
            route. There will be an extra column with names used for sorting
            purposes; this column can be ignored.

        Notes
        -----
        Some automated value correction is done on the input DataFrames, such
        as removing blank lines and extra whitespace in compound names.
        However, other common errors will raise `CostError` exceptions, and
        these problems will require manual error correction from the user.
        '''
        # We need to store the input values/DataFrames
        self.rxns = rxns
        self.materials = materials
        self.final_prod = final_prod        
        self._disp_err_df = disp_err_df

        # Correct typical input problems. 
        self._problem_correct()

        # Combine the reaction/materials sheets, add new columns
        self.rxn_data_setup()

        # Look for other common errors in the input. Throw an error if found
        # These errors must be manually fixed by the user, which is why they
        # are not covered in problem correction method
        self._sanity_check()

    def _problem_correct(self, ):
        '''Correct typical input problems in the input DataFrames.

        These are problems that tend to cause problems before merging using
        `rxn_data_setup` method, so they are automatically modified here. 
        '''
        rxn = self.rxns.copy()
        mat = self.materials.copy()
        # Drop lines that are still empty, these cause all sorts of
        # problems. Assume rxn/materials tables should have compound names
        # for valid entries
        rxn = rxn[ ~rxn[RXN_CPD].isna() ]
        mat = mat[ ~mat[RXN_CPD].isna() ]

        # Remove white space from before/after these columns
        # This is another tricky problem because trailing white space for
        # example is very hard to notice
        rxn_col = [RXN_STP, RXN_CPD, RXN_REL, RXN_CST]
        rxn[rxn_col] = rxn[rxn_col].apply(lambda x: x.str.strip())
        mat_col = [RXN_CPD, ]
        mat[mat_col] = mat[mat_col].apply(lambda x: x.str.strip())

        self.rxns = rxn.copy()
        self.materials = mat.copy()

    def rxn_data_setup(self, ):
        '''Setup the full data set for the upcoming cost calculations. 
        
        This method merges the materials and reaction sheets into a combined
        DataFrame called `fulldata`. It also serves as a "reset" of sorts
        after a costing calculation, because it uses the input materials and
        reactions DataFrames to rebuild the `fulldata` table. This is
        important if you are modifiying values, for example, and you want to
        put the system back into its starting point. 
        '''
        # Merge the materials and reaction DataFrames. A few columns are
        # dropped, which are not necessary for calculations. The merge happens
        # on the rxns DataFrame ('right'), which means that missing materials
        # will be fairly obvious (no MW, e.g.).
        mat_keeps = [RXN_CPD, MAT_MW, MAT_DEN, MAT_CST, NOTES]

        rxn_keeps = [RXN_STP, RXN_CPD, RXN_EQ, RXN_VOL, RXN_REL,
                    RXN_RCY, RXN_CST, RXN_OPX, NOTES]
        # Check if an RXN_MS column is present. This is used to calculate
        # equivalents. If this is present, add this column to our "keep" list
        amt = False
        if RXN_MS in self.rxns.columns:
            amt = True 
            rxn_keeps = rxn_keeps[:2] + [RXN_MS] + rxn_keeps[2:]

        fulldata = pd.merge(self.materials[mat_keeps], self.rxns[rxn_keeps],
                            on=RXN_CPD, how='right')\
                            .rename({'Notes_x': 'Material Notes',
                                'Notes_y': 'Reaction Notes'}, axis=1)

        # Find the step number for the final product. 
        fp_mask = fulldata.Compound == self.final_prod
        self._fp_idx = fulldata.loc[fp_mask, RXN_CST].iloc[0]

        # If the RXN_MS column is present, then you will need calculate the
        # "Equiv" based on the Amount given
        if amt:
            # Group everything by "Step" column
            grp = fulldata.groupby(RXN_STP)
            for idx, sb_grp in grp:
                # If there are no amounts for a given Step, then skip
                if ~sb_grp[RXN_MS].any():
                    continue
                # Mask out only the values that have amounts. This is
                # important in case a mixture of mass and equiv is used.
                # Also check to make sure this isn't a solvent
                mask = ~sb_grp[RXN_MS].isna() & sb_grp[RXN_VOL].isna()
                # Calculate the Equivalents
                sb_grp.loc[mask, RXN_EQ] = sb_grp.loc[mask, MAT_MW]*\
                                        sb_grp.loc[mask, RXN_MS]
                # Assume the first entry is the limiting reagent. Normalize
                # the rest of the equivalents to that value
                sb_grp.loc[mask, RXN_EQ] = sb_grp.loc[mask, RXN_EQ]/\
                                       (sb_grp.iloc[0][RXN_EQ])
                # Set the values in the full DataFrame
                mask = fulldata[RXN_STP] == idx
                fulldata.loc[mask, RXN_EQ] = sb_grp[RXN_EQ]
            # Remove the Amount column
            fulldata.drop(RXN_MS, axis=1, inplace=True)

        # Add an indexing column, which will be used for sorting purposes
        fulldata['idx_col'] = np.arange(fulldata.shape[0])
        # Set MultiIndex
        fulldata.set_index(['idx_col', RXN_STP, RXN_CPD], inplace=True)
        # This is necessary so that slices of the DataFrame are views and not
        # copies
        fulldata = fulldata.sort_index()
        # Remove sorting index, drop column
        fulldata = fulldata.reset_index('idx_col').drop('idx_col', axis=1)
        
        # Save the full data set
        self.fulldata = fulldata

        # Add a modified variable placeholder. This will store modified
        # values for the helper functions
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
        cost_recalc_mask = ~self.fulldata[RXN_CST].isna()
        self.fulldata.loc[cost_recalc_mask, MAT_CST] = np.nan

        # Create or clear a bunch of columns that will be populated during 
        # cost calculation. 
        empty_cols = [RXN_KG, RXN_RMC, RXN_RMP, PRD_KG, PRD_RMC, PRD_RMP,]
        self.fulldata[ empty_cols ] = np.nan

        # For dynamic Excel
        if excel or self.fulldata.columns.str.contains(SUFFIX).any():
            excel_cols = [DYN_CST, DYN_RKG, DYN_RRMC, DYN_RRMP, DYN_PKG,
                    DYN_PRMC, DYN_PRMP,]
            self.fulldata[ excel_cols ] = ''

        # Change modified values that were set by helper functions
        for mod in self._mod_vals:
            self._set_val(*mod)

    def _sanity_check(self, ):
        '''Run some sanity checks on the DataFrames to catch common errors.

        Any problems found here will raise custom `CostError` exceptions. (See
        `costcalc.expections` module for the `CostError` expection class and
        error text.) Some of the errors are tricky, and the typical Python
        error output is unhelpful. These errors will require user
        intervention, unlike the issues corrected by the `_problem_correct`
        method. 
        ''' 
        # Check for a missing material -- Everything should have a MW
        # If not, this may suggest that there is a line that should not be
        # empty
        mw_mask = self.fulldata[MAT_MW].isna()
        if mw_mask.any():
            err_line = err_lines['miss_mw'] 
            df = self.fulldata[mw_mask].iloc[:,:3]
            raise CostError(err_line, df, self._disp_err_df)

        # Check for duplicated materials. This will probably be a big issue
        # with two materials sheets.
        dup_cpds = self.materials[RXN_CPD].duplicated()
        if dup_cpds.any():
            err_line = err_lines['dup_cpd'] 
            df = self.materials.loc[dup_cpds,].iloc[:, :3]
            raise CostError(err_line, df, self._disp_err_df)

        # Check for a missing cost, which is not being calculated
        # This is a tricky error because the costing will run just fine with 
        # NaNs. Check for this by looking for NaN in both Cost and Calc columns
        cost_mask = (self.fulldata[RXN_CST].isna() & \
                    self.fulldata[MAT_CST].isna())
        if cost_mask.any():
            err_line = err_lines['mis_cst']
            df = self.fulldata.loc[cost_mask, [MAT_CST, RXN_CST]]
            raise CostError(err_line, df, self._disp_err_df)
        
        # Check for duplicated materials in a single reaction.
        # When you select a single value from a reaction, you'll get a series
        # and not a float, e.g.
        dup_rxn = self.fulldata.loc[(self._fp_idx, self.final_prod), MAT_MW]
        if isinstance(dup_rxn, pd.Series):
            err_line = err_lines['dup_rxn'] 
            df = []
            gb = self.fulldata.groupby([RXN_STP, RXN_CPD])
            for prod, group in gb:
                if group.shape[0] > 1:
                    df.append(group)
            df = pd.concat(df).iloc[:, :3]
            raise CostError(err_line, df, self._disp_err_df)
            
        # Check that all the solvent information is given
        sol_cols = [MAT_DEN, RXN_VOL, RXN_REL, RXN_RCY]
        sol_mask = ~self.fulldata[RXN_VOL].isna() 
        sol_chk = self.fulldata.loc[sol_mask, RXN_REL].isna() | \
                self.fulldata.loc[sol_mask, MAT_DEN].isna() | \
                self.fulldata.loc[sol_mask, RXN_RCY].isna()
        # If anything is missing, print a note
        if sol_chk.any():
            err_line = err_lines['mis_sol']
            sol_mask2 = sol_mask & sol_chk
            df = self.fulldata.loc[sol_mask2, sol_cols]
            raise CostError(err_line, df, self._disp_err_df)
        
        # Check to make sure that the "Relative" compound for a solvent
        # is acutally contained in the step
        sol_mask = ~self.fulldata[RXN_REL].isna()
        bad = []
        for cpd in self.fulldata[sol_mask].itertuples():
            new_idx = (cpd.Index[0], getattr(cpd, RXN_REL) )
            if new_idx not in self.fulldata.index:
                bad.append( cpd.Index )

        if len(bad) > 0:
            err_line = err_lines['mis_rel']
            df = self.fulldata.loc[bad, sol_cols]
            raise CostError(err_line, df, self._disp_err_df)
        
    def calc_cost(self, excel=False):
        '''Calculate the cost of the route. 

        This is the main costing method. This will fill in the `fulldata`
        DataFrame with the appropriate values, and it will set a new class
        attribute `cost`, which is the total RM cost of the final product.

        Parameters
        ----------
        excel : boolean (False)
            If True, the calculations will output Excel-formated equation
            strings. By default (False), the output will be numerical values.
        '''
        # Save a time stamp so it can be displayed later
        self._now = pd.Timestamp.now('US/Eastern').strftime('%Y-%m-%d %H:%M')
        # Prep the DataFrame
        self._column_clear(excel=excel)
        # Set a column of unique row labels
        nrows = self.fulldata.shape[0]
        nrowdig = len(str(nrows))
        if excel:
            self.fulldata[RNUM] = [f'r{i:0{nrowdig:d}d}' for i in range(nrows)]
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
            The name of the product for a particular reaction.

        step : str
            The step ID for which to calculate the cost. Because the step
            labels can be numbers (e.g. step 1) and/or letters (e.g. step 1a),
            this parameter is handled as a string.

        amp : float, optional (default = 1.0)
            This number is an amplifier that increases some of the values,
            e.g. masses of materials, based on how much material is being used
            to make 1 kg of the overall final product.   

        excel : boolean (False)
            Whether to output numerical values (False) or Excel-formatted
            equation strings (True).
        '''
        # Select out the reaction of interest from the full data set. Saves
        # some typing. (Make this a copy to enusre it isn't given as a view.)
        data = self.fulldata.loc[step].copy()

        # Kg of nonsolvent materials used per equivalent
        data[RXN_KG] = data[RXN_EQ]*data[MAT_MW]
        # And for Excel
        if excel:
            data[DYN_RKG] = '=' + ECOLS[RXN_EQ] + data[RNUM] + '*'\
                    + ECOLS[MAT_MW] + data[RNUM]
        # For Excel, normalize the data here. This will get overwritten for
        # solvents. This is a little weird for the product because it will be
        # the same in numerator and denominator, but that will keep things
        # fully interactive.
        if excel:
            data[DYN_RKG] += '/(' + ECOLS[RXN_EQ] +\
                    data.loc[prod,RNUM] + '*' + ECOLS[MAT_MW] +\
                    data.loc[prod, RNUM] + ')'

        # Amount of solvent
        # First figure out which materials are solvents 
        mask = ~data[RXN_VOL].isna()
        sols = data[mask].copy()
        if sols.size != 0:
            # Find the material that the volumes are relative to (taking only
            # the first one... This might not be great.
            cpd_rel = sols[RXN_REL][0]
            # What is the kg of relative cpd?
            amt_rel = data.loc[cpd_rel, RXN_KG]
            if excel:
                amt_rel_e = data.loc[cpd_rel, RNUM]
            # Calculate the mass of solvent. Take into account the solvent
            # recycyling 
            # kg sol = Volume*Density*(1-Recycle)*(kg SM)
            sols[RXN_KG] = sols[RXN_VOL]*sols[MAT_DEN]*\
                    (1 - sols[RXN_RCY])*amt_rel
            # And for Excel
            if excel:
                sols[DYN_RKG] = '=' + ECOLS[RXN_VOL] +\
                        sols[RNUM] + '*' + ECOLS[MAT_DEN] +\
                        sols[RNUM] + '*' + '(1 - ' + ECOLS[RXN_RCY]\
                        + sols[RNUM] + ')*' + ECOLS[RXN_KG] +\
                        amt_rel_e
            data[mask] = sols

        # Normalize the kg of reaction
        data[RXN_KG] /= data.loc[prod, RXN_KG]
        # Remove the 1 kg of product
        data.loc[prod, RXN_KG] = np.nan
        if excel:
            data.loc[prod, DYN_RKG] = ''

        # Calculate unknown costs. Looks for any empty values in the "Cost" 
        # column. Don't use the RXN_CST column directly, because some of
        # the costs may have been manually set using the `value_mod` method
        # unknown_cost = ~data[RXN_CST].isna()
        unknown_cost = data[MAT_CST].isna()
        # This is recursive. The final cost per kg of product will be
        # amplified by each subsequent step, which is where the new_amp
        # calculation comes into play
        for cpd, row in data.loc[unknown_cost].iterrows():
            # Don't do this for the product of the current reaction
            if cpd == prod: 
                continue
            # The amounts needed will be amplified by the appropriate kg ratio.
            # Set that ratio
            new_amp = data.loc[cpd, RXN_KG]
            # And for dynamic excel output
            if excel:
                new_amp_enum = data.loc[cpd, RNUM]
                kg_col = ECOLS[RXN_KG]
                # This will be '*(C#)'. 
                new_eamp = f'*({kg_col}{new_amp_enum})'

            # The new step is needed as well
            new_stp = data.loc[cpd, RXN_CST]
            # Run the cost calculation for the unknown compound
            if excel:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp, 
                                 eamp + new_eamp, excel=excel)
            else:
                cst = self._rxn_cost(cpd, new_stp, amp*new_amp)
            # Set the calculated cost
            data.loc[cpd, MAT_CST] = cst
            # And for Excel -- This cost will need to get swapped out later.
            # Also need to check if an OPEX is necessary
            if excel:
                if np.isnan(self.fulldata.loc[(new_stp, cpd), RXN_OPX]):
                    data.loc[cpd, DYN_CST] = '=' + ECOLS[MAT_CST] +\
                            self.fulldata.loc[(new_stp, cpd), RNUM]
                else:
                    data.loc[cpd, DYN_CST] = '=' + ECOLS[MAT_CST] +\
                            self.fulldata.loc[(new_stp, cpd), RNUM] + '+'\
                            + ECOLS[RXN_OPX] +\
                            self.fulldata.loc[(new_stp, cpd), RNUM]

        # Calculate the cost for each material in the reaction
        data[RXN_RMC] = data[RXN_KG]*data[MAT_CST]
        # The product cost will be the sum of all the reactant/solvent costs
        data.loc[prod, RXN_RMC] = data[RXN_RMC].sum()
        # Set the "Cost" to the calculated value
        data.loc[prod, MAT_CST] = data.loc[prod, RXN_RMC]
        # And for Excel
        if excel:
            data[DYN_RRMC] = '=' + ECOLS[RXN_KG] + data[RNUM] + '*' + \
                            ECOLS[MAT_CST] + data[RNUM]
            # And for Excel, first we need the rows that are not the product
            mask = data.index != prod
            # Then we need to make a comma-separated list of these rows
            cells = [f'{ECOLS[RXN_RMC]}{r}' for r in data.loc[mask, RNUM]]
            rs = ','.join(cells)
            # Combine them together into a sum
            data.loc[prod, DYN_RRMC] = '=SUM(' + rs + ')'
            # Set the cost to the calculated value
            data.loc[prod, DYN_CST] = '=' + ECOLS[RXN_RMC] +\
                    data.loc[prod, RNUM] 

        # Calculate % costs for individual rxn
        # = (RM cost/kg rxn)/(RM cost/kg rxn for the rxn product)
        p_rm_cost = data[RXN_RMC]*100/data.loc[prod, RXN_RMC]
        data[RXN_RMP] = p_rm_cost.values
        # And for Excel
        if excel:
            data[DYN_RRMP] = '=' + ECOLS[RXN_RMC] +\
                    data[RNUM] + '*100/' + ECOLS[RXN_RMC] +\
                    data.loc[prod, RNUM]
        # Remove the % cost for the rxn product
        data.loc[prod, RXN_RMP] = np.nan
        # And for Excel
        if excel:
            data.loc[prod, DYN_RRMP] = ''

        # These are the costs for ultimate product
        # For one reaction amp=1, so the individual rxn cost = ultimate rxn 
        # cost. However, for feeder reactions this will get amplified by 
        # each step   
        data[PRD_RMC] = data[RXN_RMC]*amp
        # And for Excel
        if excel:
            data[DYN_PRMC] = '=' + ECOLS[RXN_RMC] + data[RNUM] + eamp
        
        # This sets the number of kg of each material per kilogram of product
        # This is done by multiplying the per reaction value by the amplifier
        # This isn't necessary for costing, but we can use it for PMI
        data[PRD_KG] = data[RXN_KG]*amp
        # And for Excel
        if excel:
            data[DYN_PKG] = '=' + ECOLS[RXN_KG] + data[RNUM] + eamp
        
        # Set the values in the big DataFrame with this slice. This goofy call
        # is necessary to make sure the data are set correctly. 
        self.fulldata.loc[(step, data.index), data.columns] = data.values
        
        # Return the calculated product cost, which is required for the 
        # recursive nature of the algorithm. In addition, an optional OPEX
        # may be added to take into account production costs of the cpd
        if np.isnan(data.loc[prod, RXN_OPX]):
            return data.loc[prod, RXN_RMC]
        else:
            return data.loc[prod, RXN_RMC] + data.loc[prod, RXN_OPX]
    
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
        opex = self.fulldata.loc[(step, prod), RXN_OPX]
        if not np.isnan(opex):
            self.fulldata.loc[(step, prod), MAT_CST] = self.cost
            # And for Excel
            if excel:
                self.fulldata.loc[(step, prod), DYN_CST] += '+' + \
                        ECOLS[RXN_OPX] + \
                        self.fulldata.loc[(step, prod), RNUM]
                
        # Calculate % overall costs relative to the prod
        self.fulldata[PRD_RMP] = self.fulldata[PRD_RMC]*100/self.cost

        # And for Excel
        if excel:
            self.fulldata[DYN_PRMP] = '=' +\
                    ECOLS[PRD_RMC] + self.fulldata[RNUM] +\
                    '*100/' + ECOLS[RXN_RMC] +\
                    self.fulldata.loc[(step, prod), RNUM]
        
        # Filter out certain values to simplify full data set
        # Remove the cost and %s for cost-calculated materials
        # This is necessary so that this column adds up to 100% (w/o OPEX)
        mask = ~self.fulldata[RXN_CST].isna()
        self.fulldata.loc[mask, PRD_RMP] = np.nan
        if excel:
            self.fulldata.loc[mask, DYN_PRMP] = ''

        # This filters some of the costs which are simply the sum of raw
        # materials from earlier rxns. The sum of this column will now be
        # equal to the cost of the final product.
        self.fulldata.loc[mask, PRD_RMC] = np.nan
        if excel:
            self.fulldata.loc[mask, DYN_PRMC] = ''

        # This filters out the kg/kg prod values that were calculated, so that
        # the sum of this column is the PMI
        self.fulldata.loc[mask, PRD_KG] = np.nan
        if excel:
            self.fulldata.loc[mask, DYN_PKG] = ''

        # PMI Calculations
        # Adding a prefix for display purposes
        # There will be a funny column in this DF with these values...
        self._pre = '*'
        
        # First of all, calculate the PMI for each reaction individually
        gb = self.fulldata.groupby(RXN_STP)
        if excel:
            rxn_pmi = gb.agg({RXN_KG:'sum', RNUM:self._excel_pmi})\
                            .rename({RNUM: DYN_RKG}, axis=1)\
                            .reset_index()
        else:
            rxn_pmi = gb.agg({RXN_KG:'sum'}).reset_index()

        rxn_pmi[RXN_CPD] = self._pre + 'Step ' + rxn_pmi[RXN_STP] + ' PMI'
        
        # The full route PMI is not the sum of the above, but is the sum of
        # the PRD_KG column. We need to make this into a DataFrame to
        # merge with the per reaction values above.
        df_vals = {PRD_KG: [self.fulldata[PRD_KG].sum()], 
                   RXN_STP: [self._fp_idx],
                   RXN_CPD: [self._pre*2 + 'Full Route PMI'],
                   }
        if excel:
            mask = ~self.fulldata[PRD_KG].isna()
            all_cells = [f'{ECOLS[PRD_KG]}{i}' for i in\
                         self.fulldata.loc[mask, RNUM]]
            cell_sum = '=SUM(' + ','.join(all_cells) + ')'
            df_vals[DYN_PKG] = [cell_sum]
        full_pmi = pd.DataFrame(df_vals)

        # Merge the per-reaction and full PMI
        self.pmi = pd.concat([rxn_pmi, full_pmi], sort=False)\
                .set_index(RXN_STP)

    def _excel_pmi(self, col):
        # A function for creating the dynamic Excel cells for PMI
        cells = [f'{ECOLS[RXN_KG]}{i}' for i in col]
        return "=SUM(" + ','.join(cells) + ")"

    def results(self, style='compact'):
        '''Returns a formatted/simplified full output DataFrame.

        The `fulldata` attribute contains a number of potentially unnecessary
        columns that need to be filtered out for presentation purposes. This
        also combines the PMI data as well, so everything is in one table.

        Parameters
        ----------
        style : str ('compact')
            By default, this returns a minimal representation of the system
            ('compact'); however, a more complete set of columns can be
            returned if the value `'full'` is passed as the parameter here.
        '''
        # For compact display, these are the most important columns
        comp_col = [MAT_CST, RXN_EQ,] 
        if self.fulldata.Volumes.any():
            comp_col.extend([RXN_VOL, RXN_RCY,])
        if self.fulldata.OPEX.any():
            comp_col.append(RXN_OPX) 
        comp_col.extend([RXN_KG, RXN_RMC, RXN_RMP, PRD_KG, PRD_RMC, PRD_RMP])

        # For full display, don't show the Excel columns
        ecol_mask = ~self.fulldata.columns.str.contains(SUFFIX+'|'+RNUM)
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
        # So one needs to be dropped.
        gb = fd.groupby(RXN_STP)
        concated = gb.apply(self.__fd_pmi_concat, pmi)
        concated = concated.drop(RXN_STP, axis=1)\
                            .reset_index()\
                            .drop('level_1', axis=1)
        
        # Reset the index and return the DF. No sorting is necessary here.
        return concated.set_index([RXN_STP, RXN_CPD])                

    def __fd_pmi_concat(self, df, pmi):
        '''Concats the fulldata and pmi DataFrames.

        This is only used in a DataFrame.apply call in the _df_combine
        method.'''
        # The step number
        step = df[RXN_STP].iloc[0]
        # Mask out the correct PMI values
        pmi_mask = pmi[RXN_STP] == step
        pmi_small = pmi[pmi_mask]
        # Combine the route step info with the PMI. In this way, the pmi will
        # be at the end.
        return pd.concat([df, pmi_small], ignore_index=True)

    def excel(self, ):
        '''Return an Excel-formatted DataFrame.

        This table will have equation strings rather than numerical values,
        and it can be saved to disk using the `DataFrame.to_excel` method, if
        desired. Because the table contains equation strings, the Excel file
        will be fully dynamic, in that value changes will cause a
        recalculation of the entire sheet.
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
        for r, n in zip(fd[RNUM], np.arange(2, nrows+2)):
            if isinstance(r, float):
                continue
            for col in fd:
                if SUFFIX in col:
                    fd[col] = fd[col].str.replace(r, str(n))
                    
        # Set dynamic cost for the "Cost" of calculated products
        mask = fd[DYN_CST] != ''
        fd.loc[mask, MAT_CST] = fd.loc[mask, DYN_CST]
        
        # Drop the non-dynamic columns
        d_cols = [DYN_CST, RXN_KG, RXN_RMC, RXN_RMP, PRD_KG, PRD_RMC, PRD_RMP,
                RNUM]
        fd = fd.drop(d_cols, axis=1)
        
        # Rename columns to remove "dyn" SUFFIX
        col_str = fd.columns.str.replace(SUFFIX, '')
        new_cols = {c:cn for (c, cn) in zip(fd.columns, col_str)}
        fd = fd.rename(new_cols, axis=1)

        # Reset all "empty" cells to NaN. This is important for unit testing
        # purposes.
        fd = fd.where(fd != '', np.nan) 
        
        # Rerun the cost calculation without the excel stuff to get rid of all
        # the other columns
        self.calc_cost(excel=False)

        return fd

    def value_mod(self, cpd, val, val_type='Cost', step=None):
        '''Manually set a value for a given material.

        This is useful for creating alternative models without having to
        import completely new reactions and materials tables. The original
        input tables are not modified, so the model can be reset with the
        `rxn_data_setup` method. 

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
        separate step. This behavior is by design so that several `value_mod`
        method calls can be made, if necessary, before the costs are
        recaclulated.
        '''
        # Store the values
        self._mod_vals.append( (cpd, val, val_type, step) )
        # This will clear out the old calculation data and set the modified
        # value. Keeps folks from getting confused b/c previously calculated
        # values appear unchanged.
        self._column_clear()

    def _set_val(self, cpd, val, val_type, step):
        '''Set a modified value in the `fulldata` DataFrame.

        The modified values are set using the `value_mod` method, which adds
        the values to a `_mod_vals` list. This function cycles through that
        list, changing the associated values.
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
            self.fulldata.loc[cells, RXN_CST] = np.nan


