# Import everything from the core file to get all of the column variable
# definitions
from .core import *

import matplotlib.pyplot as plt

# Set up some plotting stuff for the notebooks
plt.style.use('ggplot')
plt.rc('figure', dpi=150)


class HelperFuncs(object):
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
            self.fulldata.loc[cells, rxn_cst] = np.nan

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
            mw = mat_vals[mat_mw]
            density = mat_vals[mat_den]
            cost = mat_vals[mat_cst]
        
        # If the new compound is a costing instance...
        elif isinstance(cpd_new, (ExcelCost, ColabCost)):
            # The value selection is a little different
            cpd_name = cpd_new.final_prod
            cpd_loc = (cpd_new._fp_idx, cpd_new.final_prod)
            mw = cpd_new.fulldata.loc[cpd_loc, mat_mw]
            density = cpd_new.fulldata.loc[cpd_loc, mat_den]
            cost = cpd_new.fulldata.loc[cpd_loc, mat_cst]
        
        # Else you've added the wrong kind of new compound
        else:
            raise ValueError("Oops!! Your new compound is not "
                            "of the correct type.")

        # Process the fulldata array
        # First, remove the MultiIndex
        fd_rst = self.fulldata.reset_index()
        # Check for the requested compound
        cpd_mask = fd_rst[rxn_cpd] == cpd_old
        # If a particular step is chosen, select only that step
        if step:
            cpd_mask = cpd_mask & (fd_rst[rxn_stp] == str(step))

        if not cpd_mask.any():
            raise ValueError("Oops! '" + cpd_old + "' isn't in your "
                            "current route.")
        # Swap out the compound names
        fd_rst.loc[cpd_mask, rxn_cpd] = cpd_name
        # Reset index and fulldata attribute
        self.fulldata = fd_rst.set_index([rxn_stp, rxn_cpd])

        # Set the values using `value_mod`
        self.value_mod(cpd_name, mw, val_type='MW', step=step)
        self.value_mod(cpd_name, density, val_type='Density', step=step)
        self.value_mod(cpd_name, cost, step=step)

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


