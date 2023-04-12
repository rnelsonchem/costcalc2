from .constants import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def value_mod(model, cpd, val, val_type='Cost', step=None):
    '''Manually set a value for a given material.

    Parameters
    ----------
    model : CoreCost or subclass
        This is the costing model instance that will be modified.
        
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
    model._mod_vals.append( (cpd, val, val_type, step) )
    # This will clear out the old calculation data and set the modified
    # value. Keeps folks from getting confused b/c calculated values are
    # unchanged.
    model._column_clear()


def value_scan(model, cpd, start, stop, npts, val_type='Cost', step=None):
    '''Scan a range of values for a given material.
    
    Parameters
    ----------
    model : CoreCost or subclass
        This is the costing model instance that will be modified.
        
    cpd : str
        This the compound name for which the value will be modified.

    start : int, float
        The starting value for the scan.

    stop : int, float
        The ending value for the scan.

    npts : int
        The numbers of points to calculate the costs between `start` and
        `stop`

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
    fd_copy = model.fulldata.copy()
    all_costs = []
    for val in vals:
        value_mod(model, cpd, val, val_type, step)
        model.calc_cost()
        all_costs.append(model.cost)
        # Reset the full data set 
        model.fulldata = fd_copy.copy()
        # Pop out the mod value, otherwise this list will get really long
        model._mod_vals.pop()

    # Reset the final cost
    model.cost = model.fulldata.loc[(model._fp_idx, model.final_prod), 
                              'RM cost/kg rxn']
    
    return pd.DataFrame({'Values':vals, 'Costs':all_costs})


def plot_scan(model, cpd, start, stop, npts, val_type='Cost', step=None, 
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
    costs = value_scan(model, cpd, start, stop, npts, val_type=val_type,
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


def swap(model, cpd_old, cpd_new, step=None):
    '''Swap one compound for another in the route.

    Parameters
    ----------
    model : CoreCost or subclass
        This is the costing model instance that will be modified.
        
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
        mat_mask = model.materials.Compound == cpd_new
        # Throw an error if it is not there
        if not mat_mask.any():
            raise ValueError("Oops! Your new compound isn't in the"
                            " materials list.")
        # There should only be one entry, so we'll select that one
        mat_vals = model.materials[mat_mask].iloc[0]
        # Set some values that will be used for updating
        cpd_name = cpd_new
        mw = mat_vals[MAT_MW]
        density = mat_vals[MAT_DEN]
        cost = mat_vals[MAT_CST]
    
    # If the new compound is a costing instance...
    elif isinstance(cpd_new, (ExcelCost, ColabCost)):
        # The value selection is a little different
        cpd_name = cpd_new.final_prod
        cpd_loc = (cpd_new._fp_idx, cpd_new.final_prod)
        mw = cpd_new.fulldata.loc[cpd_loc, MAT_MW]
        density = cpd_new.fulldata.loc[cpd_loc, MAT_DEN]
        cost = cpd_new.fulldata.loc[cpd_loc, MAT_CST]
    
    # Else you've added the wrong kind of new compound
    else:
        raise ValueError("Oops!! Your new compound is not "
                        "of the correct type.")

    # Process the fulldata array
    # First, remove the MultiIndex
    fd_rst = model.fulldata.reset_index()
    # Check for the requested compound
    cpd_mask = fd_rst[RXN_CPD] == cpd_old
    # If a particular step is chosen, select only that step
    if step:
        cpd_mask = cpd_mask & (fd_rst[RXN_STP] == str(step))

    if not cpd_mask.any():
        raise ValueError("Oops! '" + cpd_old + "' isn't in your "
                        "current route.")
    # Swap out the compound names
    fd_rst.loc[cpd_mask, RXN_CPD] = cpd_name
    # Reset index and fulldata attribute
    model.fulldata = fd_rst.set_index([RXN_STP, RXN_CPD])

    # Set the values using `value_mod`
    model.value_mod(cpd_name, mw, val_type='MW', step=step)
    model.value_mod(cpd_name, density, val_type='Density', step=step)
    model.value_mod(cpd_name, cost, step=step)


def sensitivity(model, col='Equiv', frac=0.1, decimals=2):
    '''Do a sensitivity analysis for the equivalents of reagents.

    Parameters
    ----------
    model : CoreCost or subclass
        This is the costing model instance that will be modified.
        
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
    sens = model.fulldata[[col]].dropna()
    sens['Val low'] = sens[col]*(1 - frac)
    sens['Val high'] = sens[col]*(1 + frac)
    
    # Re-run the costing under the current conditions, which resets the
    # cost and fulldata variables. 
    model.calc_cost()
    # Make copies of these values so they don't change
    cost_save = model.cost
    fd_save = model.fulldata.copy()

    # Loop through the values and calculate the cost if these values
    # increase or decease by the percent given
    for step_cpd, vals in sens.iterrows():
        step, cpd = step_cpd
        # Low values
        value_mod(model, cpd, vals['Val low'], val_type=col, step=step)
        model.calc_cost()
        cost_low = model.cost
        cost_low_per = (cost_low*100./cost_save) - 100
        # Just to be safe - Drop the modified value
        model._mod_vals.pop()

        # High values
        value_mod(model, cpd, vals['Val high'], val_type=col, step=step)
        model.calc_cost()
        cost_high = model.cost
        cost_high_per = (cost_high*100./cost_save) - 100
        model._mod_vals.pop()

        # Set the values in the sensitivity DataFrame
        sens.loc[(step, cpd), 'Cost low'] = cost_low
        sens.loc[(step, cpd), 'Cost high'] = cost_high
        sens.loc[(step, cpd), '% low'] = cost_low_per
        sens.loc[(step, cpd), '% high'] = cost_high_per

        # Reset the original values.  This is necessary so the next step
        # in the loop starts from the original position.
        model.cost = cost_save
        model.fulldata = fd_save.copy()
    
    if decimals:
        return sens.round(decimals)
    else:
        return sens


