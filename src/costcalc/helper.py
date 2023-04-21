from .constants import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def value_scan(model, cpd, start, stop, npts, val_type='$/kg', step=None):
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

    val_type : str (Default = '$/kg')
        This is the column name for the parameter that you'll be changing.
        This must be for a non-calculated column, such as '$/kg', 'Equiv',
        'OPEX', etc.

    step : None, int (Default = None)
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
        This DataFrame has 'Values' and '$/kg' columns for the input
        values and output costs, repectively. 

    Notes
    -----
    Although this method recalculates the cost for every value, it does
    not modify the original `fulldata` or `cost` attributes. 
    '''
    # Create the values array
    vals = np.linspace(start, stop, npts)
   
    # I need a copy of the full data set in order to reset for each
    # iteration. Using this rather than rxn_data_setup because this makes it
    # possible to use a modified model with this function
    fd_copy = model.fulldata.copy()
    all_costs = []
    for val in vals:
        model.value_mod(cpd, val, val_type, step)
        model.calc_cost()
        all_costs.append(model.cost)
        # Reset the full data set 
        model.fulldata = fd_copy.copy()
        # Pop out the mod value, otherwise this list will get really long
        model._mod_vals.pop()

    # Reset the final cost
    model.cost = model.fulldata.loc[(model._fp_idx, model.final_prod), 
                              'RM cost/kg rxn']
    
    return pd.DataFrame({'Values':vals, '$/kg':all_costs})


def plot_scan(model, cpd, start, stop, npts, val_type='$/kg', step=None, 
            label=None):
    '''Plot a range of values.

    Parameters
    ----------
    See `value_scan` method description, except for the following.

    label : None, bool, str (Default = None)
        This will add a legend label. If you use the value `True`, the
        compound name will be added to the legend. Otherwise, you can pass a
        custom string if you want that to be in the legend instead. If None
        (default) or False is used, then no legend entry will be given. This
        function does add the legend to the plot; use pyplot.legend as a
        separate function call.
    '''
    # Calculate all the costs for the given values
    costs = value_scan(model, cpd, start, stop, npts, val_type=val_type,
                        step=step)
    
    # If label is boolean, make sure the label is set properly
    if label == True:
        label = cpd
    elif label == False:
        label = None

    # Plot the values
    plt.plot(costs['Values'], costs['$/kg'], 'o', label=label)


def sensitivity(model, col='Equiv', frac=0.1, decimals=2):
    '''Do a sensitivity analysis for the equivalents of reagents.

    Parameters
    ----------
    model : CoreCost or subclass
        This is the costing model instance that will be modified.
        
    col : str (Default = 'Equiv')
        Which column from the `fulldata` DataFrame should be used for the
        sensitivity analysis. 

    frac : float (Default = 0.1)
        Fractional percentage to increase/decrease the values by before
        recosting. The default is 0.1, which is the same as +/- 10%. 

    decimals : int or None (Default = 2)
        How many decimal places to display. If `None`, the full precision
        DataFrame will be displayed.
    '''
    # Make a new DF for sensitivity analysis. If material costs are being
    # screened, then remove the compounds being costed
    sens = model.fulldata[[col]].dropna()
    if col == MAT_CST:
        mask = model.fulldata[RXN_CST].isna()
        sens = sens[mask]
    # Make values that are a certain percent above and below the current
    # numbers
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
        model.value_mod(cpd, vals['Val low'], val_type=col, step=step)
        model.calc_cost()
        cost_low = model.cost
        cost_low_per = (cost_low*100./cost_save) - 100
        # Just to be safe - Drop the modified value
        model._mod_vals.pop()

        # High values
        model.value_mod(cpd, vals['Val high'], val_type=col, step=step)
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


