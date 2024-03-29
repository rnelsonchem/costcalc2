from .constants import *

# Setting the display function
# This gets changed for Jupyter Notebook/IPython sessions, so that DataFrames
# are displayed in a fancier format
try:
    from IPython.display import display as disp
except:
    disp = print

class CostError(Exception):
    '''A custom exception for catching common input errors.

    In addition to error-type specific messages, this also captures and
    optionally displays the offending portions of the input DataFrames. 

    Parameters
    ----------
    message : str
        The error message to display. See `err_lines` dictionary for the
        custom messages.

    df : DataFrame (None)
        The bad portion of the input DataFrame.

    disp_df : Boolean (False)
        A flag that determines whether the bad DataFrame information is
        displayed along with the error message. Although this is useful for
        debugging, this parameter is set to False by default, because this
        could print sensitive information to the error logs.


    Notes
    -----
    See the `err_lines` dictionary defined below for the text of the various
    error messages. See the `CoreCost._sanity_check` method for the checks
    that are run.
    '''
    def __init__(self, message, df=None, disp_df=False):
        self.message = message
        self.df = df
        self.disp_df = disp_df
        super().__init__(message)

    def __str__(self, ):
        if self.disp_df:
            self.message += '\nSee the following entries:\n'
            self.message += self.df.to_string()
        return self.message


# A dictionary containing custom error messages for common input errors
err_lines = {
    'miss_mw': 'Missing MW error! \n'\
            'This most commonly happens for 1 of 2 reasons:\n'\
            '1. A reaction compound is not defined in the '\
            'materials table.\n'\
            '2. The compound names in reaction/materials '\
            'tables do not match *exactly* (even white space).',

    'dup_cpd': 'Duplicated material error!\n'\
            'Check that you do not have two or more compounds of the same\n'\
            'name in your material sheet(s).',

    'mis_cst': 'Missing material cost error!\n'\
           'This may be due to a missing material and/or its cost.\n',

    'mis_sol': 'Missing solvent info error!\n'\
           'One of your solvent entries is missing critical info.',
           
    }
    
