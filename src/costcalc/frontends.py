import urllib.parse as parse
import warnings

import numpy as np
import pandas as pd

from .core import CoreCost
from .constants import *

# Filter warnings about data validation in Excel files
warnings.filterwarnings('ignore', category=UserWarning, 
                module='openpyxl')

# Setting the display function
# This gets changed for Jupyter Notebook/IPython sessions, so that DataFrames
# are displayed in a fancier format
try:
    from IPython.display import display as disp
    from IPython.display import Javascript
except:
    disp = print

# Set pandas to display lots of DataFrame rows so things don't get cut out
pd.options.display.max_rows = 1000
# For future reference -- Can't use this in conjunction with "fillna" becuase
# the DF doesn't display correctly 
# Set Pandas precision
#pd.set_option('precision', 2)

class ExcelCost(CoreCost):
    '''Costing class designed for reading data from Excel/csv spreadsheets.

    The Excel/csv files require specific columns to be present. See the
    `CoreCost` class doc string for more information.

    Parameters
    ----------
    materials_file : str
        The name/path of the materials list file. Can be xlsx or csv.

    rxn_file : str
        The name/path of the file defining the reactions. Can be xlsx or csv.
        Both `rxn_file` and `materials_file` can be the same file.

    final_prod : str
        The name of the final product for route. 

    materials_sheet : int, str (default = 0)
        The sheet in the Excel file that contains the materials information.
        The default (0) is the first sheet in the spreadsheet. You could use a
        different number or a name (str) if you want a different sheet from the
        spreadsheet.

    rxn_sheet : int, str (default = 0)
        See `materials_sheet` description, except this is for the sheet
        defining the reactions.

    alt_mat_file : str (default = None)
        The name/path of an optional, secondary materials sheet. This is
        useful if you have separate master and route-specific materials
        sheets, for example. Can be an xlsx or csv.

    alt_mat_sheet : int, str (default = 0)
        The sheet number/name for the secondary materials sheet. See
        `materials_sheet` description. 

    dis_err_df : boolean (default = False)
        If a `CostError` exception is raised, this controls whether the
        offending section of the DataFrame will be printed along with the
        error. This is useful for debugging, but is typically set to False, so
        that potentially sensitive information is not printed to the error
        logs.

    Notes
    -----
    See the `CoreCost` class doc string for additional information. 
    '''
    def __init__(self, materials_file, rxn_file, final_prod,
            materials_sheet=0, rxn_sheet=0, 
            alt_mat_file=None, alt_mat_sheet=0,
            disp_err_df = True):
        # Set up the reaction DataFrame
        self._rxn_file = rxn_file
        self._rxn_sheet = rxn_sheet
        rxns = self._rxn_read()

        # Create the Materials DataFrame from a main sheet and an optional
        # alternate sheet.
        self._materials_file = materials_file
        self._materials_sheet = materials_sheet
        self._alt_mat_file = alt_mat_file
        self._alt_mat_sheet = alt_mat_sheet
        materials = self._materials_build()                

        # Run the __init__ method from the CoreCost class
        super(ExcelCost, self).__init__(materials, rxns, final_prod,
                disp_err_df)

    def _rxn_read(self, ):
        '''Read an Excel sheet that defines the reactions.
        '''
        rxns = self._get_sheet_vals(self._rxn_file, self._rxn_sheet,
                            dtypes={RXN_STP:str,})
        return rxns

    def _get_sheet_vals(self, fname, fsheet, dtypes=None):
        '''A simple Excel/CSV reader function for both reaction and materials
        files.
        '''
        # Read the file, drop NaN-only and commented rows.
        if isinstance(fname, str) and fname[-3:].lower() == 'csv':
            df = pd.read_csv(fname, dtype=dtypes, comment='#')\
                            .dropna(how='all')
        else:
            df = pd.read_excel(fname, fsheet, dtype=dtypes, comment='#')\
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
        return materials

    def _materials_read(self, mat_file, wsheet):
        '''Read an Excel file sheet that defines the materials used in
        costing.

        This is a separate method so it can be overwritten in other classes.
        '''
        mats = self._get_sheet_vals(mat_file, wsheet,)
        return mats

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

        fill : str or None, optional (Default = '-')
            Fill NaN values with this string. This makes the table a little
            easier to read. Set this to `None` if you want to see the table
            with the typical NaN labels.
        '''
        # Get the prepared DataFrame for printing 
        fd = super(ExcelCost, self).results(style)

        # Print the time the calculation was run
        utc_datetime = self._now.strftime('%Y-%m-%d %H:%M UTC')
        print(f'As of {utc_datetime} --')
        
        # Print a string about the final cost of the product
        cost_str = f'The final cost of {self.final_prod} is '\
                f'${self.cost:.2f}/kg.'
        print(cost_str)
            
        # Display the correct format of data based on the kwargs
        if isinstance(decimals, int):
            disp(fd.round(decimals).fillna(fill))
        else:
            disp(fd.fillna(fill))

    def excel(self, fname, ):
        '''Save the costing DataFrame as a dynamic Excel file.

        Parameters
        ----------
        fname : str, ExcelWriter
            The name you want to give to the Excel file, or an ExcelWriter
            object to send multiple costing sheets to the same file. See
            `pandas.DataFrame.to_excel` documentation for details.

        '''
        # Get the process DataFrame
        fd = super(ExcelCost, self).excel()
       
        # Can't use a ":" in an Excel sheet name for time
        utc_datetime = self._now.strftime('%Y-%m-%d %Hh%Mm UTC')
        # Create the excel file. Can only save with the date and not the time
        fd.to_excel(fname, sheet_name=f'As of {utc_datetime}')
            

class ColabCost(ExcelCost):
    '''Costing class designed for the Colab Python environment.

    The materials and reaction information must be in Google Sheet format; the
    urls below are the URL values for those sheets. The user must have read
    access for the sheets to be able to use them. 

    Parameters
    ----------
    materials_url : str
        The URL to the Google sheet containing the materials list. 

    rxn_url : str
        The URL to the Google sheet containing the reaction information. 

    alt_mat_url : str, optional (default = None)
        The URL to the Google sheet containing an optional, secondary
        materials list.  This is useful if you have separate master and user
        materials sheets, for example.

    Notes
    -----
    This is subclass of `ExcelCost` and `CoreCost`, so see those class
    docstrings for more information.

    '''
    def __init__(self, materials_url, rxn_url, final_prod, materials_sheet=0,
            rxn_sheet=0, alt_mat_url=None, alt_mat_sheet=0):
        # Do some imports that are only possible in the Colab environment
        # This should prevent these from running in a non-Colab environment
        from google.auth import default
        from google.colab import auth
        from google.colab import files
        import gspread
        # These will have to be made global
        global default
        global auth
        global files
        global gspread

        # Authenticate the Colab environment 
        auth.authenticate_user()
        creds, _ = default()
        self._gc = gspread.authorize(creds)
        
        # Fix the final product and setup a mod variable
        super(ColabCost, self).__init__(materials_url, rxn_url, final_prod,
                materials_sheet, rxn_sheet, alt_mat_url, alt_mat_sheet)
        
    def _materials_read(self, mat_url, wsheet):
        '''Read a Google sheet that defines the materials used in costing.

        Parameters
        ----------
        mat_key : str
            The unique key for a Google spreadsheet that defines the
            materials. 

        wsheet : str or int
            The specific sheet to extract from the Google spreadsheet. 
        '''
        mats = self._get_sheet_vals(mat_url, wsheet)

        # Warning for bad column namen
        self._cost_warn(mats)

        # Convert numeric/date columns. Everything is read from a Google sheet
        # as strings
        num_cols = [MAT_MW, MAT_DEN, MAT_CST] 
        for nc in num_cols:
            mats[nc] = pd.to_numeric(mats[nc])

        return mats
        
    def _rxn_read(self, ):
        '''Read a Google Sheet containing reaction info and convert it to a
        DataFrame.
        '''
        rxns = self._get_sheet_vals(self._rxn_file,
                                     self._rxn_sheet)

        # Warning for bad column names
        self._col_warn(rxns)

        # Set some rxns columns to numeric values. Everything is read from a
        # Google sheet as strings
        num_cols = [RXN_EQ, RXN_VOL, RXN_RCY, RXN_OPX, RXN_MS]
        for nc in num_cols:
            if nc in rxns.columns:
                rxns[nc] = pd.to_numeric(rxns[nc])
        
        return rxns
        
    def _get_sheet_vals(self, key, sheet):
        '''General code for obtaining Google Sheet values and returning a 
        DataFrame.
        '''
        # Check if the `key` is a URL. If so, pull apart the url and grab the
        # part of the path that is the key to the spreadsheet.
        if key.startswith('http'):
            url_frag = parse.urlparse(key)
            key = url_frag.path.split('/')[3]

        # Using the Google sheet handle, pull down all values and make a 
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
        # Drop empty rows, these are particularly difficult to error check
        # sometimes. Assuming that both materials/rxns sheets should have
        # compound names for valid lines.
        cpd_mask = val_df[RXN_CPD].isna()
        val_df = val_df[~cpd_mask]

        # Drop lines with comment markers
        mask = ~val_df.iloc[:, 0].str.startswith('#')
        val_df = val_df[mask]
        
        return val_df

    def results(self, style='compact', decimals=2, fill='-'):
        '''Print the results of the costing calculation.

        See the `ExcelCost.results` method doc string for details on the
        function parameters.
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

    def excel(self, fname):
        '''Download the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str
            The name you want to give to the Excel file.

        '''
        super(ColabCost, self).excel(fname, )

        # Use the Colab function for downloading the file from the server
        files.download(fname)
        

class WebAppCost(ExcelCost):
    '''Costing class designed for the Streamlit Web App environment.

    This is largely the same as the `ExcelCost` class, so see that method for
    a more complete docstring. 
    '''
    def __init__(self, materials_file, rxn_file, final_prod,
            materials_sheet=0, rxn_sheet=0, 
            alt_mat_file=None, alt_mat_sheet=0,
            disp_err_df=False):

        super(WebAppCost, self).__init__(materials_file, rxn_file, final_prod,
                materials_sheet, rxn_sheet, alt_mat_file, alt_mat_sheet,
                disp_err_df)

    def results(self, style='compact', decimals=2, fill=np.nan):
        '''Returns the results of the costing calculation.

        See the `ExcelCost.results` method for additional documentation. This
        is largely the same as that method except for the optional `fill`
        keyword argument below.

        Parameters
        ----------
        fill : str, None, or Numpy.nan (Default = Numpy.nan)
            Fill NaN values with this string. This makes the table a little
            easier to read. Set this to `None` if you want to see the table
            with the typical NaN labels.
        '''
        # Use the CoreCost.results method to get the initial results DataFrame
        fd = super(ExcelCost, self).results(style)

        # Display the DataFrames for different permutations of kwargs
        if isinstance(decimals, int):
            fd = fd.round(decimals).fillna(fill)
        else:
            fd = fd.fillna(fill)

        return fd


