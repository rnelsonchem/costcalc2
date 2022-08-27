import urllib.parse as parse
from io import BytesIO

from .core import CoreCost

import numpy as np
import pandas as pd

# Set pandas to display lots of DataFrame rows so things don't get cut out
pd.options.display.max_rows = 1000
# For future reference -- Can't use this in conjunction with "fillna" becuase
# the DF doesn't display correctly 
# Set Pandas precision
#pd.set_option('precision', 2)


class ExcelCost(CoreCost):
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
        super(ExcelCost, self).__init__(rxns, materials, final_prod)

    def _rxn_read(self, ):
        '''Read an Excel sheet that defines the reactions.
        '''
        rxns = self._excel_csv_reader(self._rxn_file, self._rxn_sheet,
                            dtypes={'Step':str, 'Cost calc':str})
        return rxns

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
        # Drop lines that are still empty, these cause all sorts of problems
        # Assume rxn/materials sheets should have Cpd names for valid entries
        cpd_mask = df['Compound'].isna()
        df = df[~cpd_mask]
        
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
        mats = self._excel_csv_reader(mat_file, wsheet,)
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
        print('As of', self._now, '--')
        
        # Print a string about the final cost of the product
        if decimals:
            dec_str = ':.{:d}f'.format(decimals)
        else:
            dec_str = ':f'
        cost_str = f'The final cost of {self.final_prod} is '\
                f'${self.cost:.2f}/kg.'
        print(cost_str)
            
        # Display the correct format of data based on the kwargs
        if decimals:
            print(fd.round(decimals).fillna(fill))
        else:
            print(fd.fillna(fill))

    def excel(self, fname, ):
        '''Save the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str, ExcelWriter
            The name you want to give to the Excel file, or an ExcelWriter
            object to send multiple to the same file. See
            `pandas.DataFrame.to_excel` documentation for details.

        '''
        # Get the process DataFrame
        fd = super(ExcelCost, self).excel()

        # Create the excel file. Can only save with the date and not the time
        fd.to_excel(fname, 
                    sheet_name='As of ' + self._now.split()[0], )
            

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
        from IPython.display import display as disp
        from IPython.display import Javascript
        from google.auth import default
        from google.colab import auth
        from google.colab import files
        import gspread
        # These will have to be made global
        global disp
        global Javascript
        global default
        global auth
        global files
        global gspread

        # Authenticate the Colab environment 
        auth.authenticate_user()
        creds, _ = default()
        self._gc = gspread.authorize(creds)
        
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
        # Check if "Amount" column is present, and add it if so
        if "Amount" in rxns.columns:
            num_cols.append("Amount")
        for nc in num_cols:
            rxns[nc] = pd.to_numeric(rxns[nc])
        
        self.rxns = rxns
        
    def _get_gsheet_vals(self, key, sheet):
        '''General code for getting Google Sheet values and returning a 
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
        cpd_mask = val_df['Compound'].isna()
        val_df = val_df[~cpd_mask]

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

        fill : str or None, optional (Default = '-')
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


    def excel_save(self, fname, decimals=None):
        '''Download the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str
            The name you want to give to the Excel file.

        decimals : str or None, optional (default = None)
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
        
class WebAppCost(ExcelCost):
    '''Costing class designed for the Streamlit Web App environment.

    This is largely the same as the `ExcelCost` class, so see that method for
    a more complete docstring. 
    '''
    def _excel_csv_reader(self, fname, fsheet, dtypes=None):
        '''A simple Excel reader function for both reaction and materials
        files.

        For the WebApp, it is assumed that only Excel files will be passed
        into the app. In the app, the Excel file is converted into a binary
        stream, so you can not easily tell the difference between Excel/CSV
        files. The read_excel pandas method can still read the bytestream just
        fine.
        '''
        # Read the file, drop NaN-only and commented rows.
        df = pd.read_excel(fname, fsheet, dtype=dtypes, comment='#')\
                        .dropna(how='all')

        # Drop lines that are still empty, these cause all sorts of problems
        # Assume rxn/materials sheets should have Cpd names for valid entries
        cpd_mask = df['Compound'].isna()
        df = df[~cpd_mask]
        
        return df
        
    def results(self, style='compact', decimals=2, fill=np.nan):
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

        fill : str or None, optional (Default = np.nan)
            Fill NaN values with this string. This makes the table a little
            easier to read. Set this to `None` if you want to see the table
            with the typical NaN labels.
        '''
        fd = self._prep_results(style)

        # Display the DataFrames for different permutations of kwargs
        if decimals:
            fd = fd.round(decimals).fillna(fill)
        else:
            fd = fd.fillna(fill)

        return fd

    def excel_save(self, fname, decimals=None):
        '''Save the costing DataFrame as an Excel file.

        Parameters
        ----------
        fname : str
            The name you want to give to the Excel file.

        decimals : str or None, optional (default = None)
            The number of decimal places to display in the Excel sheet. If
            `None`, then the full precision will be saved. 

        Note
        ----
        In some cases, this function will throw an error. In that case, try
        running this again in order to get it to work. 
        '''
        # Get the process DataFrame
        fd = self._prep_excel(decimals)

        # Can set some keyword arguments here
        kwargs = {}
        # If decimals is given, set that value to the rounding for float
        # formatting in the output
        if decimals:
            kwargs['float_format'] = '%.{:d}f'.format(decimals)
           
        # Create the excel file. Can only save with the date and not the time
        # This must be done as a Bytes object, as described in the refs
        output = BytesIO()
        writer = pd.ExcelWriter(output)
        # Dump the data to the ExcelWriter, which in turn sends it to a
        # bytestream object
        fd.to_excel(writer, 
                    sheet_name='As of ' + self._now.split()[0], 
                    **kwargs)
        writer.save()
        
        # Get the byte string for output and return this for saving
        proc_excel = output.getvalue()

        return proc_excel

### References:
# https://discuss.streamlit.io/t/download-button-for-csv-or-xlsx-file/17385/2
