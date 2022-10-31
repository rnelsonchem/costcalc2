from os.path import dirname
from pathlib import Path

import pandas as pd

from costcalc.constants import *

# The data file directory
ddir = Path(dirname(__file__), 'data')

# Load the test data for comparisons
# For some sheets, the dtype must be specified for string columns
# Materials file
mat_clean = pd.read_csv(ddir / 'clean_mat.csv')

# Linear and convergent rxn sheets
lin_clean = pd.read_csv(ddir / 'clean_lin_route.csv', 
        dtype={RXN_STP:str, RXN_CST:str})
con_clean = pd.read_csv(ddir / 'clean_conv_route.csv', 
        dtype={RXN_STP:str, RXN_CST:str})

# Linear and convergent fulldata DataFrames
lin_clean_fd = pd.read_csv(ddir / 'clean_lin_fd.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_fd = pd.read_csv(ddir / 'clean_conv_fd.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear and convergent compact results DataFrames
lin_clean_res_c = pd.read_csv(ddir / 'clean_lin_res_com.csv',
        dtype={RXN_STP:str, })\
                .set_index([RXN_STP, RXN_CPD])
con_clean_res_c = pd.read_csv(ddir / 'clean_con_res_com.csv',
        dtype={RXN_STP:str, })\
                .set_index([RXN_STP, RXN_CPD])

# Linear and convergent full results DataFrames
lin_clean_res_f = pd.read_csv(ddir / 'clean_lin_res_ful.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_res_f = pd.read_csv(ddir / 'clean_con_res_ful.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear and convergent Excel DataFrames
# The "Cost" column has multiple dtypes, that are not set correctly when the
# data is imported via Pandas. This function fixes the dyptes
def dtype_fix(x):
    try:
        return float(x)
    except:
        return x

lin_clean_excel = pd.read_csv(ddir / 'clean_lin_excel.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_clean_excel['Cost'] = lin_clean_excel['Cost']\
        .apply(dtype_fix)
con_clean_excel = pd.read_csv(ddir / 'clean_con_excel.csv',
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_excel['Cost'] = con_clean_excel['Cost']\
        .apply(dtype_fix)



