import io
from os.path import dirname
from pathlib import Path

import pandas as pd

from costcalc.constants import *

# The data file directory
ddir = Path(dirname(__file__), 'data')

# Read all of the csv files here as strings and format out column names
# Resave the formated strings as StringIO objects to be passed to the read_csv
# function. This allows the csv files to use generic column names

csvs = ddir.glob('*.csv')
IOs = {}
for csv in csvs:
    print(csv)
    with open(csv) as txt:
        txt_file = txt.read().format(**globals())
        IOs[csv.parts[-1]] = io.StringIO(txt_file)

# Load the test data for comparisons
# For some sheets, the dtype must be specified for string columns
# Materials file
mat_clean = pd.read_csv(IOs['clean_mat.csv'])

# Linear and convergent rxn sheets
lin_clean = pd.read_csv(IOs['clean_lin_route.csv'], 
        dtype={RXN_STP:str, RXN_CST:str})
con_clean = pd.read_csv(IOs['clean_conv_route.csv'], 
        dtype={RXN_STP:str, RXN_CST:str})

# Linear and convergent fulldata DataFrames
lin_clean_fd = pd.read_csv(IOs['clean_lin_fd.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_fd = pd.read_csv(IOs['clean_conv_fd.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear and convergent compact results DataFrames
lin_clean_res_c = pd.read_csv(IOs['clean_lin_res_com.csv'],
        dtype={RXN_STP:str, })\
                .set_index([RXN_STP, RXN_CPD])
con_clean_res_c = pd.read_csv(IOs['clean_con_res_com.csv'],
        dtype={RXN_STP:str, })\
                .set_index([RXN_STP, RXN_CPD])

# Linear and convergent full results DataFrames
lin_clean_res_f = pd.read_csv(IOs['clean_lin_res_ful.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_res_f = pd.read_csv(IOs['clean_con_res_ful.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear rxn sheets with duplicated, split materials
lin_split = pd.read_csv(IOs['split_solv_lin_route.csv'],
        dtype={RXN_STP: str, RXN_CST: str})
lin_split_fd = pd.read_csv(IOs['split_solv_lin_fd.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_split_res_f = pd.read_csv(IOs['split_solv_lin_res_ful.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear rxn sheets with OPEX added to each step
lin_opex = pd.read_csv(IOs['opex_lin_route.csv'],
        dtype={RXN_STP: str, RXN_CST: str})
lin_opex_fd = pd.read_csv(IOs['opex_lin_fd.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_opex_res_f = pd.read_csv(IOs['opex_lin_res_ful.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])

# Linear rxn sheets with Mass column
lin_mass = pd.read_csv(IOs['mass_lin_route.csv'],
        dtype={RXN_STP: str, RXN_CST: str})
lin_mass_fd = pd.read_csv(IOs['mass_lin_fd.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_mass_res_f = pd.read_csv(IOs['mass_lin_res_ful.csv'],
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

lin_clean_excel = pd.read_csv(IOs['clean_lin_excel.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_clean_excel[MAT_CST] = lin_clean_excel[MAT_CST]\
        .apply(dtype_fix)

con_clean_excel = pd.read_csv(IOs['clean_con_excel.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
con_clean_excel[MAT_CST] = con_clean_excel[MAT_CST]\
        .apply(dtype_fix)

lin_split_excel = pd.read_csv(IOs['split_solv_lin_excel.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_split_excel[MAT_CST] = lin_split_excel[MAT_CST]\
        .apply(dtype_fix)

lin_opex_excel = pd.read_csv(IOs['opex_lin_excel.csv'],
        dtype={RXN_STP:str, RXN_CST:str})\
                .set_index([RXN_STP, RXN_CPD])
lin_opex_excel[MAT_CST] = lin_opex_excel[MAT_CST]\
        .apply(dtype_fix)



