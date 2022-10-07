from os.path import dirname
from pathlib import Path

import pytest
import pandas as pd

# Using some of the testing routines from third-party libraries
from numpy.testing import *
from pandas.testing import *

from costcalc import CoreCost
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


### Tests - Working reactions ###

class Test_CoreFunctions(object):
    @pytest.mark.parametrize(
            "rxn, mat, fd",
            [(lin_clean, mat_clean, lin_clean_fd),
            (con_clean, mat_clean, con_clean_fd),
            ]
    )
    def test_rxn_data_setup(self, rxn, mat, fd):
        coster = CoreCost(mat, rxn, 'Product')
        assert_frame_equal(coster.fulldata, fd)

    @pytest.mark.parametrize(
            "rxn, mat, cost",
            [(lin_clean, mat_clean, 66.40818831355247,),
            (con_clean, mat_clean, 62.9051883439257,),
            ]
    )
    def test_calc_cost_clean(self, rxn, mat, cost):
        coster = CoreCost(mat, rxn, 'Product')
        coster.calc_cost()
        assert_allclose(coster.cost, cost)

    @pytest.mark.parametrize(
            "rxn, mat, res, style",
            [(lin_clean, mat_clean, lin_clean_res_c, 'compact'),
            (con_clean, mat_clean, con_clean_res_c, 'compact'),
            (lin_clean, mat_clean, lin_clean_res_f, 'full'),
            (con_clean, mat_clean, con_clean_res_f, 'full'),
            ]
    )
    def test_results(self, rxn, mat, res, style):
        coster = CoreCost(mat, rxn, 'Product')
        coster.calc_cost()
        assert_frame_equal(coster.results(style=style), res)


    @pytest.mark.parametrize(
            "rxn, mat, excel",
            [(lin_clean, mat_clean, lin_clean_excel,),
            (con_clean, mat_clean, con_clean_excel,),
            ]
    )
    def test_excel(self, rxn, mat, excel):
        coster = CoreCost(mat, rxn, 'Product')
        coster.calc_cost()
        assert_frame_equal(coster.excel(), excel)


