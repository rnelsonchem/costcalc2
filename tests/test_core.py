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

# Use Pandas to read in the test Excel files 
# For the reactions sheets, the dtype must be specified for a couple of
# columns
mats = pd.read_excel(ddir / 'clean_mat.xlsx')
br_clean = pd.read_excel(ddir / 'clean_rxn.xlsx', 
        sheet_name='Bromine Route',
        dtype={RXN_STP:str, RXN_CST:str})
cl_clean = pd.read_excel(ddir / 'clean_rxn.xlsx', 
        sheet_name='Chlorine Route',
        dtype={RXN_STP:str, RXN_CST:str})
br_clean_fd = pd.read_pickle(ddir / 'br_clean_fd.pickle')
cl_clean_fd = pd.read_pickle(ddir / 'cl_clean_fd.pickle')


### Tests ###

class Test_CoreFunctions(object):
    @pytest.mark.parametrize(
            "rxn, mats, fd",
            [(br_clean, mats, br_clean_fd),
            (cl_clean, mats, cl_clean_fd),
            ]
    )
    def test_rxn_data_setup(self, rxn, mats, fd):
        coster = CoreCost(mats, rxn, 'Product')
        assert_frame_equal(coster.fulldata, fd)

    @pytest.mark.parametrize(
            "rxn, mats, cost",
            [(br_clean, mats, 43.05580071473825,),
            (cl_clean, mats, 61.25477348169882,),
            ]
    )
    def test_calc_cost_clean(self, rxn, mats, cost):
        coster = CoreCost(mats, rxn, 'Product')
        coster.calc_cost()
        assert_allclose(coster.cost, cost)


