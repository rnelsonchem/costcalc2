from os.path import dirname
from pathlib import Path

import pytest
import pandas as pd

# Using some of the numerical testing routines in Numpy
from numpy.testing import *

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


### Tests ###

@pytest.mark.parametrize(
        "rxn, mats, cost",
        [(br_clean, mats, 43.05580071473825),
        (cl_clean, mats, 61.25477348169882),
        ]
)
class Test_CoreFunctions(object):
    def test_calc_cost_clean(self, rxn, mats, cost):
        coster = CoreCost(mats, df, 'Product')
        coster.calc_cost()
        assert_allclose(coster.cost, cost)


