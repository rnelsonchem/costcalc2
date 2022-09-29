from os.path import dirname
from pathlib import Path

import pytest
import pandas as pd

# Using some of the numerical testing routines in Numpy
from numpy.testing import *

from costcalc import CoreCost
from costcalc.core import rxn_stp, rxn_cst

# The data file directory
ddir = Path(dirname(__file__), 'data')

class Test_CoreFunctions(object):

    @pytest.mark.parametrize(
            "name, cost",
            [('Bromine Route', 43.05580071473825),
            ('Chlorine Route', 61.25477348169882),
            ]
    )
    def test_calc_cost_clean(self, name, cost):
        # Use Pandas to read in the test Excel files
        # For the reactions sheets, the dtype must be specified for a couple of
        # columns
        mats = pd.read_excel(ddir / 'clean_mat.xlsx')
        rxns = pd.read_excel(ddir / 'clean_rxn.xlsx',
                sheet_name = name,
                dtype={rxn_stp:str, rxn_cst:str})

        coster = CoreCost(mats, rxns, 'Product')
        coster.calc_cost()
        assert_allclose(coster.cost, cost)


