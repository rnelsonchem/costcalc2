import pytest
import pandas as pd
import numpy as np

# Using some of the testing routines from third-party libraries
from numpy.testing import *
from pandas.testing import *

from costcalc import CoreCost
from costcalc.constants import *
from costcalc.exceptions import CostError
from .data_build import *

### Tests - Working reactions ###

class Test_CoreFunctions(object):
    '''Tests for the public functions in the Core class.

    The test names here incorporate the public function name. These checks are
    parameterized with generic linear (lin) and convergent (con) routes, which
    is the minimal types of routes.  
    '''
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


class Test_CoreErrors(object):
    '''Test for known errors in the Core class.

    These errors are all captured in the code, so these are forcing the known
    errors to make sure they are correctly found/diagnosed.
    '''
    def test_miss_mw(self, ):
        '''Test for a missing material.
        '''
        mat_clean_miss = mat_clean.copy().drop(0)
        txt_match = 'Missing MW error!'
        with pytest.raises(CostError, match=txt_match) as err:
            coster = CoreCost(mat_clean_miss, lin_clean, 'Product')
        assert err.value.df.index == ('1', 'Starting Material')

    def test_missing_cost(self, ):
        '''Test for a missing cost or cost calculation.
        '''
        mat_clean_miss = mat_clean.copy()
        mask = mat_clean_miss[RXN_CPD] == 'Reagent C'
        mat_clean_miss.loc[mask, MAT_CST] = np.nan 
        txt_match = 'Missing material cost'
        with pytest.raises(CostError, match=txt_match) as err:
            coster = CoreCost(mat_clean_miss, lin_clean, 'Product')
        assert err.value.df.index == ('3', 'Reagent C')

    def test_duplicate_mat(self, ):
        '''Test for duplicated compound in the materials table.
        '''
        mat_clean_miss = mat_clean.copy()
        first_cpd = mat_clean_miss.iloc[:1,:]
        mat_clean_miss = pd.concat([mat_clean_miss, first_cpd])
        txt_match = 'Duplicated material error'
        with pytest.raises(CostError, match=txt_match) as err:
            coster = CoreCost(mat_clean_miss, lin_clean, 'Product')
        assert err.value.df.iloc[0, 0] == 'Starting Material'

    def test_duplicate_rxn(self, ):
        '''Test for duplicated compounds in a single rxn step.
        '''
        lin_clean_dup = lin_clean.copy()
        bromine = lin_clean_dup.iloc[1:2,:]
        lin_clean_dup = pd.concat([lin_clean_dup.iloc[:1,:], bromine, 
                            lin_clean_dup.iloc[1:,:]]).reset_index(drop=True)
        txt_match = 'Duplicated material in a reaction step'
        with pytest.raises(CostError, match=txt_match) as err:
            # This also throws a Pandas warning because of the index is no
            # longer sorted. This captures the warning as well.
            with pytest.warns(pd.errors.PerformanceWarning):
                coster = CoreCost(mat_clean, lin_clean_dup, 'Product')
        assert err.value.df.index[0] == ('1', 'Bromine')

    def test_missing_sol_density(self, ):
        '''Test for a missing density for solvent.
        '''
        mat_clean_dup = mat_clean.copy()
        mask = (mat_clean[RXN_CPD] == 'Water')
        mat_clean_dup.loc[mask, MAT_DEN] = np.nan
        txt_match = 'Missing solvent info error!'
        with pytest.raises(CostError, match=txt_match) as err:
            coster = CoreCost(mat_clean_dup, lin_clean, 'Product')
        assert err.value.df.index == ('2', 'Water')


