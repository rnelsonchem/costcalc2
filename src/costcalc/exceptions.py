from .constants import *

class CostError(Exception):
    def __init__(self, message, df=None, disp_df=False):
        self.message = message
        self.df = df
        self.disp_df = disp_df
        super().__init__(message)

    def __str__(self, ):
        if self.disp_df:
            print(self.message)
            print('See the following entries:')
            disp(self.df)
        return self.message


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

    'dup_rxn': 'Duplicated material in a reaction step error!\n'\
            'Please combine multiple uses of a material into one entry.',

    'mis_cst': 'Missing material cost error!\n'\
           'This may be due to a missing material and/or its cost.\n'\
           f'Otherwise, you may be missing a Step in the "{RXN_CST}" column.',

    'mis_sol': 'Missing solvent info error!\n'\
           'One of your solvent entries is missing critical info.',
           
    'mis_rel': f'"{RXN_REL}" solvent entry error!\n'\
           f'A "{RXN_REL}" compound for one of your solvents is incorrect.\n'\
           f'This may be due to a name mismatch with the "{RXN_CPD}" column.',
           
    }
    
