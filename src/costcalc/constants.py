'''Constants that are used across different modules.
'''
# Column name mapper: This is setting variable names for many of the important
# DataFrame columns. 
# Reactions Table Columns
RXN_STP = 'Step' # The step labels
RXN_CPD = 'Compound' # The compound names column, same for materials table
RXN_EQ = 'Equiv' # Molar equivalents of reagent
RXN_MS = 'Mass' # Mass of material used, typically kg, just be consistent
RXN_VOL = 'Volumes' # Volumes of solvent L/kg of limiting reagent (usually)
RXN_RCY = 'Sol Recyc' # Fractional amount of solvent that is recycled
RXN_CST = 'Cost step' # Step where cost is calculated
RXN_OPX = 'OPEX' # Operational expenditures in $/kg, only for rxn product
NOTES = 'Notes' # The notes column name, same for materials table
# Materials Table Columns
MAT_MW = 'MW' # Molecular weight
MAT_DEN = 'Density' # Density, only used for solvents
MAT_CST = '$/kg' # The estimated material price/cost ($/kg)
# Calculated columns
RXN_KG = 'kg/kg rxn' # kg of materials used per kg reaction/step product
RXN_RMC = 'RM cost/kg rxn' # Cost of a material to make 1 kg rxn product
RXN_RMP = '% RM cost/kg rxn' # % material cost per 1 kg rxn product
PRD_KG = 'kg/kg prod' # kg of material used per kg of route product
PRD_RMC = 'RM cost/kg prod' # Cost of a material to make 1 kg route prod
PRD_RMP = '% RM cost/kg prod' # % material cost per 1 kg route product
# Excel temporary column names, these are renamed on Excel export
# These are the same as the calculated columns except for the SUFFIX
SUFFIX = ' dyn'
RNUM = 'rnum' # Row number column for indexing Excel cells
DYN_CST = MAT_CST + SUFFIX # Dynamic product cost
DYN_RKG = RXN_KG + SUFFIX # dynamic kg/kg rxn
DYN_RRMC = RXN_RMC + SUFFIX # dynamic RM cost/kg rxn
DYN_RRMP = RXN_RMP + SUFFIX # dynamic %RM cost/kg rxn
DYN_PKG = PRD_KG + SUFFIX # dynamic kg/kg prod
DYN_PRMC = PRD_RMC + SUFFIX # dynamic RM cost/kg prod
DYN_PRMP = PRD_RMP + SUFFIX # dynamic %RM cost/kg prod


# Excel Column dictionary
# This dictionary will be used for creating dynamic excel sheets, for which
# I'll need to know the DataFrame<->Excel column name mapping
ECOLS = {RXN_STP:'A', RXN_CPD:'B', MAT_MW:'C', MAT_DEN:'D', MAT_CST:'E',
        RXN_EQ:'F', RXN_VOL:'G', RXN_RCY:'H', RXN_CST:'I',
        RXN_OPX:'J', RXN_KG:'K', RXN_RMC:'L', RXN_RMP:'M', PRD_KG:'N',
        PRD_RMC:'O', PRD_RMP:'P', }

### References ###
# https://stackoverflow.com/questions/1383239/can-i-use-init-py-to-define-global-variables
