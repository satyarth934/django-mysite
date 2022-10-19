import numpy as np

def create_twoway_dict(input_dict):
    temp_dict = dict()
    for k, v in input_dict.items():
        temp_dict[k] = v
        temp_dict[v] = k

    return temp_dict


MOLECULE_PROPERTIES = {
    "Cetane Number": {
        "min": -np.inf,
        "max": np.inf,
    },
    "Research Octane Number": {
        "min": -np.inf,
        "max": np.inf,
    },
    "Melting Point": {
        "min": -np.inf,
        "max": np.inf,
    },
    "Flash Point": {
        "min": -np.inf,
        "max": np.inf,
    },
    "Yield Sooting Index": {
        "min": -np.inf,
        "max": np.inf,
    },
    "H1 Receptor pKd": {
        "min": -np.inf,
        "max": np.inf,
    },
    "M2 Receptor pKd": {
        "min": -np.inf,
        "max": np.inf,
    },
}

PROPERTY_CODES = {
    "CN": "Cetane Number",
    "FP": "Flash Point",
    "H1": "H1 Receptor pKd",
    "M2": "M2 Receptor pKd",
    "MP": "Melting Point",
    "RON": "Research Octane Number",
    "YSI": "Yield Sooting Index",
}
PROPERTY_CODES = create_twoway_dict(PROPERTY_CODES)

DB_COLUMN_NAMES = {
    "Cetane Number": "Cetane_number",
    "Research Octane Number": "Research_octane_number",
    "Melting Point": "Melting_point",
    "Flash Point": "Flash_point",
    "Yield Sooting Index": "Yield_sooting_index",
    "H1 Receptor pKd": "H1_receptor_pKd",
    "M2 Receptor pKd": "M2_receptor_pKd",

    "Cetane_Number": "Cetane_number",
    "Research_Octane_Number": "Research_octane_number",
    "Melting_Point": "Melting_point",
    "Flash_Point": "Flash_point",
    "Yield_Sooting_Index": "Yield_sooting_index",
    "H1_Receptor_pKd": "H1_receptor_pKd",
    "M2_Receptor_pKd": "M2_receptor_pKd",

    "CN": "Cetane_number",
    "FP": "Research_octane_number",
    "H1": "Melting_point",
    "M2": "Flash_point",
    "MP": "Yield_sooting_index",
    "RON": "H1_receptor_pKd",
    "YSI": "M2_receptor_pKd",
}

# NO_SORT = "No sorting"
# ASCENDING = "Lower is better"
# DESCENDING = "Higher is better"

# SORTING_OPTIONS = [
#     ASCENDING,
#     DESCENDING,
#     NO_SORT,
# ]

ASCENDING = 0
DESCENDING = 1
NO_SORT = 2

SORTING_OPTIONS = [
    (ASCENDING, "Lower is better"),
    (DESCENDING, "Higher is better"),
    (NO_SORT, "No sorting"),
]

# SORTING_OPTIONS_DB = (
#     ("ASC", ASCENDING),
#     ("DESC", DESCENDING),
#     ("NOSO", NO_SORT),
# )
