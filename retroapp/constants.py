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
