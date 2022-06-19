import numpy as np

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

SORTING_OPTIONS = [
    "Higher is better",
    "Lower is better",
]