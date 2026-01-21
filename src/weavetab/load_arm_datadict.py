#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits
import numpy as np

RED     = "\033[91m"
GREEN   = "\033[92m"
RESET   = "\033[0m"

def load_arm_datadict(filepath):
    """Load WEAVE cube extensions for a given arm into a data dictionary.
    The function automatically detects whether the cube corresponds to the RED or BLUE arm.
    Parameters
    ----------
    filepath : Path
        Path to the WEAVE FITS data cube file.
    Returns
    -------
    cube : dict
        Dictionary containing the data and headers for the detected arm.
    arm : str
        Detected arm ("RED" or "BLUE").
    """

    cube = {
        "arm": None,
        "primary_header": None,
        "data_header": None,
        "data": None,
        "ivar": None,
        "data_noss": None,
        "ivar_noss": None,
        "sensfunc": None,
    }

    with fits.open(filepath, memmap=False) as hdul:

        # First, detect which arm we are working with
        arms_found = set()
        for hdu in hdul:
            name = hdu.name.upper()
            if name.endswith("_DATA"):
                arm = name.replace("_DATA", "")
                if arm in {"RED", "BLUE"}:
                    arms_found.add(arm)

        if len(arms_found) == 0:
            raise ValueError("No RED or BLUE arm found in FITS file.")
        if len(arms_found) > 1:
            raise ValueError(f"Multiple arms found: {arms_found}")

        arm = arms_found.pop()
        cube["arm"] = arm
        print(f"{GREEN}INFO:{RESET} Detected arm: {arm}")

        # Secondly, load the extensions into a data dictionary
        for hdu in hdul:
            name = hdu.name.upper()

            if name == "PRIMARY":
                cube["primary_header"] = hdu.header

            elif name == f"{arm}_DATA":
                cube["data_header"] = hdu.header
                cube["data"] = hdu.data.astype(np.float32)

            elif name == f"{arm}_IVAR":
                cube["ivar"] = hdu.data.astype(np.float32)

            elif name == f"{arm}_DATA_NOSS":
                cube["data_noss"] = hdu.data.astype(np.float32)

            elif name == f"{arm}_IVAR_NOSS":
                cube["ivar_noss"] = hdu.data.astype(np.float32)

            elif name == f"{arm}_SENSFUNC":
                cube["sensfunc"] = hdu.data.astype(np.float32)

    return cube, arm
