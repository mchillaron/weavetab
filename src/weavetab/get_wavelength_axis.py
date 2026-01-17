#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def get_wavelength_axis(cube):
    """Generates the wavelength axis from the FITS header information.
    Parameters
    ----------
    cube : dict
        Dictionary containing the data and headers for the WEAVE cube.
    Returns
    -------
    wavelength : ndarray
        1D array containing the wavelength values corresponding to each spectral slice.
    """

    header = cube["data_header"]

    crval3 = header["CRVAL3"]  # reference value
    crpix3 = header["CRPIX3"]  # reference pixel
    naxis3 = header["NAXIS3"]  # number of elements along the wavelength axis
    cdelt3 = header.get("CD3_3", header.get("CDELT3", None))
    if cdelt3 is None:
        raise KeyError("CD3_3 nor CDELT3 not found in header.")

    # Standard definition of WCS in FITS
    indices = np.arange(1, naxis3 + 1)  # FITS counts from 1
    wavelength = crval3 + (indices - crpix3) * cdelt3

    return wavelength