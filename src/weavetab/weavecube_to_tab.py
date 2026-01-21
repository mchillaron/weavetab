#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy import units as u
from astropy.io.fits import Header
from astropy.table import Table
from pathlib import Path
from tqdm import tqdm
import numpy as np

GREEN   = "\033[92m"
RESET   = "\033[0m"

def weavecube_to_tab(cube_dict, cube_filename, working_dir, number_simulations, region=None,
                            save_tab=False):
    
    """Extract (x,y) coordinates, calibrated spectra and their inverse variance from a WEAVE cube.
    The information is stored in an Astropy Table, and saved as a FITS file in a subdirectory
    within the working directory. WCS information is also included in the FITS header.

    Parameters
    ----------
    cube_dict : dict
        Dictionary containing the WEAVE cube data and metadata.
    cube_filename : str
        Filename of the WEAVE cube (without path).
    working_dir : str
        Path to the working directory where the output FITS file will be saved.
    region : tuple of int, optional
        Region of interest in the cube specified as (x1, x2, y1, y2). FITS indices.
        If None, the full cube is processed.
    number_simulations : int, optional
        Number of simulations to perform. Default is 0.
    """

    data = cube_dict["data"]
    ivar = cube_dict["ivar"]
    sensfunc = cube_dict["sensfunc"]  # allows flux calibration
    header = cube_dict["data_header"]

    naxis1 = header["NAXIS1"]  # spatial axis 1, FITS x axis
    naxis2 = header["NAXIS2"]  # spatial axis 2, FITS y axis
    naxis3 = header["NAXIS3"]  # spectral
    print(f'{GREEN}INFO:{RESET} Cube dimensions - NAXIS1: {naxis1}, NAXIS2: {naxis2}, NAXIS3: {naxis3}')

    if region is not None:
        x1_fits, x2_fits, y1_fits, y2_fits = region
    else:
        x1_fits, x2_fits, y1_fits, y2_fits = 1, naxis1, 1, naxis2

    # Convert to NumPy indices
    x_start = x1_fits - 1
    x_stop  = x2_fits          # inclusive FITS → exclusive NumPy
    y_start = y1_fits - 1
    y_stop  = y2_fits

    # Safety clipping
    x_start = max(0, x_start)
    y_start = max(0, y_start)
    x_stop = min(naxis1, x_stop)
    y_stop = min(naxis2, y_stop)

    rows = []
    num_zero_spectra = 0
    total_pixels = (x_stop - x_start) * (y_stop - y_start)

    for y in tqdm(range(y_start, y_stop), desc="Extracting spectra"):
        for x in range(x_start, x_stop):
            spectrum_adus = data[:, y, x]
            spectrum_ivar = ivar[:, y, x]

            spectrum_flux = spectrum_adus * sensfunc     # Flux calibration
            if np.all(spectrum_flux == 0):               # Skip empty spectra
                num_zero_spectra += 1
                continue

            # sigma(ADU)= 1/sqrt(spectrum_ivar), avoiding division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                sigma_adus = np.where(spectrum_ivar > 0,
                                    1.0 / np.sqrt(spectrum_ivar),
                                    0)

            #sigma_adus = np.where(spectrum_ivar > 0, 1.0 / np.sqrt(spectrum_ivar), 0)  #CHECK
            sigma_flux = sigma_adus * sensfunc           # Flux calibration of sigma
            
            # Store FITS-style coordinates (1-based)
            rows.append((x+1, y+1, spectrum_adus, spectrum_flux, spectrum_ivar, sigma_adus, sigma_flux))

    print(f"{GREEN}INFO:{RESET} Found and discarded {num_zero_spectra} zero spectra")
    print(f"{GREEN}INFO:{RESET} Kept {len(rows)} spectra out of {total_pixels}")

    
    table = Table(rows=rows, names=("x", "y", "specADU", "spec", "ivarADU", "sigmaADU", "sigma"))
    table.meta["Header"] = header
    table["x"].unit = u.pixel
    table["y"].unit = u.pixel
    table["specADU"].unit = u.adu
    table["spec"].unit = u.erg / (u.s * u.cm**2 * u.AA)
    table["ivarADU"].unit = (u.erg / (u.s * u.cm**2 * u.AA))**-2
    table["sigmaADU"].unit = u.adu
    table["sigma"].unit = u.erg / (u.s * u.cm**2 * u.AA)

    # Add WCS information to the FITS header
    #wcs_header = Header()
    #for key in header.keys():
    #    if key.startswith(('CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT', 'CD')):
    #        wcs_header[key] = header[key]
    #        table.meta.update(wcs_header)

    cols_to_save = ["x", "y", "spec", "sigma"]
    table_to_save = table[cols_to_save]
    print('This is a preview of the first 10 rows of the table:')
    table_to_save[:10].pprint(max_width=120)

    # Create output directory if it doesn't exist
    if save_tab:
        output_dir = Path(working_dir) / f"{cube_filename}_spectra_tabs"
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f'{GREEN}INFO:{RESET} Created output directory: {output_dir}')
        else:
            print(f'{GREEN}INFO:{RESET} Output directory already exists: {output_dir}')
            print(f'{GREEN}INFO:{RESET} Existing files will be overwritten.')
    
        # The output filename depends on the case:
        # 1) if region is None, output filename is "extracted_spectra.fits"
        # 2) if region is specified, output filename is "extracted_spectra_x1x2_y1y2.fits"
        # 3) all the above + if number_simulations > 0, append "_000{number_simulations}" 
        # to the filename before the extension, be the number zero-padded to 3 digits.

        if region is None:
            base_name = "extracted_spectra"
        else:
            base_name = (f"extracted_spectra_{x1_fits}-{x2_fits}_{y1_fits}-{y2_fits}")
        if number_simulations > 0:
            base_name += f"_{number_simulations:04d}"
        output_path = output_dir / f"{base_name}.fits"

        table_to_save.write(output_path, format="fits", overwrite=True)
        print(f'{GREEN}INFO:{RESET} Extracted spectra saved to {output_path}')

    return table_to_save

    