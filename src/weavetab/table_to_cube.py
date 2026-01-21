#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#


from astropy.io import fits
from pathlib import Path
import numpy as np

GREEN   = "\033[92m"
CYAN    = "\033[96m"
RESET   = "\033[0m"

def table_to_cube(table, cube_filename, working_dir, region, number_simulations):
    """
    Rebuild a FITS cube from an Astropy Table produced by extract_spectra_to_tab().
    Creates:
      - Primary HDU with full header
      - Extension 1: data cube (spec)
      - Extension 2: sigma cube (sigma)
    """

    header = table.meta["Header"]

    nw = header["NAXIS3"]

    if region is not None:
        x1_fits, x2_fits, y1_fits, y2_fits = region
        nx = x2_fits - x1_fits + 1 # +1 so we don't lose the last pixel
        ny = y2_fits - y1_fits + 1
    else:
        nx = header["NAXIS1"]
        ny = header["NAXIS2"]
        x1_fits, x2_fits, y1_fits, y2_fits = 1, nx, 1, ny


    # Prepare empty cubes
    cube = np.zeros((nw, ny, nx), dtype=float)
    sigma_cube = np.zeros((nw, ny, nx), dtype=float)

    # Fill cubes
    for row in table:
        # Convert global FITS coords → local region coords → numpy coords
        x_local = row["x"] - x1_fits
        y_local = row["y"] - y1_fits

        cube[:, y_local, x_local] = row["spec"]
        sigma_cube[:, y_local, x_local] = row["sigma"]

        #x = row["x"] - 1   # FITS → numpy
        #y = row["y"] - 1
        #cube[:, y, x] = row["spec"]
        #sigma_cube[:, y, x] = row["sigma"]

    # Build HDUs
    primary_hdu = fits.PrimaryHDU(header=header)
    cube_hdu = fits.ImageHDU(data=cube, name="DATA")
    sigma_hdu = fits.ImageHDU(data=sigma_cube, name="SIGMA")

    hdul = fits.HDUList([primary_hdu, cube_hdu, sigma_hdu])

    output_dir = Path(working_dir) / f"{cube_filename}_spectra_cubes"
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f'{GREEN}INFO:{RESET} Created output directory: {output_dir}')
    else:
        print(f'{GREEN}INFO:{RESET} Output directory already exists: {output_dir}')
        print(f'{GREEN}INFO:{RESET} Existing files will be overwritten.')

    if region is None:
        base_name = "cube_spectra"
    else:
        base_name = (f"cube_spectra_{x1_fits}-{x2_fits}_{y1_fits}-{y2_fits}")
    if number_simulations > 0:
        base_name += f"_{number_simulations:04d}"

    print(f"{CYAN}Dimensions of the new cube:{RESET} (nw={nw}, ny={ny}, nx={nx})")

    output_path = output_dir / f"{base_name}.fits"
    hdul.writeto(output_path, overwrite=True)
    print(f"{GREEN}INFO:{RESET} Cube reconstructed and saved to {output_path}")
