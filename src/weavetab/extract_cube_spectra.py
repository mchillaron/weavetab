#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.table import Table
from pathlib import Path

import argparse
import os
import re
import numpy as np

from weavetab.load_arm_datadict import load_arm_datadict
from weavetab.get_wavelength_axis import get_wavelength_axis

def extract_cube_spectra(working_dir, cube_path, number_simulations):
    """Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.
    Parameters
    ----------
    working_dir : Path
        Current working directory.
    cube_path : Path
        Path to the WEAVE FITS data cube file.
    number_simulations : int
        Number of simulations to perform.
    """

    cube_dict, arm = load_arm_datadict(cube_path)
    print('INFO: Loaded data cube successfully.')

    # Check if a wavelength_{arm} CSV file exists in the working directory, if not, create one.
    wavelength_csv = working_dir / f'wavelength_{arm}.csv'

    if not wavelength_csv.is_file():
        wavelength = get_wavelength_axis(cube_dict)
        np.savetxt(wavelength_csv, wavelength, delimiter=",")
        print(f'INFO: Created wavelength axis CSV file: {wavelength_csv}')
    else:
        print(f'INFO: Wavelength axis CSV file already exists: {wavelength_csv}')
        wavelength = np.loadtxt(wavelength_csv, delimiter=',') # and load it as 1d array

    naxis1 = cube_dict["data_header"]["NAXIS1"]  # spatial axis 1, FITS x axis
    naxis2 = cube_dict["data_header"]["NAXIS2"]  # spatial axis 2, FITS y axis
    naxis3 = cube_dict["data_header"]["NAXIS3"]  # spectral
    print(f'INFO: Cube dimensions - NAXIS1: {naxis1}, NAXIS2: {naxis2}, NAXIS3: {naxis3}')
    
    # Extract spectra and spatial coordinates from the cube
    # make an astropy table with columns: x, y, flux (1d array of length naxis3), ivar (1d array of length naxis3)
    



    


def main():
    parser = argparse.ArgumentParser(description="Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.")
    parser.add_argument('-F', '--cube-file', type=str, help="Name of the WEAVE FITS data cube file.")
    parser.add_argument('-n', '--number-simulations', type=int, default=0, help="Number of simulations to perform.")
    args = parser.parse_args()

    cube_file = args.cube_file
    number_simulations = args.number_simulations

    # cube filename must have a .fit, .fits, .FIT, or .FITS extension
    if not re.match(r'.*\.(fit|fits|FIT|FITS)$', cube_file):
        raise ValueError("Cube file must have a .fit, .fits, .FIT, or .FITS extension.")
    if not os.path.isfile(cube_file):
        raise FileNotFoundError(f"Cube file '{cube_file}' does not exist.")
    
    working_dir = Path.cwd()
    print(f"INFO: Current working directory: {working_dir}")
    cube_path = working_dir / cube_file
    cube_path = cube_path.resolve()                  # get the absolute path of the cube file
    print(f"INFO: processing cube file: {cube_path}")
    
    if number_simulations < 0:
        raise ValueError("Number of simulations must be a non-negative integer.")
    if number_simulations > 0:
        print(f"INFO: Number of simulations to perform: {number_simulations}")
        print(f"INFO: A total number of {number_simulations} FITS tables will be created.")
    if number_simulations == 0:
        print("INFO: No simulations will be performed. Only the original spectra will be extracted.")

    extract_cube_spectra(working_dir, cube_path, number_simulations)


if __name__ == "__main__":
    main()