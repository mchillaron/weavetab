#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from pathlib import Path

import argparse
import os
import re
import numpy as np

from weavetab.load_arm_datadict import load_arm_datadict
from weavetab.get_wavelength_axis import get_wavelength_axis
from weavetab.weavecube_to_tab import weavecube_to_tab
from weavetab.table_to_cube import table_to_cube
from weavetab.simulate_cubes import simulate_cubes

RED     = "\033[91m"
GREEN   = "\033[92m"
YELLOW  = "\033[93m"
BLUE    = "\033[94m"
MAGENTA = "\033[95m"
CYAN    = "\033[96m"
BOLD = "\033[1m"
RESET   = "\033[0m"

def extract_cube_spectra(working_dir, cube_path, number_simulations, region=None, save_tab=False):
    """Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.
    Parameters
    ----------
    working_dir : Path
        Current working directory.
    cube_path : Path
        Path to the WEAVE FITS data cube file.
    number_simulations : int
        Number of simulations to perform.
    region : tuple of int, optional
        Region of interest in the cube specified as (x1, x2, y1, y2). 
        If None, the full cube is processed.
    save_tab : 
    -----------
    """

    print(f"{BOLD}{MAGENTA} Loading cube and detecting WEAVE arm{RESET}")
    cube_dict, arm = load_arm_datadict(cube_path)
    print(f'{GREEN}INFO:{RESET} Loaded data cube successfully.')

    cube_filename = cube_path.stem  # Get the cube filename without extension
    output_dir = Path(working_dir) / f"{cube_filename}_spectra_cubes"
    # Create output directory if it doesn't exist
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f'{GREEN}INFO:{RESET} Created output directory: {output_dir}')
    else:
        print(f'{GREEN}INFO:{RESET} Output directory already exists: {output_dir}')
        # Alert the user that existing files may be overwritten
        existing_files = list(output_dir.glob('*'))
        if existing_files:
            print(f'{YELLOW}WARNING:{RESET} The following files already exist in the output directory and may be overwritten:')
            for file in existing_files:
                print(f'  - {file.name}')
            # ask the user if they want to proceed; press enter to proceed, or type 'n' to abort
            proceed = input(f"{YELLOW}WARNING:{RESET} Do you want to proceed and overwrite these files? (y/n): ").strip().lower()
            if proceed != 'y':
                print(f"{RED}ABORTED:{RESET} Operation cancelled by the user. No files were overwritten.")
                return
        else:
            print(f'{GREEN}INFO:{RESET} No existing files found in the output directory. Safe to proceed.')

    print(f"{BOLD}{MAGENTA} Obtaining the wavelength range and saving to CSV{RESET}")
    # Check if a wavelength_{arm} CSV file exists in the working directory, if not, create one. If it exists, load it as a 1D array.
    wavelength_csv = output_dir / f'wavelength_{arm}.csv'         
    if not wavelength_csv.is_file():
        wavelength = get_wavelength_axis(cube_dict)
        np.savetxt(wavelength_csv, wavelength, delimiter=",")
        print(f'{GREEN}INFO:{RESET} Created wavelength axis CSV file: {wavelength_csv}')
    else:
        print(f'{GREEN}INFO:{RESET} Wavelength axis CSV file already exists: {wavelength_csv}')
        wavelength = np.loadtxt(wavelength_csv, delimiter=',') # and load it as 1d array

    print(f'{BOLD}{MAGENTA} Extracting spectra from the original cube (no simulations) to Table.{RESET}')
    table, table_header = weavecube_to_tab(cube_dict, output_dir, number_simulations, region, save_tab)

    print(f'{BOLD}{MAGENTA} Transforming the FITS Table to Cube.{RESET}')
    new_cube, new_cube_path = table_to_cube(table, table_header, cube_filename, output_dir, region)
    
    if number_simulations > 0:
        print(f'{BOLD}{MAGENTA} Simulating WEAVE cubes.{RESET}')
        simulate_cubes(new_cube, output_dir, region, number_simulations)

    print(f"Goodbye!")
    
    


def main():
    parser = argparse.ArgumentParser(description="Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.")
    parser.add_argument('-F', '--cube-file', type=str, help="Name of the WEAVE FITS data cube file.")
    parser.add_argument('-reg', '--region', type=int, nargs=4, metavar=("x1", "x2", "y1", "y2"), help="Region of interest in the cube: x1 x2 y1 y2. FITS indices.")
    parser.add_argument('-n', '--number-simulations', type=int, default=0, help="Number of simulations to perform.")
    parser.add_argument('-s', '--save-table', type=bool, default=False, help="Save a FITS table containing spatial coordinates and spectrum per pixel. Default False.")
    args = parser.parse_args()

    print(f"{BOLD} Welcome to weavetab {RESET}")

    cube_file = args.cube_file
    number_simulations = args.number_simulations
    save_tab = args.save_table

    if args.region is None:
        # process full cube
        region = None
    else:
        region = args.region
        print(f"{GREEN}INFO:{RESET} Region of interest specified: {region}")
        x1, x2, y1, y2 = region
        if x2 is not None and x2 <= x1:
            raise ValueError("x2 must be greater than x1.")
        if y2 is not None and y2 <= y1:
            raise ValueError("y2 must be greater than y1.")

    # cube filename must have a .fit, .fits, .FIT, or .FITS extension
    if not re.match(r'.*\.(fit|fits|FIT|FITS)$', cube_file):
        raise ValueError("Cube file must have a .fit, .fits, .FIT, or .FITS extension.")
    if not os.path.isfile(cube_file):
        raise FileNotFoundError(f"Cube file '{cube_file}' does not exist.")
    
    working_dir = Path.cwd()
    print(f"{GREEN}INFO:{RESET} Current working directory: {working_dir}")
    cube_path = working_dir / cube_file
    cube_path = cube_path.resolve()                  # get the absolute path of the cube file
    print(f"{GREEN}INFO:{RESET} processing cube file: {cube_path}")
    
    if number_simulations < 0:
        raise ValueError("Number of simulations must be a non-negative integer.")
    if number_simulations > 0:
        print(f"{GREEN}INFO:{RESET} Number of simulations to perform: {number_simulations}")
        print(f"{GREEN}INFO:{RESET} A total number of {number_simulations+1} FITS cubes will be created.")
    if number_simulations == 0:
        print(f"{GREEN}INFO:{RESET} No simulations will be performed. Only the original spectra will be extracted.")

    if save_tab:
        print(f"{GREEN}INFO:{RESET} A FITS table will be saved containing spatial coordinates and spectrum per pixel.")
    else:
        print(f"{GREEN}INFO:{RESET} No FITS table will be saved. Only a FITS cube.")

    extract_cube_spectra(working_dir, cube_path, number_simulations, region, save_tab)


if __name__ == "__main__":
    main()