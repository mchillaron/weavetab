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
import numpy as np
import time

from weavetab.load_arm_datadict import load_arm_datadict
from weavetab.get_wavelength_axis import get_wavelength_axis
from weavetab.weavecube_to_tab import weavecube_to_tab
from weavetab.table_to_cube import table_to_cube
from weavetab.simulate_cubes import simulate_cubes
from weavetab.read_continuum_reg_file import read_continuum_file


def info(msg):
    GREEN   = "\033[92m"
    RESET   = "\033[0m"
    print(f"{GREEN}INFO:{RESET} {msg}")

def infostep(msg):
    BOLD = "\033[1m"
    MAGENTA = "\033[95m"
    RESET   = "\033[0m"
    print(f"{MAGENTA}{BOLD} {msg} {RESET}")

def extract_cube_spectra(working_dir, cube_path, number_simulations, debug_level, start, 
                         region, save_tab, cont_regions, seed):
    """Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.
    Parameters
    ----------
    working_dir : Path
        Current working directory.
    cube_path : Path
        Path to the WEAVE FITS data cube file.
    number_simulations : int
        Number of simulations to perform.
    debug_level : int
        Debug level. 0: no debug, 1: print more detailed intermediate results. Default 0.
    region : tuple of int, optional
        Region of interest in the cube specified as (x1, x2, y1, y2). 
        If None, the full cube is processed.
    save_tab : bool, optional
        Whether to save a FITS table containing spatial coordinates and spectrum per pixel. Default is False.
    cont_regions : array-like, optional
        Continuum regions to compute the scaling factor between flux rms and sigma cube. The array should have shape (N, 2) where each row is (lmin, lmax) in Angstroms. If not provided, no scaling will be applied to the sigma cube.
    seed : int
        Seed for initializing cube simulations.  
    -----------
    """

    infostep("Loading cube and detecting WEAVE arm")
    cube_dict, arm = load_arm_datadict(cube_path)
    info('Loaded data cube successfully.')

    cube_filename = cube_path.stem  # Get the cube filename without extension
    output_dir = Path(working_dir) / f"{cube_filename}_spectra_cubes"

    # Create output directory if it doesn't exist

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        info(f'Created output directory: {output_dir}')
    else:
        info(f'Output directory already exists: {output_dir}')
        # Alert the user that existing files may be overwritten
        existing_files = list(output_dir.glob('*'))
        if existing_files:
            info('The following files already exist in the output directory and may be overwritten:')
            for file in existing_files:
                print(f'  - {file.name}')
            # ask the user if they want to proceed; press enter to proceed, or type 'n' to abort
            proceed = input(f"Do you want to proceed and overwrite these files? (y/n): ").strip().lower()
            if proceed != 'y':
                info("Operation cancelled by the user. No files were overwritten.")
                return
        else:
            info('No existing files found in the output directory. Safe to proceed.')


    # Check if a wavelength_{arm} CSV file exists in the working directory, 
    # if not, create one. If it exists, load it as a 1D array.

    infostep("Obtaining the wavelength range and saving to CSV")
    wavelength_csv = output_dir / f'wavelength_{arm}.csv'         
    if not wavelength_csv.is_file():
        wavelength = get_wavelength_axis(cube_dict)
        np.savetxt(wavelength_csv, wavelength, delimiter=",")
        info(f'Created wavelength axis CSV file: {wavelength_csv}')
    else:
        info(f'Wavelength axis CSV file already exists: {wavelength_csv}')
        wavelength = np.loadtxt(wavelength_csv, delimiter=',') # and load it as 1d array

    infostep('Extracting spectra from the original cube (no simulations) to Table')
    table, table_header = weavecube_to_tab(cube_dict, output_dir, number_simulations, region, save_tab)

    infostep('Transforming the FITS Table to Cube')
    new_cube, new_cube_path = table_to_cube(table, table_header, cube_filename, output_dir, region)
    
    if number_simulations > 0:
        infostep('Simulating WEAVE cubes')
        simulate_cubes(new_cube, wavelength, output_dir, region, number_simulations, cont_regions, seed, debug_level)

    print(f"Goodbye!")
    end = time.perf_counter()
    elapsed = end - start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    print(f"⏱️ Execution time: {elapsed:.3f} s ({hours:02d}:{minutes:02d}:{seconds:06.3f})")
    
    

def main():
    start = time.perf_counter()
    parser = argparse.ArgumentParser(description="Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.")
    parser.add_argument('-F', '--cube-file', required=True, help="Name of the WEAVE FITS data cube file.")
    parser.add_argument('-reg', '--region', type=int, nargs=4, metavar=("x1", "x2", "y1", "y2"), help="Region of interest in the cube: x1 x2 y1 y2. FITS indices.")
    parser.add_argument('-n', '--number-simulations', type=int, default=0, help="Number of simulations to perform.")
    parser.add_argument('-s', '--save-table', action='store_true', help="Save a FITS table containing spatial coordinates and spectrum per pixel. Default False.")
    parser.add_argument('-cont', '--cont-regions', type=str, default=None, help="Path to a .txt file containing the continuum regions to compute the scaling factor between flux rms and sigma fron inverse invariance. " \
                                                                                "The .txt file should have two columns: lmin lmax, which are the limits of the continuum regions in Angstroms. " \
                                                                                "If not provided, only the inverse variance extension of the cube will be considered as noise for simulations.")
    parser.add_argument('-z', '--redshift', type=float, default=None, help='Object redshift. Required if continuum regions file is given.')
    parser.add_argument('--debug', type=int, default=0, help="Debug level. 0: no debug, 1: print more detailed intermediate results. Default 0.")
    parser.add_argument('--seed', type=int, default=12345, help='Seed to initialize simulations.')
    args = parser.parse_args()

    print(f"\033[1m Welcome to weavetab \033[0m")

    cube_file = args.cube_file
    number_simulations = args.number_simulations
    save_tab = args.save_table
    cont_regions_file = args.cont_regions
    redshift = args.redshift
    debug_level = args.debug
    seed = args.seed

    if args.region is None:
        # process full cube
        region = None
        info("No region given. Processing the complete cube")
    else:
        region = args.region
        info(f"Region of interest specified: {region}")
        x1, x2, y1, y2 = region
        if x2 <= x1 or y2 <= y1:
            raise ValueError("Invalid region: ensure x2 > x1 and y2 > y1.")
        if any(v < 0 for v in region):
            raise ValueError("Region indices must be non-negative.")

    # cube filename must have a .fit, .fits, .FIT, or .FITS extension
    if not cube_file.lower().endswith(('.fit', '.fits')):
        raise ValueError("Cube file must have a FITS extension (.fit/.fits).")
    
    cube_path = Path(cube_file).resolve()
    if not cube_path.is_file():
        raise FileNotFoundError(f"Cube file '{cube_path}' does not exist.")

    # Extracting root for working directory
    working_dir = Path.cwd()
    info(f"Current working directory: {working_dir}")
    info(f"Processing cube file: {cube_path}")

    # Number of simulations: has to be positive number
    if number_simulations < 0:
        raise ValueError("Number of simulations must be a non-negative integer.")
    
    if number_simulations == 0:
        info("No simulations will be performed. Only the original spectra will be extracted.")
    else:
        info(f"Number of simulations: {number_simulations}")
        info(f"A total number of {number_simulations+1} FITS cubes will be created.")

    # Save or not a table with results of processing
    if save_tab:
        info("A FITS table will be saved containing spatial coordinates and spectrum per pixel.")
    else:
        info("No FITS table will be saved. Only a FITS cube.")

    # continuum regions of spectra to consider to scale sigma for simulation noise
    #xfits, yfits = None, None
    if cont_regions_file is not None:
        cont_regions_path = working_dir / cont_regions_file
        if not cont_regions_path.is_file():
            raise FileNotFoundError(f"Continuum regions file '{cont_regions_path}' does not exist in the working directory")
        
        info(f"Continuum regions file provided: {cont_regions_path}")

        # redshift
        if redshift is None:
            raise ValueError("Object redshift required")
        else:
            info(f"Object redshift is {redshift}")

        cont_regions = read_continuum_file(cont_regions_path, redshift)
        info(f'Loaded continuum regions from {cont_regions_file}: {cont_regions}')

        #if xfits is None or yfits is None:
        #   raise ValueError("No (x, y) coordinates provided in file")
        #info(f'Using spaxel (x={xfits}, y={yfits}) for scaling factor computation.')
    else:
        info("No continuum regions file provided. No scaling will be applied to inverse invariance for simulation noise")
        cont_regions = None

    # Debug level
    if debug_level < 0 or debug_level > 1:
        raise ValueError(f"Invalid debug number of '{debug_level}': must be 0 or 1")
    
    # seed
    if args.seed != 12345:
        if not (0 <= seed < 2**10):
            raise ValueError("Seed must be between 0 and 2**10 - 1")
        seed = args.seed
        info(f"Using provided seed: {seed}")
    else:
        info(f"No seed provided. Default seed is: {seed}")
        

    extract_cube_spectra(working_dir, cube_path, number_simulations, debug_level, start, region, save_tab, cont_regions, seed)


if __name__ == "__main__":
    main()