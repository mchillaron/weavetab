#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from astropy.io import fits
from tqdm import tqdm

import numpy as np
import re

GREEN   = "\033[92m"
CYAN    = "\033[96m"
RESET   = "\033[0m"

def simulate_cubes(new_cube, output_dir, region, number_simulations):
    """
    Generate Monte-Carlo simulated cubes based on the cube created by table_to_cube().

    Parameters
    ----------
    new_cube : fits.HDUList
        HDUList containing Primary, DATA, SIGMA extensions.
    new_cube_path : str or Path
        Path of the reconstructed cube (for naming consistency).
    cube_filename : str
        Base name of the cube.
    working_dir : str
        Directory where simulated cubes will be written.
    region : tuple or None
        Region used in table_to_cube (for naming).
    number_simulations : int
        Number of simulated cubes to generate.
    """

    header_total = new_cube[0].header
    data_cube = new_cube[1].data
    sigma_cube = new_cube[2].data

    nw, ny, nx = data_cube.shape
    print(f"{GREEN}INFO:{RESET} Starting {number_simulations} simulations...")
    print(f"Cube dimensions: (nw={nw}, ny={ny}, nx={nx})")

    # create a new directory in output_dir to store the simulated cubes
    new_cube_path = output_dir / "simulated_cubes"
    print("Creating directory for simulated cubes...")
    new_cube_path.mkdir(parents=True, exist_ok=True)
    # Check if the directory is empty and alert the user if it is not, that files will be overwritten.
    if any(new_cube_path.iterdir()):
        print(f'{GREEN}INFO:{RESET} Warning: The directory {new_cube_path} is not empty. Existing files may be overwritten.')
    else:
        print(f'{GREEN}INFO:{RESET} The directory {new_cube_path} is empty. Safe to proceed.')

    # Check which was the last simulation number in the directory and ask the user if they want to continue from that number or overwrite all simulations.
    existing_simulations = sorted(new_cube_path.glob("*.fits"))
    if existing_simulations:
        last_simulation = existing_simulations[-1].stem
        match = re.search(r"_(\d{4})$", last_simulation)
        if match:
            last_number = int(match.group(1))
            print(f'{GREEN}INFO:{RESET} Found existing simulations up to number {last_number:04d}.')
            user_input = input(f"Do you want to continue from simulation {last_number + 1:04d}? (y/n): ")
            if user_input.lower() == 'y':
                start_number = last_number + 1
            else:
                start_number = 1
                print(f'{GREEN}INFO:{RESET} All existing simulations will be overwritten.')
        else:
            print(f'{GREEN}INFO:{RESET} No valid simulation files found. Starting from simulation 0001.')
            start_number = 1
    else:
        print(f'{GREEN}INFO:{RESET} No existing simulations found. Starting from simulation 0001.')
        start_number = 1

    # Region naming
    if region is None:
        region_tag = "cube_spectra_simul"
    else:
        x1_fits, x2_fits, y1_fits, y2_fits = region
        region_tag = f"cube_spectra_{x1_fits}-{x2_fits}_{y1_fits}-{y2_fits}_simul"

    for i in tqdm(range(start_number, start_number + number_simulations), desc="Simulating WEAVE cubes"): #bar_format=bar_format,
        print(i)
        # Draw random cube: Normal(mean=data, sigma=sigma)
        simulated_cube = np.random.normal(loc=data_cube, scale=sigma_cube)

        # Build the new cube
        primary_hdu = fits.PrimaryHDU(header=header_total)
        data_hdu = fits.ImageHDU(data=simulated_cube, name="DATA")
        hdul = fits.HDUList([primary_hdu, data_hdu])

        sim_name = f"{region_tag}_{i:04d}.fits"
        output_path = new_cube_path / sim_name
        hdul.writeto(output_path, overwrite=True)

        print(f"Simulation {i-start_number+1}/{number_simulations} saved to {output_path}")

    print(f"{GREEN}INFO:{RESET} All simulations completed.")
