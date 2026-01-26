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

GREEN   = "\033[92m"
CYAN    = "\033[96m"
RESET   = "\033[0m"

def simulate_cubes(new_cube, new_cube_path, region, number_simulations):
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

    # Region naming
    if region is None:
        region_tag = "cube_spectra_simul"
    else:
        x1_fits, x2_fits, y1_fits, y2_fits = region
        region_tag = f"cube_spectra_{x1_fits}-{x2_fits}_{y1_fits}-{y2_fits}_simul"

    for i in tqdm(range(1, number_simulations + 1), desc="Simulating WEAVE cubes"): #bar_format=bar_format,

        # Draw random cube: Normal(mean=data, sigma=sigma)
        simulated_cube = np.random.normal(loc=data_cube, scale=sigma_cube)

        # Build the new cube
        primary_hdu = fits.PrimaryHDU(header=header_total)
        data_hdu = fits.ImageHDU(data=simulated_cube, name="DATA")
        hdul = fits.HDUList([primary_hdu, data_hdu])

        sim_name = f"{region_tag}_{i:04d}.fits"
        output_path = new_cube_path / sim_name
        hdul.writeto(output_path, overwrite=True)

        print(f"Simulation {i}/{number_simulations} saved to {output_path}")

    print(f"{GREEN}INFO:{RESET} All simulations completed.")
