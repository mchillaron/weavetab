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

from .factor_sigma_rms import process_all_cube, rescale_sigma_global, rescale_sigma_per_spectrum
from .debug_plots import plot_spectrum_with_regions, plot_spectrum_with_regions_and_info, plot_rms_sigma_and_factors, plot_simulations, plot_distributions, plot_factor_map

def info(msg):
    GREEN   = "\033[92m"
    RESET   = "\033[0m"
    print(f"{GREEN}INFO:{RESET} {msg}")

def simulate_cubes(new_cube, wavelength, output_dir, region, number_simulations, 
                   cont_regions, seed, debug_level):
    """
    Generate Monte Carlo simulated data cubes from an input cube using Gaussian noise.

    Parameters
    ----------
    new_cube : fits.HDUList
        Input cube containing:
        - Primary HDU (header)
        - DATA extension (flux cube)
        - SIGMA extension (noise cube)
    
    wavelength : array-like
        Wavelength axis corresponding to the cube.

    output_dir : pathlib.Path
        Directory where simulated cubes will be stored.

    region : tuple or None
        Region used for naming outputs.

    number_simulations : int
        Number of simulations to generate in this run.

    cont_regions : array-like or None
        Continuum wavelength regions used to compute a scaling factor for the sigma cube.
        If provided, the scaling is computed from the *collapsed spectrum* of the cube.

    seed : int
        Base random seed ensuring reproducibility.

    debug_level : int
        Debug flag (0 = off, 1 = verbose with plots).

    Notes
    -----
    - Simulations are generated as:
        simulated_cube ~ Normal(mean=data_cube, sigma=sigma_cube)
    - Random numbers are generated using NumPy's Generator (default_rng).
    - The sigma cube can optionally be rescaled using continuum regions.
    - This function supports full reproducibility and resumable simulations. Each simulation
    is deterministically associated with a unique random seed derived from a base seed
    and the simulation index. This guarantees that:
    1. Running the same command with the same seed reproduces identical simulations.
    2. Additional simulations can be appended without altering previous ones.
    3. Deleted simulations can be regenerated exactly.
    """

    header_total = new_cube[0].header
    data_cube = new_cube[1].data
    sigma_cube = new_cube[2].data

    nw, ny, nx = data_cube.shape
    info(f"Starting {number_simulations} simulations...")
    print(f"Cube dimensions: (nw={nw}, ny={ny}, nx={nx})")

    # create a new directory in output_dir to store the simulated cubes
    new_cube_path = output_dir / "simulated_cubes"
    print("Creating directory for simulated cubes...")
    new_cube_path.mkdir(parents=True, exist_ok=True)

    # Check if the directory is empty and alert the user if it is not, that files will be overwritten.
    if any(new_cube_path.iterdir()):
        info(f'Warning: The directory {new_cube_path} is not empty. Existing files may be overwritten.')
    else:
        info(f'The directory {new_cube_path} is empty. Safe to proceed.')

    # Check which was the last simulation number in the directory and ask the user if they want 
    # to continue from that number or overwrite all simulations.
    # seed will be adapted depending on the case

    existing_simulations = sorted(new_cube_path.glob("*.fits"))
    if existing_simulations:
        last_simulation = existing_simulations[-1].stem
        match = re.search(r"_(\d{4})$", last_simulation)

        if match:
            last_number = int(match.group(1))
            info(f'Found existing simulations up to {last_number:04d}.')
            user_input = input(f"Continue from {last_number + 1:04d}? (y/n): ")

            if user_input.lower() == 'y':
                start_number = last_number + 1
            else:
                start_number = 1
                info("Overwriting previous simulations.")
        else:
            info('No valid simulation files found. Starting from simulation 0001.')
            start_number = 1
    else:
        info(f'No existing simulations found. Starting from simulation 0001.')
        start_number = 1

    # Region naming
    if region is None:
        region_tag = "cube_spectra_simul"
    else:
        # the x and y indexes are fits format
        x1, x2, y1, y2 = region
        region_tag = f"cube_spectra_{x1}-{x2}_{y1}-{y2}_simul"

    # Rescaling the sigma cube using factor_sigma_rms.py on collapsed spectrum
    if cont_regions is not None:
    
        info("Plotting the continuum regions over collapsed spectrum")
        flux_total = np.nanmean(data_cube, axis=(1, 2))
        sigma_total = np.nanmean(sigma_cube, axis=(1, 2))

        plot_spectrum_with_regions(wavelength, flux_total, cont_regions, new_cube_path, debug_level)

        info("Computing scaling for all spectra in the cube")
        coords, lambda_centers, all_factors, all_rms_robust, all_rms_std, all_sigma_median, all_sigma_mean = process_all_cube(data_cube, sigma_cube, wavelength, 
                                                                                                                    cont_regions, debug_level=debug_level)

        plot_distributions(all_rms_robust, all_rms_std, all_sigma_mean,
                       all_sigma_median, all_factors, new_cube_path)
        
        global_median_factors = plot_rms_sigma_and_factors(lambda_centers, all_rms_robust, all_rms_std, all_sigma_mean, all_sigma_median, all_factors, 
                                  wavelength, flux_total, cont_regions, new_cube_path, debug_level)

        plot_spectrum_with_regions_and_info(wavelength, flux_total, cont_regions, lambda_centers,
                                            all_factors, all_rms_robust, all_sigma_median,
                                            new_cube_path, debug_level)
        
        # option 1: global rescaling with only one median factor for the whole cube
        # sigma_cube, global_factor = rescale_sigma_global(sigma_cube, all_factors, debug_level)

        # option 2: rescaling spaxel by spaxel with its median factor
        sigma_cube, factors_map = rescale_sigma_per_spectrum(sigma_cube, all_factors, coords)
        plot_factor_map(factors_map, new_cube_path, debug_level)

    else:
        print('No sigma rescaling applied. Simulations will be based on the original sigma cube.')

    # setting the seed properly depending on the case:
    total_needed = start_number + number_simulations - 1
    base_seq = np.random.SeedSequence(seed)
    child_seeds = base_seq.spawn(total_needed) 

    for i in tqdm(range(start_number, start_number + number_simulations), desc="Simulating WEAVE cubes"): #bar_format=bar_format,
        
        rng = np.random.default_rng(child_seeds[i - 1])  # index from 0

        # Draw random cube: Normal(mean=data, sigma=sigma)
        simulated_cube = rng.normal(loc=data_cube, scale=sigma_cube)

        header = header_total.copy()
        header["SEED"] = seed     # save seed in header
        header["SIM"] = i         # save simulation index

        primary_hdu = fits.PrimaryHDU(header=header)
        data_hdu = fits.ImageHDU(data=simulated_cube, name="DATA")
        sigma_hdu = fits.ImageHDU(data=sigma_cube, name="SIGMA")

        hdul = fits.HDUList([primary_hdu, data_hdu, sigma_hdu])

        sim_name = f"{region_tag}_{i:04d}.fits"
        output_path = new_cube_path / sim_name
        hdul.writeto(output_path, overwrite=True)

        print(f"Simulation {i-start_number+1}/{number_simulations} saved to {output_path}")

    
    # Plot with random simulations
    y0, x0 = sigma_cube.shape[1] // 2, sigma_cube.shape[2] // 2
    flux_spec = data_cube[:, y0, x0]
    sigma_spec = sigma_cube[:, y0, x0]

    sims = []
    for j in range(20):
        sim = rng.normal(loc=flux_spec, scale=sigma_spec)
        sims.append(sim)

    plot_simulations(wavelength, flux_spec, np.array(sims), new_cube_path, debug_level)

    info("All simulations completed.")
