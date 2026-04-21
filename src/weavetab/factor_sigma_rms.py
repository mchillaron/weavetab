#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
import teareduce as tea

from scipy.interpolate import interp1d

def process_all_cube(data_cube, sigma_cube, wavelength, cont_regions, debug_level):
    """
    Compute scaling statistics across all spaxels in the cube.
    Returns mean, median, standard and robust deviation of all quantities per region.

    Returns
    -------
    lambda_centers : ndarray
    mean_factors : ndarray
    std_factors : ndarray
    mean_rms : ndarray
    std_rms : ndarray
    mean_sigma : ndarray
    std_sigma : ndarray
    """

    nw, ny, nx = data_cube.shape

    coords = []
    all_lambda_centers = []
    all_factors = []
    all_rms_robust = []
    all_rms_std = []
    all_sigma_median = []
    all_sigma_mean = []

    for y in range(ny):
        for x in range(nx):

            flux_spec = data_cube[:, y, x]
            sigma_spec = sigma_cube[:, y, x]

            if np.all(np.isnan(flux_spec)):
                continue

            lambda_c, factors, rms_rob, rms_std, sig_mean, sig_med = compute_scaling_factors(
                wavelength,
                flux_spec,
                sigma_spec,
                cont_regions,
                debug=debug_level
            )

            if len(factors) == 0:
                continue
            
            coords.append((y, x)) # python indexes
            all_lambda_centers.append(lambda_c)
            all_factors.append(factors)
            all_rms_robust.append(rms_rob)
            all_rms_std.append(rms_std)
            all_sigma_median.append(sig_med)
            all_sigma_mean.append(sig_mean)

    coords = np.array(coords)
    all_lambda_centers = np.array(all_lambda_centers)
    all_factors = np.array(all_factors)  # shape: (Npix, Nregions)
    all_rms_robust = np.array(all_rms_robust)
    all_rms_std=np.array(all_rms_std)
    all_sigma_median = np.array(all_sigma_median)
    all_sigma_mean = np.array(all_sigma_mean)

    return coords, all_lambda_centers, all_factors, all_rms_robust, all_rms_std, all_sigma_median, all_sigma_mean



def robust_rms(flux, debug=False):
    """
    Compute different RMS estimators for a flux array.

    Parameters
    ----------
    flux : ndarray
        Flux values within a continuum region.
    debug : bool
        If True, print diagnostic values.

    Returns
    -------
    rms_robust : float
        Robust RMS estimate (e.g., using robust_std).
    rms_std : float
        Standard deviation.
    """

    median = np.nanmedian(flux)
    rms_std = np.nanstd(flux)
    rms_robust = tea.robust_std(flux)

    if debug:
        mad = np.nanmedian(np.abs(flux - median))
        rms_mad = 1.4826 * mad
        print(f"rms_std={rms_std}, rms_robust={rms_robust}, rms_mad={rms_mad}")

    return rms_robust, rms_std



def compute_scaling_factors(wave, flux, sigma, cont_regions, debug):
    """
    Compute scaling factors between flux RMS and sigma values in predefined continuum regions.

    For each wavelength region, this function estimates:
    - Robust RMS of the flux
    - Standard RMS
    - Mean and median of the sigma spectrum
    - Scaling factor = robust RMS / sigma_median

    Parameters
    ----------
    wave : array
        Wavelength array.
    flux : array
        Flux spectrum.
    sigma : array
        Sigma (noise) spectrum.
    cont_regions : list of tuples
        Continuum regions defined as (lmin, lmax).
    debug : bool
        Enable verbose output.

    Returns
    -------
    lambda_centers : ndarray
        Central wavelength of each region.
    factors : ndarray
        Scaling factors per region.
    rms_robust_list : ndarray
        Robust RMS values.
    rms_std_list : ndarray
        Standard RMS values.
    sigma_mean_list : ndarray
        Mean sigma values.
    sigma_median_list : ndarray
        Median sigma values.
    """

    factors = []
    lambda_centers = []

    rms_robust_list = []
    rms_std_list = []
    sigma_mean_list = []
    sigma_median_list = []

    for i, (lmin, lmax) in enumerate(cont_regions):

        mask = (wave >= lmin) & (wave <= lmax)
        n_pixels = np.sum(mask)
        #print(f"Pixels in region: {n_pixels}")

        if n_pixels < 5:
            print("WARNING: Too few pixels, skipping region")
            lambda_centers.append(0.5 * (lmin + lmax))
            rms_robust = np.nan
            rms_std = np.nan
            factor = np.nan
            sigma_mean = np.nan
            sigma_median=np.nan
            continue

        flux_region = flux[mask]
        sigma_region = sigma[mask]

        rms_robust, rms_std = robust_rms(flux_region, debug=debug)
        sigma_mean = np.nanmean(sigma_region)
        sigma_median = np.nanmedian(sigma_region)

        if sigma_median == 0:
            sigma_median = np.nan
            factor = np.nan
            continue

        factor = rms_robust / sigma_median

        lambda_centers.append(0.5 * (lmin + lmax))
        factors.append(factor)

        rms_robust_list.append(rms_robust)
        rms_std_list.append(rms_std)
        sigma_mean_list.append(sigma_mean)
        sigma_median_list.append(sigma_median)

    return (
        np.array(lambda_centers),
        np.array(factors),
        np.array(rms_robust_list),
        np.array(rms_std_list),
        np.array(sigma_mean_list),
        np.array(sigma_median_list),
    )



def rescale_sigma_global(sigma_cube, all_factors, debug):
    """
    Rescale the entire sigma cube using a single global factor.

    The factor is computed as:
        median( mean(factors per region) )

    Only applied if factor > 1.

    Parameters
    ----------
    sigma_cube : ndarray (nw, ny, nx)
    all_factors : ndarray (npix, nregions)
    verbose : bool

    Returns
    -------
    sigma_rescaled : ndarray
    global_factor : float
    """

    # mean factor per region
    mean_factors = np.nanmean(all_factors, axis=0)

    # global median
    global_factor = np.nanmedian(mean_factors)

    if debug > 0:
        print(f"INFO: Global scaling factor = {global_factor:.3f}")

    if global_factor > 1:
        sigma_rescaled = sigma_cube * global_factor
        if debug > 0:
            print("INFO: Sigma cube rescaled globally.")
    else:
        sigma_rescaled = sigma_cube
        if debug > 0:
            print("INFO: Global factor <= 1 → no rescaling applied.")

    return sigma_rescaled, global_factor



def rescale_sigma_per_spectrum(sigma_cube, all_factors, coords):
    """
    Rescale sigma cube per spaxel using individual scaling factors.

    For each spectrum:
        factor_i = median( mean(factors across regions) )

    Only spectra with factor_i > 1 are rescaled.

    Parameters
    ----------
    sigma_cube : ndarray (nw, ny, nx)
    all_factors : ndarray (npix, nregions)
    verbose : bool

    Returns
    -------
    sigma_rescaled : ndarray
    factors_map : ndarray (ny, nx)
    """

    nw, ny, nx = sigma_cube.shape
    npix = ny * nx

    # median per region, for each spaxel
    median_factors_per_spaxel = np.nanmedian(all_factors, axis=1)

    factors_map = np.full((ny, nx), np.nan)                         # complete map, initially NaN
    for (y, x), f in zip(coords, median_factors_per_spaxel):        # filling only valid pixels
        factors_map[y, x] = f

    sigma_rescaled = sigma_cube.copy()
    mask = factors_map > 1

    n_rescaled = np.sum(mask)
    n_total = np.sum(~np.isnan(factors_map))

    print(f"INFO: Spaxels rescaled: {n_rescaled}/{n_total} ({100*n_rescaled/n_total:.2f}%)")
    sigma_rescaled[:, mask] *= factors_map[mask]

    return sigma_rescaled, factors_map