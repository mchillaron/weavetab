#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
from scipy.interpolate import interp1d

# robust rms based on MAD
def robust_rms(flux):
    median = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - median))
    rms_std = np.nanstd(flux)
    rms_robust = 1.4826 * mad
    print(f"robust_rms: rms_std={rms_std}, rms_robust={rms_robust}")
    return rms_robust, rms_std

def compute_scaling_factors(wave, flux, sigma, cont_regions):
    print("\n[compute_scaling_factors] Starting...")

    factors = []
    lambda_centers = []
    rms_robust_list = []
    rms_std_list = []
    sigma_mean_list = []
    sigma_median_list = []

    print(f"Input shapes -> wave: {wave.shape}, flux: {flux.shape}, sigma: {sigma.shape}")

    for i, (lmin, lmax) in enumerate(cont_regions):
        print(f"\nRegion {i+1}: {lmin} - {lmax}")

        mask = (wave >= lmin) & (wave <= lmax)
        n_pixels = np.sum(mask)

        print(f"Pixels in region: {n_pixels}")

        if n_pixels < 5:
            print("WARNING: Too few pixels, skipping region")
            continue

        flux_region = flux[mask]
        sigma_region = sigma[mask]

        print(f"flux_region shape: {flux_region.shape}")
        print(f"sigma_region shape: {sigma_region.shape}")

        rms_robust, rms_std = robust_rms(flux_region)
        sigma_mean = np.nanmean(sigma_region)
        sigma_median = np.nanmedian(sigma_region)

        print(f"RMS (robust): {rms_robust}")
        print(f"RMS (standard): {rms_std}")
        print(f"sigma_mean: {sigma_mean}")
        print(f"sigma_median: {sigma_median}")

        if sigma_median > 0:
            factor = rms_robust / sigma_median
            print(f"  Scaling factor: {factor}")

            factors.append(factor)
            lambda_centers.append(0.5 * (lmin + lmax))

            # save for debug plots
            rms_robust_list.append(rms_robust)
            rms_std_list.append(rms_std)
            sigma_mean_list.append(sigma_mean)
            sigma_median_list.append(sigma_median)
        else:
            print("WARNING: sigma_mean <= 0, skipping")

    print("\n[compute_scaling_factors] Done.\n")

    return (
        np.array(lambda_centers),
        np.array(factors),
        np.array(rms_robust_list),
        np.array(rms_std_list),
        np.array(sigma_mean_list),
        np.array(sigma_median_list),
    )


def build_scaling_function(lambda_centers, factors, kind='linear'):
    return interp1d(
        lambda_centers,
        factors,
        bounds_error=False,
        fill_value="extrapolate",
        kind=kind
    )


def rescale_sigma(wave, sigma, scaling_func):
    print("\n[rescale_sigma] Starting sigma rescaling...")

    factor_lambda = scaling_func(wave)
    print(f"[rescale_sigma] factor_lambda shape: {factor_lambda.shape}")

    # Expand dimensions for cube broadcasting
    factor_lambda = factor_lambda[:, None, None]
    print(f"[rescale_sigma] reshaped factor_lambda: {factor_lambda.shape}")
    print(f"[rescale_sigma] sigma shape: {sigma.shape}")

    sigma_rescaled = sigma * factor_lambda

    print("[rescale_sigma] Rescaling done.\n")
    return sigma_rescaled