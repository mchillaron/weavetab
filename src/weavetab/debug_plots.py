#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np
import matplotlib.pyplot as plt

def plot_rms_vs_sigma(lambda_c, rms_robust, rms_std, sigma_mean, sigma_median, output_dir, debug_level):
    plt.figure()

    plt.scatter(lambda_c, rms_robust, marker="D", label="RMS (robust)")
    plt.scatter(lambda_c, rms_std, marker="s", label="RMS (std)")
    plt.scatter(lambda_c, sigma_mean, marker="*", label="Sigma mean")
    plt.scatter(lambda_c, sigma_median, marker="o",label="Sigma median")
    for x, y1, y2 in zip(lambda_c, rms_robust, sigma_median):
        plt.plot([x, x], [y1, y2], color="gray", linewidth=1, alpha=0.6)

    plt.xlabel("Wavelength")
    plt.ylabel("Value")
    plt.legend()
    plt.title("RMS vs Sigma")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(f"{output_dir}/debug_rms_vs_sigma.pdf")
    if debug_level == 1:
        plt.show()
    plt.close()

def plot_scaling_function(wave, lambda_c, factors, scaling_func, output_dir, debug_level):
    plt.figure()

    plt.scatter(lambda_c, factors, label="Measured factors")

    wave_dense = np.linspace(wave.min(), wave.max(), 1000)
    plt.plot(wave_dense, scaling_func(wave_dense), label="Interpolated")

    plt.xlabel("Wavelength")
    plt.ylabel("Scaling factor")
    plt.legend()
    plt.title("Sigma rescaling function")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(f"{output_dir}/debug_scaling_function.pdf")
    if debug_level == 1:
        plt.show()
    plt.close()

def plot_spectrum_with_regions(wave, flux, cont_regions, output_dir, debug_level):
    plt.figure(figsize=(10,5))

    plt.plot(wave, flux, color='black', lw=1)

    for (lmin, lmax) in cont_regions:
        plt.axvspan(lmin, lmax, alpha=0.2)

    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Continuum regions selection")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(f"{output_dir}/debug_regions.pdf")
    if debug_level == 1:
        plt.show()
    plt.close()

def plot_simulations(wave, flux, simulated_fluxes, output_dir, n_show=20):
    plt.figure(figsize=(10,5))

    # original
    plt.plot(wave, flux, color='black', lw=2, label="Original")

    # subset simulaciones
    for i in range(min(n_show, simulated_fluxes.shape[0])):
        plt.plot(wave, simulated_fluxes[i], alpha=0.2)

    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Simulated spectra spread")
    plt.legend()
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(f"{output_dir}/debug_simulations.pdf")
    plt.show()
    plt.close()

def plot_spectrum_with_regions_and_info(
    wave, 
    flux, 
    cont_regions, 
    lambda_centers, 
    factors,
    output_dir,
    debug_level,
    rms_list=None,
    sigma_list=None
):

    plt.figure(figsize=(12,5))
    plt.plot(wave, flux, color='black', lw=1)
    y_text = np.nanpercentile(flux, 95)

    for i, (lmin, lmax) in enumerate(cont_regions):

        # color shade of the region
        plt.axvspan(lmin, lmax, alpha=0.2, color='green' if factors[i] < 1.2 else 'red')

        if i < len(lambda_centers):
            xc = lambda_centers[i]
            factor = factors[i]
            text = f"{factor:.2f}"

            # extra info
            if rms_list is not None and sigma_list is not None:
                text += f"\nRMS={rms_list[i]:.2e}"
                text += f"\nσ={sigma_list[i]:.2e}"

            # dibujar texto
            plt.text(
                xc, 
                y_text, 
                text,
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=90,
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')
            )

    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Continuum regions + scaling factors")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()

    plt.tight_layout()
    plt.savefig(f"{output_dir}/debug_regions_with_factors.pdf")
    if debug_level == 1:
        plt.show()
    plt.close()