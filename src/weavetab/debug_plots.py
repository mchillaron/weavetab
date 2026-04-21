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

def plot_spectrum_with_regions(wave, flux, cont_regions, output_dir, debug_level):
    """
    Plot a reference spectrum with highlighted continuum regions.

    Parameters
    ----------
    wave : ndarray
        Wavelength array.
    flux : ndarray
        Flux spectrum (e.g., collapsed cube spectrum).
    cont_regions : list of tuples
        Continuum regions defined as (lmin, lmax).
    output_dir : Path
        Output directory.
    debug_level : int
        If 1, display the plot interactively.
    """

    plt.figure(figsize=(10,5))
    plt.plot(wave, flux, color='black', lw=1)

    for (lmin, lmax) in cont_regions:
        plt.axvspan(lmin, lmax, alpha=0.2)

    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Continuum regions selection")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(output_dir / "debug_regions.pdf")

    if debug_level == 1:
        plt.show() 
    plt.close()



def plot_distributions(all_rms_robust, all_rms_std, all_sigma_mean,
                       all_sigma_median, all_factors, output_dir):

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()

    datasets = [
        ("RMS robust", all_rms_robust.flatten()),
        ("RMS std", all_rms_std.flatten()),
        ("Sigma mean", all_sigma_mean.flatten()),
        ("Sigma median", all_sigma_median.flatten()),
        ("Scaling factor", all_factors.flatten()),
    ]

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

    for i, (title, data) in enumerate(datasets):
        ax = axes[i]

        data = data[~np.isnan(data)]

        mean = np.mean(data)
        median= np.median(data)
        std  = np.std(data)

        bins = np.histogram_bin_edges(data, bins='auto')
        ax.hist(data, bins=bins, color=colors[i])

        ax.axvline(mean, linestyle='--', c='black')
        ax.axvline(median, linestyle='--', c='blue')
        ax.axvline(mean - std, linestyle=':', c='hotpink')
        ax.axvline(mean + std, linestyle=':', c='hotpink')
        ax.set_title(title)

    if len(datasets) < 6:
        axes[-1].axis('off')

        axes[-1].text(
            0.1, 0.7,
            "Dashed line: mean\n"
            "Dotted lines: ±1σ",
            fontsize=10
        )
        #"Mean → dashed\n±1σ → dotted", fontsize=11

    plt.tight_layout()
    plt.savefig(output_dir / "debug_distributions.pdf")
    plt.close()



def plot_rms_sigma_and_factors(
    lambda_c,
    all_rms_robust,
    all_rms_std,
    all_sigma_mean,
    all_sigma_median,
    all_factors,
    wave,
    flux_collapsed,
    cont_regions,
    output_dir,
    debug_level
):
    """
    Two-panel plot:
    Left: RMS vs sigma per region
    Right: scaling factors per region

    Includes:
    - error bars with caps
    - shaded continuum regions
    - background spectrum (faint)
    - mean & median lines (right panel)
    """

    mean_lambda_c = np.nanmean(lambda_c, axis=0)

    mean_rms_rob = np.nanmean(all_rms_robust, axis=0)
    std_rms_rob  = np.nanstd(all_rms_robust, axis=0)

    mean_rms_std = np.nanmean(all_rms_std, axis=0)
    std_rms_std  = np.nanstd(all_rms_std, axis=0)

    mean_sigma_med = np.nanmean(all_sigma_median, axis=0)
    std_sigma_med  = np.nanstd(all_sigma_median, axis=0)

    mean_sigma_mean = np.nanmean(all_sigma_mean, axis=0)
    std_sigma_mean  = np.nanstd(all_sigma_mean, axis=0)

    mean_factors = np.nanmean(all_factors, axis=0)
    std_factors  = np.nanstd(all_factors, axis=0)

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)

    ax1, ax2 = axes

    lmin = min(r[0] for r in cont_regions)
    lmax = max(r[1] for r in cont_regions)

    mask = (wave >= lmin) & (wave <= lmax)

    for ax in axes:
        # continuum bands
        for (l1, l2) in cont_regions:
            ax.axvspan(l1, l2, color='lightgray', alpha=0.3, zorder=0)

    # Subplot 1: rms vs. sigma
    ax1.errorbar(mean_lambda_c, mean_rms_rob, yerr=std_rms_rob, fmt='D', capsize=4, label="RMS (robust)")
    ax1.errorbar(mean_lambda_c, mean_rms_std, yerr=std_rms_std, fmt='s', capsize=4, label="RMS (std)")
    ax1.errorbar(mean_lambda_c, mean_sigma_med, yerr=std_sigma_med, fmt='o', capsize=4, label="Sigma median")
    ax1.errorbar(mean_lambda_c, mean_sigma_mean, yerr=std_sigma_mean, fmt='*', capsize=4, label="Sigma mean")

    ymin1 = np.nanmin([
    mean_rms_rob - std_rms_rob,
    mean_rms_std - std_rms_std,
    mean_sigma_mean - std_sigma_mean,
    mean_sigma_med - std_sigma_med
    ])

    ymax1 = np.nanmax([
        mean_rms_rob + std_rms_rob,
        mean_rms_std + std_rms_std,
        mean_sigma_mean + std_sigma_mean,
        mean_sigma_med + std_sigma_med
    ])

    # spectrum in the background
    # normalize the spectrum to make it visible in the background
    flux_norm1 = flux_collapsed[mask] / np.nanmax(flux_collapsed[mask])
    flux_norm1 = flux_norm1 * (ymax1 - ymin1) + ymin1

    ax1.plot(wave[mask], flux_norm1, color='gray', alpha=0.5, zorder=0)
    margin1 = 0.15 * (ymax1 - ymin1)
    ax1.set_ylim(ymin1 - margin1, ymax1 + margin1)

    # líneas verticales RMS ↔ sigma
    for x, y1, y2 in zip(mean_lambda_c, mean_rms_rob, mean_sigma_med):
        ax1.plot([x, x], [y1, y2], color="gray", alpha=0.5)

    ax1.set_xlabel("Wavelength")
    ax1.set_ylabel("Value")
    ax1.set_title("RMS vs Sigma (mean ± 1σ)")
    ax1.legend()

    # Subplot 2: factors    
    ax2.errorbar(mean_lambda_c, mean_factors, yerr=std_factors, fmt='o', capsize=4, label="Scaling factor")

    ymin2 = np.nanmin(mean_factors - std_factors)
    ymax2 = np.nanmax(mean_factors + std_factors)

    # normalize the spectrum to make it visible in the background
    flux_norm = flux_collapsed[mask] / np.nanmax(flux_collapsed[mask])
    flux_norm = flux_norm * (ymax2 - ymin2) + ymin2

    ax2.plot(wave[mask], flux_norm, color='gray', alpha=0.5, zorder=0)

    margin2 = 0.15 * (ymax2 - ymin2)
    ax2.set_ylim(ymin2 - margin2, ymax2 + margin2)

    # mean and median of global data
    global_mean = np.nanmean(mean_factors)
    global_median = np.nanmedian(mean_factors)

    ax2.axhline(global_mean, linestyle='--', label="Mean")
    ax2.axhline(global_median, linestyle=':', label="Median")

    ax2.set_xlabel("Wavelength")
    ax2.set_ylabel("Scaling factor")
    ax2.set_title("Scaling factors per region")
    ax2.legend()

    for ax in axes:
        ax.tick_params(direction='in', which='both', top=True, right=True)
        ax.minorticks_on()

    plt.tight_layout()
    plt.savefig(output_dir / "debug_rms_sigma_factors.pdf")

    if debug_level == 1:
        plt.show()

    plt.close()
    return global_median



def plot_spectrum_with_regions_and_info(
    wave,
    flux,
    cont_regions,
    all_lambda_centers,
    all_factors,
    all_rms_robust,
    all_sigma_median,
    output_dir,
    debug_level
):

    """
    Plot collapsed spectrum with continuum regions and global scaling statistics.

    This version shows:
    - Mean scaling factor per region
    - 1σ dispersion across all spaxels
    - Optional RMS and sigma statistics
    """

    plt.figure(figsize=(12, 5))
    plt.plot(wave, flux, color='black', lw=1)

    lambda_centers = np.nanmean(all_lambda_centers, axis=0)
    mean_factors = np.nanmean(all_factors, axis=0)
    std_factors = np.nanstd(all_factors, axis=0)

    rms_mean = np.nanmean(all_rms_robust, axis=0)
    rms_std = np.nanstd(all_rms_robust, axis=0)

    sigma_mean = np.nanmean(all_sigma_median, axis=0)
    sigma_std = np.nanstd(all_sigma_median, axis=0)

    y_text = np.nanpercentile(flux, 95)

    for i, (lmin, lmax) in enumerate(cont_regions):

        color = 'green' if mean_factors[i] < 1.0 else 'red'
        plt.axvspan(lmin, lmax, alpha=0.2, color=color)

        xc = lambda_centers[i]

        text = f"{mean_factors[i]:.2f} ± {std_factors[i]:.2f}"

        if rms_mean is not None:
            text += f"\nRMS={rms_mean[i]:.2e} ± {rms_std[i]:.2e}"

        if sigma_mean is not None:
            text += f"\nσ={sigma_mean[i]:.2e} ± {sigma_std[i]:.2e}"

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
    plt.title("Continuum regions + global scaling factors")
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()

    plt.tight_layout()
    plt.savefig(output_dir / "debug_regions_with_factors.pdf")

    if debug_level == 1:
        plt.show()

    plt.close()



def plot_factor_map(factors_map, output_dir, debug_level):
    plt.figure(figsize=(6,5))

    im = plt.imshow(factors_map, origin='lower', aspect='auto')
    plt.colorbar(im, label="Scaling factor")

    plt.title("Scaling factor map (per spaxel)")
    plt.xlabel("X")
    plt.ylabel("Y")

    plt.tight_layout()
    plt.savefig(output_dir / "debug_factor_map.pdf")

    if debug_level == 1:
        plt.show()
    plt.close()


    # Second plot
    factors = factors_map.flatten()
    factors = factors[np.isfinite(factors)]
    mask = factors > 1

    plt.figure()

    bins = np.histogram_bin_edges(factors, bins=50)

    plt.hist(
        [factors[~mask], factors[mask]],
        bins=bins,
        stacked=True,
        label=["≤1", ">1"]
    )

    plt.axvline(1, linestyle='--', color='black', label="Threshold")

    plt.xlabel("Scaling factor")
    plt.ylabel("Density")
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_dir / "debug_factor_hist.pdf")

    if debug_level == 1:
        plt.show()
    plt.close()



def plot_simulations(wave, flux, simulated_fluxes, output_dir, debug_level):
    plt.figure(figsize=(10,5))

    # subset simulaciones
    for i in range(min(20, simulated_fluxes.shape[0])):
        plt.plot(wave, simulated_fluxes[i], alpha=0.2)

    plt.plot(wave, flux, color='black', lw=0.5, label="Original")
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.title("Simulated spectra")
    plt.legend()
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.minorticks_on()
    
    plt.savefig(f"{output_dir}/debug_simulations.pdf")
    if debug_level == 1:
        plt.show()
    plt.close()
