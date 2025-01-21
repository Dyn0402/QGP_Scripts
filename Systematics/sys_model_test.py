#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 20 3:29 PM 2025
Created in PyCharm
Created as QGP_Scripts/sys_model_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    # Define systematic parameter ranges
    sys_def = 1.0
    sys_range = 0.5
    # variation_range = 0.5
    variation_range = sys_range * 2. / 3.
    # prob_func = 'gaussian'
    prob_func = 'uniform'
    num_samples = 10000

    def val_func(sys):
        return 2 * sys - 4

    range_mult = 2  # Go this many ranges on each side of default for plotting purposes

    # Plot systematic parameter distribution
    sys_vals = np.linspace(sys_def - sys_range * range_mult, sys_def + sys_range * range_mult, 1000)
    if prob_func == 'gaussian':
        sys_probs = np.exp(-0.5 * ((sys_vals - sys_def) / sys_range) ** 2) / (sys_range * np.sqrt(2 * np.pi))
    elif prob_func == 'uniform':
        sys_probs = np.where(np.abs(sys_vals - sys_def) < sys_range, 1 / (2 * sys_range), 0)
    else:
        sys_probs = np.zeros_like(sys_vals)
    # sys_probs /= np.sum(sys_probs)

    fig, ax = plt.subplots()
    ax.axvline(sys_def, color='green', linestyle='--', label='Default Parameter')
    ax.axvline(sys_def - variation_range, color='orange', linestyle='--', label='Variation Low')
    ax.axvline(sys_def + variation_range, color='red', linestyle='--', label='Variation High')
    ax.plot(sys_vals, sys_probs)
    ax.set_xlabel('Systematic Parameter Value')
    ax.set_ylabel('Probability Density')
    ax.set_title('Systematic Parameter Distribution')
    ax.set_ylim(bottom=0)
    ax.set_xlim(sys_def - sys_range * range_mult, sys_def + sys_range * range_mult)
    ax.legend()
    fig.tight_layout()

    # Plot val_func
    fig_val_func, ax_val_func = plt.subplots()
    # ax_val_func.axvline(sys_def, color='green', linestyle='--', label='Default Parameter')
    # ax_val_func.axhline(val_func(sys_def), color='orange', linestyle='--', label='Default Value')
    ax_val_func.plot(sys_vals, val_func(sys_vals), zorder=1)
    ax_val_func.plot(sys_def, val_func(sys_def), color='green', marker='^', markersize=10, label='Default Parameter', zorder=5)
    ax_val_func.plot(sys_def - variation_range, val_func(sys_def - variation_range), color='orange', marker='^', markersize=10, label='Variation Low', zorder=5)
    ax_val_func.plot(sys_def + variation_range, val_func(sys_def + variation_range), color='red', marker='^', markersize=10, label='Variation High', zorder=5)
    ax_val_func.set_xlabel('Systematic Parameter Value')
    ax_val_func.set_ylabel('Value')
    ax_val_func.set_title('Value Function')
    ax_val_func.set_xlim(sys_def - sys_range * range_mult, sys_def + sys_range * range_mult)
    ax_val_func.set_ylim(bottom=val_func(sys_def - sys_range * range_mult), top=val_func(sys_def + sys_range * range_mult))
    ax_val_func.legend()
    fig_val_func.tight_layout()

    # Sample sys parameters from distribution and calculate values
    if prob_func == 'gaussian':
        sys_samples = np.random.normal(sys_def, sys_range, num_samples)
    elif prob_func == 'uniform':
        sys_samples = np.random.uniform(sys_def - sys_range, sys_def + sys_range, num_samples)
    else:
        print('Invalid probability function')
        return
    val_samples = val_func(sys_samples)
    # Get the 1 sigma percentile locations
    val_1_sigma_low = np.percentile(val_samples, 16)
    val_1_sigma_high = np.percentile(val_samples, 84)

    # Calculate systematic error as I would
    sys_err = np.abs(val_func(sys_def + variation_range) - val_func(sys_def))

    # Plot value distribution
    fig_val_dist, ax_val_dist = plt.subplots()
    ax_val_dist.hist(val_samples, bins=100, density=True, orientation='horizontal')
    ax_val_dist.plot(0.1, val_func(sys_def), color='yellow', marker='o', markersize=10, zorder=5, label='Original')
    ax_val_dist.errorbar(0.1, val_func(sys_def), xerr=0, yerr=sys_err, color='yellow', marker='', alpha=0.4, elinewidth=5, zorder=5)
    ax_val_dist.plot(0.15, val_func(sys_def), color='orange', marker='o', markersize=10, zorder=5,
                     label=r'$2/\sqrt{12}$')
    ax_val_dist.errorbar(0.15, val_func(sys_def), xerr=0, yerr=sys_err * 2 / np.sqrt(12), color='orange', marker='',
                         alpha=0.4, elinewidth=5, zorder=5)
    ax_val_dist.plot(0.2, val_func(sys_def), color='pink', marker='o', markersize=10, zorder=5, label=r'$1/\sqrt{12}$')
    ax_val_dist.errorbar(0.2, val_func(sys_def), xerr=0, yerr=sys_err / np.sqrt(12), color='pink', marker='',
                         alpha=0.4, elinewidth=5, zorder=5)
    ax_val_dist.axhline(val_1_sigma_low, color='red', linestyle='--', label='68% of Values', zorder=5)
    ax_val_dist.axhline(val_1_sigma_high, color='red', linestyle='--', zorder=5)
    ax_val_dist.set_ylabel('Value')
    ax_val_dist.set_xlabel('Probability Density')
    ax_val_dist.set_title('Value Distribution')
    ax_val_dist.legend()
    fig_val_dist.tight_layout()

    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
