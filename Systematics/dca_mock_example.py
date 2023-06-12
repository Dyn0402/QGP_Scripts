#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 08 11:49 AM 2023
Created in PyCharm
Created as QGP_Scripts/dca_mock_example

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import rayleigh


def main():
    mu = 0
    sigma_primary, sigma_secondary = 0.2, 0.5
    secondary_scale = 0.3
    r = np.linspace(0, 3, 1000)
    n_primary = 1000
    binning = np.arange(0, 3, 0.1)
    def_cut, chosen_sys = 1.0, 0.8
    sys_cuts = np.arange(0.1, 3.0, 0.05)
    primary_val, secondary_val = 2, 1

    # primary_dist = rayleigh.pdf(r, mu, sigma_primary)
    primary_dist = rayleigh(mu, sigma_primary)
    primary_pdf = primary_dist.pdf(r)
    # secondary_dist = secondary_scale * rayleigh.pdf(r, mu, sigma_secondary)
    secondary_dist = rayleigh(mu, sigma_secondary)
    secondary_pdf = secondary_scale * secondary_dist.pdf(r)
    total_pdf = primary_pdf + secondary_pdf

    purity = primary_pdf / total_pdf

    plt.figure()
    plt.grid()
    plt.axhline(0, color='black')
    plt.plot(r, primary_pdf, label='primaries')
    plt.plot(r, secondary_pdf, label='secondaries')
    plt.plot(r, total_pdf, label='total')
    plt.legend()

    plt.figure()
    plt.plot(r, purity)

    primary_rvs = primary_dist.rvs(n_primary)
    secondary_rvs = secondary_dist.rvs(int(secondary_scale * n_primary))

    plt.figure()
    plt.hist(primary_rvs, bins=binning, density=False, histtype='step', label='primary')
    plt.hist(secondary_rvs, bins=binning, density=False, histtype='step', label='secondary')
    plt.hist(np.concatenate((primary_rvs, secondary_rvs)), bins=binning, density=False, histtype='step', label='total')
    plt.legend()

    avgs = []
    all_cuts = [def_cut, *sys_cuts]
    for cut in all_cuts:
        prim_under_cut, second_under_cut = np.count_nonzero(primary_rvs < cut), np.count_nonzero(secondary_rvs < cut)
        total_under_cut = prim_under_cut + second_under_cut
        avg_val = (prim_under_cut * primary_val + second_under_cut * secondary_val) / total_under_cut
        avgs.append(avg_val)
        if np.isclose(cut, def_cut):
            def_avg = avg_val
        elif np.isclose(cut, chosen_sys):
            chosen_avg = avg_val
    plt.figure()
    plt.axhline(def_avg, color='black', alpha=0.5, ls='--')
    plt.axhline(chosen_avg, color='red', alpha=0.5, ls='--')
    plt.axhline(primary_val, color='blue', label='primary_val')
    plt.axhline(secondary_val, color='orange', label='secondary_val')
    plt.scatter(all_cuts, avgs, marker='o')

    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
