#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 28 12:50 PM 2021
Created in PyCharm
Created as QGP_Scripts/dists_with_bin_width

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from AzimuthBinData import AzimuthBinData


def main():
    base_path = '/home/dylan/Research/Data/default/'
    set_name = 'rapid05_n1ratios_dca1_nsprx1_m2r6_m2s0_nhfit20_'
    set_num = 0
    energy = 7
    divs = [60, 180, 300]
    cent = 8
    ratio_bin_step = 0.02
    ratio_binning = np.arange(0 - ratio_bin_step / 2, 1 + ratio_bin_step * 3 / 2, ratio_bin_step)
    ratio_bin_centers = (ratio_binning[1:] + ratio_binning[:-1]) / 2

    fig, ax = plt.subplots()
    fig_norm, ax_norm = plt.subplots()
    ax.axvline(0.5, color='black', ls='dotted', alpha=0.4)
    ax_norm.axvline(0.5, color='black', ls='dotted', alpha=0.4)

    for div in divs:
        path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

        az_data = AzimuthBinData(div=div, path=path)
        ratio_hist = az_data.get_ratio_dist(ratio_binning)
        ax.bar(ratio_bin_centers, ratio_hist, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5, label=f'{div}$^\circ$')
        ratio_hist_norm = az_data.get_ratio_dist(ratio_binning, True)
        ax_norm.bar(ratio_bin_centers, ratio_hist_norm, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5,
                    label=f'{div}$^\circ$')

    fig.legend()
    ax.grid()
    fig_norm.legend()
    ax_norm.grid()
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
