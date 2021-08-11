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
import seaborn as sns


def main():
    plt.rcParams['figure.autolayout'] = True
    # raw_mix_comp()
    raw_dists()
    print('donzo')


def raw_dists():
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
    fig_kde, ax_kde = plt.subplots()
    ax.axvline(0.5, color='black', ls='dotted', alpha=0.4)
    ax_norm.axvline(0.5, color='black', ls='dotted', alpha=0.4)
    colors = dict(zip(divs, ['b', 'g', 'r', 'c', 'm', 'y', 'k'][:len(divs)]))

    for div in divs:
        path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

        az_data = AzimuthBinData(div=div, path=path)
        ratio_hist = az_data.get_ratio_dist(ratio_binning)
        ax.bar(ratio_bin_centers, ratio_hist, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5, label=f'{div}$^\circ$')
        ratio_hist_norm = az_data.get_ratio_dist(ratio_binning, True)
        ax_norm.bar(ratio_bin_centers, ratio_hist_norm, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5,
                    label=f'{div}$^\circ$')
        ratio_unbin = az_data.get_ratio_dist()
        values = []
        for x, num in ratio_unbin.items():
            for i in range(num):
                values.append(x)
        sns.histplot(ax=ax_kde, x=values, stat='probability', bins=ratio_binning, kde=True, color=colors[div])
        # sns.rugplot(ax=ax_kde, x=values)

    fig_kde.legend(labels=[f'{div}$^\circ$' for div in divs])
    fig.legend()
    ax.grid()
    fig_norm.legend()
    ax_norm.grid()
    plt.show()


def raw_mix_comp():
    base_raw_path = '/home/dylan/Research/Data/default/'
    base_mix_path = '/home/dylan/Research/Data_Mix/default/'
    set_name = 'rapid05_n1ratios_dca1_nsprx1_m2r6_m2s0_nhfit20_'
    set_num = 0
    energy = 7
    divs = [120, 180, 240]
    cent = 8
    ratio_bin_step = 0.02
    ratio_binning = np.arange(0 - ratio_bin_step / 2, 1 + ratio_bin_step * 3 / 2, ratio_bin_step)
    ratio_bin_centers = (ratio_binning[1:] + ratio_binning[:-1]) / 2

    for div in divs:
        fig_diff, ax_diff = plt.subplots()
        fig_norm, ax_norm = plt.subplots()
        fig_diff_div, ax_diff_div = plt.subplots()
        fig_kde, ax_kde = plt.subplots()
        ax_diff.axhline(0, color='red', ls='--', alpha=1.0)
        ax_diff_div.axhline(0, color='red', ls='--', alpha=1.0)
        ax_norm.axvline(0.5, color='black', ls='dotted', alpha=0.4, label='0.5')
        ax_norm.axvline(div / 360, color='red', ls='dashed', alpha=0.8, label='mean')

        hists = {}
        for data, data_path in [('raw', base_raw_path), ('mix', base_mix_path)]:
            path = f'{data_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
            az_data = AzimuthBinData(div=div, path=path)
            ratio_hist_norm = az_data.get_ratio_dist(ratio_binning, True)
            hists.update({data: np.asarray(ratio_hist_norm)})
            ax_norm.bar(ratio_bin_centers, ratio_hist_norm, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5,
                        label=f'{data} {div}$^\circ$')
            ratio_unbin = az_data.get_ratio_dist()
            values = []
            for x, num in ratio_unbin.items():
                for i in range(num):
                    values.append(x)
            sns.kdeplot(ax=ax_kde, x=values)

        ax_diff.scatter(ratio_bin_centers, hists['raw'] - hists['mix'], label=f'Raw - Mixed {div}$^\circ$')
        diff_div = [2 * (r - m) / (r + m) if r > 0 and m > 0 else float('nan')
                    for r, m in zip(hists['raw'], hists['mix'])]
        ax_diff_div.scatter(ratio_bin_centers, diff_div, label=f'2*(Raw - Mixed) / (Raw + Mixed) {div}$^\circ$')

        ax_diff.legend()
        ax_diff.grid()
        ax_diff_div.grid()
        ax_diff_div.legend()
        ax_norm.legend()
        ax_norm.grid()
    plt.show()


if __name__ == '__main__':
    main()
