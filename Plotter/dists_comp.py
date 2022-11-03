#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 26 12:12 PM 2021
Created in PyCharm
Created as QGP_Scripts/dists_comp

@author: Dylan Neff, Dyn04
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from AzimuthBinData import AzimuthBinData
from check_azdata_stats import print_stats


def main():
    plt.rcParams['figure.autolayout'] = True
    dists = set_dists()
    show_dists(dists)
    print('donzo')


def set_dists():
    return [dist1(), dist2()]


def dist1():
    dist_type = 'flat80_anticlmulti_spread1_amp01_resample'
    dist_dict = {'base_path': f'F:/Research/Data_Sim/{dist_type}/',
                 'set_name': f'Sim_spread1_amp01_flat80_anticlmulti_norotate_resample_',
                 'set_num': 0,
                 'energy': 62,
                 'div': 60,
                 'cent': 8,
                 'color': 'r',
                 'label': 's1a01'}

    return dist_dict


def dist2():
    dist_type = 'flat80_anticlmulti_spread1_amp1_resample'
    dist_dict = {'base_path': f'F:/Research/Data_Sim/{dist_type}/',
                 'set_name': f'Sim_spread1_amp1_flat80_anticlmulti_norotate_resample_',
                 'set_num': 0,
                 'energy': 62,
                 'div': 60,
                 'cent': 8,
                 'color': 'b',
                 'label': 's1a1'}

    return dist_dict


def show_dists(dists):
    ratio_bin_step = 0.02
    ratio_binning = np.arange(0 - ratio_bin_step / 2, 1 + ratio_bin_step * 3 / 2, ratio_bin_step)
    ratio_bin_centers = (ratio_binning[1:] + ratio_binning[:-1]) / 2

    fig, ax = plt.subplots()
    fig_norm, ax_norm = plt.subplots()
    fig_kde, ax_kde = plt.subplots()
    ax.axvline(0.5, color='black', ls='dotted', alpha=0.4)
    ax_norm.axvline(0.5, color='black', ls='dotted', alpha=0.4)
    # colors = dict(zip(divs, ['b', 'g', 'r', 'c', 'm', 'y', 'k'][:len(divs)]))

    for d in dists:
        path = f'{d["base_path"]}{d["set_name"]}{d["set_num"]}/{d["energy"]}GeV/' \
               f'ratios_divisions_{d["div"]}_centrality_{d["cent"]}_local.txt'

        az_data = AzimuthBinData(div=d['div'], path=path)
        print_stats(az_data)
        ratio_hist = az_data.get_ratio_dist(ratio_binning)
        ax.bar(ratio_bin_centers, ratio_hist, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5,
               label=f'{d["div"]}$^\circ$')
        ratio_hist_norm = az_data.get_ratio_dist(ratio_binning, True)
        ax_norm.bar(ratio_bin_centers, ratio_hist_norm, ratio_binning[1:] - ratio_binning[:-1], alpha=0.5,
                    label=f'{d["div"]}$^\circ$')
        ratio_unbin = az_data.get_ratio_dist()
        values = []
        for x, num in ratio_unbin.items():
            for i in range(num):
                values.append(x)
        sns.histplot(ax=ax_kde, x=values, stat='probability', bins=ratio_binning, kde=False, color=d["color"])
        # sns.rugplot(ax=ax_kde, x=values)

    fig_kde.legend(labels=[d['label'] for d in dists])
    fig.legend()
    ax.grid()
    fig_norm.legend()
    ax_norm.grid()
    plt.show()


if __name__ == '__main__':
    main()
