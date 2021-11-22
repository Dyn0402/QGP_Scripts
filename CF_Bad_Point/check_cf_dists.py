#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 21 11:26 PM 2021
Created in PyCharm
Created as QGP_Scripts/check_cf_dists

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from AzimuthBinData import AzimuthBinData
from DistStats import DistStats


def main():
    # For CF default 60 degrees, mixed medians are 0 systematics are nan. Need to find source of issue
    base_path = '/home/dylan/Research/Data_CF_Mix/default/'
    set_name = 'CF_rapid05_n1ratios_'
    set_nums = np.arange(0, 60)
    energy = 7
    div = 120
    cent = 8

    for set_num in set_nums:
        path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

        az_data = AzimuthBinData(div=div, path=path)

        ratio_dist = az_data.get_ratio_dist()
        az_ratio_stats = DistStats(dist=ratio_dist)
        # print(f'Ratio: {path}')
        # print(f'mean: {az_ratio_stats.get_mean()}')
        # print(f'sd: {az_ratio_stats.get_sd()}')
        # print(f'skew: {az_ratio_stats.get_skewness()}')
        # print(f'kurt: {az_ratio_stats.get_kurtosis()}')
        # print(f'kurt*var: {az_ratio_stats.get_kurt_var()}')
        # az_data.plot_ratio_dist(binning=20)

        # pull_dist = az_data.get_pull_dist()
        # az_pull_stats = DistStats(dist=pull_dist)
        # print(f'Pull: {path}')
        # print(f'mean: {az_pull_stats.get_mean()}')
        # print(f'sd: {az_pull_stats.get_sd()}')
        # print(f'skew: {az_pull_stats.get_skewness()}')
        # print(f'kurt: {az_pull_stats.get_kurtosis()}')
        # print(f'kurt*var: {az_pull_stats.get_kurt_var()}')
        # az_data.plot_pull_dist(binning=20)

        print(f'sd: {az_ratio_stats.get_sd()}')

        if (az_ratio_stats.get_sd() <= 0.01):
            print(f'Ratio: {path}')
            print(f'mean: {az_ratio_stats.get_mean()}')
            print(f'sd: {az_ratio_stats.get_sd()}')
            print(f'skew: {az_ratio_stats.get_skewness()}')
            print(f'kurt: {az_ratio_stats.get_kurtosis()}')
            print(f'kurt*var: {az_ratio_stats.get_kurt_var()}')
            az_data.plot_ratio_dist(binning=20)

    print('donzo')


if __name__ == '__main__':
    main()
