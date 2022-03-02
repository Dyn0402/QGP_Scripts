#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 5:28 PM 2020
Created in PyCharm
Created as QGP_Scripts/check_azdata_stats.py

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats


def main():
    # base_path = '/home/dylan/Research/Data/default/'
    base_path = 'D:/Research/Data_Ampt_Old/default_resample/'
    base_mix_path = 'D:/Research/Data_Ampt_Old_Mix/default_resample/'
    # set_name = 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_'
    set_name = 'Ampt_rapid05_resample_norotate_'
    set_num = 0
    energy = 62
    div = 120
    cent = 8

    total_protons = [4]

    path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_mix = f'{base_mix_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

    az_data = AzimuthBinData(div=div, path=path)
    az_data_mix = AzimuthBinData(div=div, path=path_mix)

    # print_stats(az_data)
    print_binom_stats(az_data, az_data_mix, total_protons)

    # az_data.plot_ratio_dist(show=False)
    # az_data.plot_pull_dist(show=True)

    print('donzo')


def print_stats(az_data):
    ratio_dist = az_data.get_ratio_dist()
    az_ratio_stats = DistStats(dist=ratio_dist)
    print(f'Path: {az_data.path}')
    print(f'Ratio:')
    print(f'mean: {az_ratio_stats.get_mean()}')
    print(f'sd: {az_ratio_stats.get_sd()}')
    print(f'skew: {az_ratio_stats.get_skewness()}')
    print(f'kurt: {az_ratio_stats.get_kurtosis()}')
    print(f'kurt*var: {az_ratio_stats.get_kurt_var()}\n')

    pull_dist = az_data.get_pull_dist()
    az_pull_stats = DistStats(dist=pull_dist)
    print(f'Pull:')
    print(f'mean: {az_pull_stats.get_mean()}')
    print(f'sd: {az_pull_stats.get_sd()}')
    print(f'skew: {az_pull_stats.get_skewness()}')
    print(f'kurt: {az_pull_stats.get_kurtosis()}')
    print(f'kurt*var: {az_pull_stats.get_kurt_var()}\n')


def print_binom_stats(az_data, az_data_mix, total_protons):
    dist = az_data.get_dist()
    dist_mix = az_data_mix.get_dist()

    for total_proton in total_protons:
        tp_stats = DistStats(dist=dist[total_proton])
        print(f'Binomial Slice Raw Total Protons {total_proton}:')
        print(f'mean: {tp_stats.get_mean()}')
        print(f'sd: {tp_stats.get_sd()}')
        print(f'skew: {tp_stats.get_skewness()}')
        print(f'kurt: {tp_stats.get_kurtosis()}')
        print(f'kurt*var: {tp_stats.get_kurt_var()}\n')

        tp_stats_mix = DistStats(dist=dist_mix[total_proton])
        print(f'Binomial Slice Mix Total Protons {total_proton}:')
        print(f'mean: {tp_stats_mix.get_mean()}')
        print(f'sd: {tp_stats_mix.get_sd()}')
        print(f'skew: {tp_stats_mix.get_skewness()}')
        print(f'kurt: {tp_stats_mix.get_kurtosis()}')
        print(f'kurt*var: {tp_stats_mix.get_kurt_var()}\n')

        print(f'Binomial Slice Divide Total Protons {total_proton}:')
        print(f'mean: {tp_stats.get_mean() / tp_stats_mix.get_mean()}')
        print(f'sd: {tp_stats.get_sd() / tp_stats_mix.get_sd()}')
        print(f'skew: {tp_stats.get_skewness() / tp_stats_mix.get_skewness()}')
        print(f'kurt: {tp_stats.get_kurtosis() / tp_stats_mix.get_kurtosis()}')
        print(f'kurt*var: {tp_stats.get_kurt_var() / tp_stats_mix.get_kurt_var()}\n')


if __name__ == '__main__':
    main()
