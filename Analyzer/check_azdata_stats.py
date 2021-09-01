#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 5:28 PM 2020
Created in PyCharm
Created as QGP_Scripts/check_azdata_stats.py

@author: Dylan Neff, dylan
"""

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats


def main():
    base_path = '/home/dylan/Research/Data_Ampt/default/'
    set_name = 'Ampt_rapid05_n1ratios_'
    set_num = 0
    energy = 7
    div = 120
    cent = 3

    path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

    az_data = AzimuthBinData(div=div, path=path)

    # az_data.print_dist()
    pull_dist = az_data.get_pull_dist()

    az_pull_stats = DistStats(dist=pull_dist)
    print(path)
    print(f'mean: {az_pull_stats.get_mean()}')
    print(f'sd: {az_pull_stats.get_sd()}')
    print(f'skew: {az_pull_stats.get_skewness()}')
    print(f'kurt: {az_pull_stats.get_kurtosis()}')
    print(f'kurt*var: {az_pull_stats.get_kurt_var()}')

    print('donzo')


if __name__ == '__main__':
    main()
