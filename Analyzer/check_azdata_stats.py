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
    base_path = '/home/dylan/Research/Data/'
    set_name = 'eta05_n1ratios_dca3'
    set_num = 0
    energy = 7
    div = 120
    cent = 8

    path = f'{base_path}{set_name}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'

    az_data = AzimuthBinData(div=div, path=path)

    # az_data.print_dist()
    ratio_dist = az_data.get_ratio_dist()

    az_ratio_stats = DistStats(dist=ratio_dist)
    print(f'mean: {az_ratio_stats.get_mean()}')
    print(f'sd: {az_ratio_stats.get_sd()}')
    print(f'skew: {az_ratio_stats.get_skewness()}')
    print(f'kurt: {az_ratio_stats.get_kurtosis()}')
    print(f'kurt*var: {az_ratio_stats.get_kurt_var()}')

    print('donzo')


if __name__ == '__main__':
    main()
