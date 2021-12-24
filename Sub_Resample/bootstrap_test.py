#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 04 5:17 PM 2021
Created in PyCharm
Created as QGP_Scripts/bootstrap_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Bootstrap_Az_Bin import BootstrapAzBin as Bab
from DistStats import DistStats as Ds


def main():
    num_protons = 10
    div = 60
    path = '/home/dylan/Research/Data_Sim/single10_anticlmulti_spread0_amp0_test/' \
           'Sim_spread0_amp0_single10_anticlmulti_n1ratios_0/62GeV/ratios_divisions_60_centrality_8_local.txt'
    # path = '/home/dylan/Research/Data_Sim/single10_anticlmulti_resample_spread0_amp0_test/' \
    #        'Sim_spread0_amp0_single10_anticlmulti_resample_n1ratios_0/62GeV/ratios_divisions_60_centrality_8_local.txt'
    # path = 'D:/Research/Data_Sim/single10_anticlmulti_resample_spread0_amp0_test/' \
    #        'Sim_spread0_amp0_single10_anticlmulti_resample_n1ratios_0/62GeV/ratios_divisions_60_centrality_8_local.txt'
    bs_dists = Bab(div, path)
    bs_dists.data.print_dist()
    rs_stats = Ds(bs_dists.data.data[num_protons])
    print(rs_stats.get_sd())
    print('\nBootstraps:')
    sds = []
    n = num_protons
    p = div / 360.0
    q = 1 - p
    binom_sd = np.sqrt(n * p * q)
    print(binom_sd)
    print(f'Num bootstraps: {len(bs_dists.data_bs)}')
    for bs in bs_dists.data_bs:
        bs_stats = Ds(bs.data[num_protons])
        sds.append(bs_stats.get_sd().val)
    sns.displot(sds, kde=True)
    plt.axvline(rs_stats.get_sd().val, color='red', ls='--')
    plt.axvline(binom_sd, color='green', ls='--')
    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
