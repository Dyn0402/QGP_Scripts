#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 21 8:32 PM 2023
Created in PyCharm
Created as QGP_Scripts/bes_v2_plot

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from analyze_binom_slices import read_flow_values


def main():
    main_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
               'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}
    energies = [7, 11, 39]

    fig_v2, ax_v2 = plt.subplots(dpi=144)
    plt.xticks(rotation=30)
    ax_v2.grid()
    ax_v2.axhline(0, color='black')
    fig_res, ax_res = plt.subplots(dpi=144)
    plt.xticks(rotation=30)
    ax_res.grid()
    ax_res.axhline(0, color='black')

    for energy in energies:
        v2_file_path = f'{main_dir}{energy}GeV/v2.txt'
        cent_bins, v2_vals, v2_errs, res_vals, res_errs = [], [], [], [], []
        with open(v2_file_path, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]:
                line = line.strip().split('\t')
                cent_bins.append(cent_map[int(line[0])])
                v2_vals.append(float(line[1].split()[0]))
                v2_errs.append(float(line[1].split()[1]))
                res_vals.append(float(line[2].split()[0]))
                res_errs.append(float(line[2].split()[1]))
        cent_bins, v2_vals, v2_errs, res_vals, res_errs = [list(reversed(x)) for x in
                                                           [cent_bins, v2_vals, v2_errs, res_vals, res_errs]]
        ax_v2.errorbar(cent_bins, v2_vals, v2_errs, ls='none', marker='o', alpha=0.8, label=f'{energy}GeV')
        ax_res.errorbar(cent_bins, res_vals, res_errs, ls='none', marker='o', alpha=0.8, label=f'{energy}GeV')
    ax_v2.legend()
    ax_res.legend()
    ax_v2.set_xlabel('Centrality')
    ax_res.set_xlabel('Centrality')
    ax_v2.set_ylabel('v2')
    ax_res.set_ylabel('Event Plane Resolution')
    fig_v2.tight_layout()
    fig_res.tight_layout()

    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
