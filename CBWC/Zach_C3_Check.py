#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 20 9:43 AM 2023
Created in PyCharm
Created as QGP_Scripts/Zach_C3_Check

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

import uproot
import awkward as ak
import vector


def main():
    base_path = 'C:/Users/Dylan/Desktop/Zach_C3_Check/'
    with uproot.open(f'{base_path}3p5.root') as file:
        ref3_hist = file['hFXTMult3']
        # print(ref3_hist)
        # print(dir(ref3_hist))
        # print(ref3_hist.to_numpy())
        # print(ref3_hist.values())
        ref3_numpy = ref3_hist.to_numpy()
        bin_cents = (ref3_numpy[1][1:] + ref3_numpy[1][:-1]) / 2
        bin_widths = (ref3_numpy[1][1:] - ref3_numpy[1][:-1])
        # print(len(bin_cents))
        # print(len(bin_widths))
        # print(len(ref3_numpy[0]))
        plt.figure()
        plt.bar(bin_cents, ref3_numpy[0], width=bin_widths)
        ref3_counts = {ref3: count for ref3, count in zip(bin_cents, ref3_numpy[0])}
    with uproot.open(f'{base_path}out_3p5_n0p5_0_norm.root') as file:
        c1_hist = file['C1'].tojson()
        c3_hist = file['C3'].tojson()
        plt.figure()
        plt.errorbar(c3_hist['fX'], c3_hist['fY'], c3_hist['fEY'], c3_hist['fEX'], ls='none', label='C3')
        plt.plot(c1_hist['fX'], np.array(c1_hist['fY']), color='red', label='Poisson')
        k3 = [c3 * ref3_counts[ref3]**2 / ((ref3_counts[ref3] - 1) * (ref3_counts[ref3] - 2))
              for ref3, c3 in zip(c3_hist['fX'], c3_hist['fY'])]
        plt.errorbar(c3_hist['fX'], k3, c3_hist['fEY'], c3_hist['fEX'], ls='none', label='k3')
        plt.xlabel('FXTMult3')
        plt.ylim(-50, 100)
        plt.legend()
        plt.tight_layout()
    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
