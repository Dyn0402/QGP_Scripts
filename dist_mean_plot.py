#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 15 2:03 PM 2020
Created in PyCharm
Created as QGP_Scripts/dist_mean_plot.py

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt


def main():
    x = [7, 11, 19, 27, 39, 62]
    bes_mean = [22.57, 17.7, 14.33, 13.39, 11.07, 10.46]
    bes_mean_err = [1.154, 1.084, 1.029, 1.007, 0.9675, 0.9577]
    ampt_mean = [24.91, 23.05, 18.16, 15.75, 13.92, 12.06]
    ampt_mean_err = [1.238, 1.159, 1.088, 1.05, 1.016, 0.9855]

    plt.errorbar(x, bes_mean, yerr=bes_mean_err, fmt='bo', label='BES1')
    plt.errorbar(x, ampt_mean, yerr=ampt_mean_err, fmt='go', label='AMPT')
    plt.legend()
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Mean of Proton Distribution')
    plt.title('Proton Distribution Means: AMPT vs BES1')
    plt.show()

    ampt_bes = [x/bes_mean[i] for i, x in enumerate(ampt_mean)]
    ampt_bes_err = [abs(x)*((ampt_mean_err[i]/ampt_mean[i])**2 + (bes_mean_err[i]/bes_mean[i])**2)**0.5 for i, x
                    in enumerate(ampt_bes)]

    plt.errorbar(x, ampt_bes, yerr=ampt_bes_err, fmt='ro')
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Proton Mean AMPT/BES1')
    plt.title('Proton Distribution Means: AMPT / BES1')
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
