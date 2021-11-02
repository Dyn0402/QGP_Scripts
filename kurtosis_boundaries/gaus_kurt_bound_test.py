#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 01 9:57 PM 2021
Created in PyCharm
Created as QGP_Scripts/gaus_kurt_bound_test

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import operator

from DistStats import DistStats


def main():
    p = 0.6
    q = 1 - p
    n = 50
    dist = binom(n, p)
    vals = dist.rvs(size=1000000)
    bins = np.arange(-0.5, 50.5)

    counts = np.histogram(vals, bins)[0]
    dist_dict = dict(zip((bins[:-1] + bins[1:]) / 2, counts))
    stats = DistStats(dist_dict)
    print(stats.get_kurtosis())
    print((1 - 6 * p * q) / (n * p * q))

    fig, ax = plt.subplots()
    ax.hist(vals, bins)
    # plt.show()

    bounds_r = np.arange(30, 50)
    bounds_l = np.arange(0, 30)
    for bounds, op in zip([bounds_r, bounds_l], [operator.lt, operator.gt]):
        kurts = []
        kurt_errs = []
        for bound in bounds:
            vals_bound = [val if op(val, bound) else bound for val in vals]
            counts_bnd = np.histogram(vals_bound, bins)[0]
            dist_dict_bnd = dict(zip((bins[:-1] + bins[1:]) / 2, counts_bnd))
            stats_bnd = DistStats(dist_dict_bnd)
            print(f'{bound}: {stats_bnd.get_kurtosis()}')
            kurts.append(stats_bnd.get_kurtosis().val)
            kurt_errs.append(stats_bnd.get_kurtosis().err)

        kurts = np.asarray(kurts)
        kurt_errs = np.asarray(kurt_errs)
        fig, ax = plt.subplots()
        ax.fill_between(bounds, kurts + kurt_errs, kurts - kurt_errs)
        ax.axhline((1 - 6 * p * q) / (n * p * q), color='red', ls='--')
        ax.set_ylabel('kurtosis')
        ax.set_xlabel('Upper x bound')

    # plt.hist(vals_bound, bins)
    plt.show()


    print('donzo')


if __name__ == '__main__':
    main()
