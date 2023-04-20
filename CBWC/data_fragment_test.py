#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 20 10:45 AM 2023
Created in PyCharm
Created as QGP_Scripts/data_fragment_test

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import poisson
from DistStats import DistStats


def main():
    data = poisson.rvs(mu=10, size=10000)
    data2 = poisson.rvs(mu=20, size=100)
    data = np.concatenate([data, data2])
    np.random.shuffle(data)
    print(data)
    stats = DistStats(data, unbinned=True)
    c3 = stats.get_cumulant(3)
    k3 = stats.get_k_stat(3)
    print(c3, k3)
    n_subsets = 100
    sub_c3s, sub_k3s = [], []
    for sub_data in np.array_split(data, n_subsets):
        sub_stats = DistStats(sub_data, unbinned=True)
        sub_c3s.append(sub_stats.get_cumulant(3))
        sub_k3s.append(sub_stats.get_k_stat(3))
        print(sub_c3s[-1], sub_k3s[-1])

    plt.figure()
    plt.title('C3')
    plt.axhline(c3.val, color='red')
    plt.axhspan(c3.val - c3.err, c3.val + c3.err, color='red', alpha=0.3)
    plt.errorbar(np.arange(n_subsets), [sub_c3.val for sub_c3 in sub_c3s], [sub_c3.err for sub_c3 in sub_c3s],
                 ls='none', marker='o')
    sub_mean = np.mean(sub_c3s)
    plt.axhline(sub_mean.val, color='blue')
    plt.axhspan(sub_mean.val - sub_mean.err, sub_mean.val + sub_mean.err, color='blue', alpha=0.3)

    plt.figure()
    plt.title('K3')
    plt.axhline(k3.val, color='red')
    plt.axhspan(k3.val - k3.err, k3.val + k3.err, color='red', alpha=0.3)
    plt.errorbar(np.arange(n_subsets), [sub_k3.val for sub_k3 in sub_k3s], [sub_k3.err for sub_k3 in sub_k3s],
                 ls='none', marker='o')
    sub_mean = np.mean(sub_k3s)
    plt.axhline(sub_mean.val, color='blue')
    plt.axhspan(sub_mean.val - sub_mean.err, sub_mean.val + sub_mean.err, color='blue', alpha=0.3)

    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
