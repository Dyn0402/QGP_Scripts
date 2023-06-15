#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 15 11:57 AM 2023
Created in PyCharm
Created as QGP_Scripts/sum_of_uniform_gaus_dists

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    # check_uniform()
    check_gaus()
    print('donzo')


def check_uniform():
    num_rands = 10000
    num_dists = 2
    low_bound, up_bound = -1, 1
    dist_range = up_bound - low_bound

    binning = np.linspace(low_bound + dist_range / 100, up_bound + dist_range / 100, 100)

    rvs = np.random.random((num_dists, num_rands)) * dist_range + low_bound
    # print(rvs)
    rvs_sum = rvs.sum(axis=0)
    # print(rvs_sum)
    print(f'Standard Deviation of Distribution: {np.std(rvs_sum)}')
    print(f'Square of Sum of Variances: {np.sqrt(num_dists * dist_range**2 / 12)}')

    fig, ax = plt.subplots()
    print(rvs.shape)
    ax.hist(rvs[0], bins=binning)

    fig2, ax2 = plt.subplots()
    ax2.hist(rvs_sum, bins=binning*2)

    plt.show()


def check_gaus():
    num_rands = 10000
    num_dists = 5
    mean, sigma = 0, 1
    # dist_range = up_bound - low_bound

    binning = np.linspace(mean - 5 * sigma, mean + 5 * sigma, 100)

    rvs = np.random.normal(mean, sigma, size=(num_dists, num_rands))
    # print(rvs)
    rvs_sum = rvs.sum(axis=0)
    # print(rvs_sum)
    print(f'Standard Deviation of Distribution: {np.std(rvs_sum)}')
    print(f'Square of Sum of Variances: {np.sqrt(num_dists * sigma**2)}')

    fig, ax = plt.subplots()
    print(rvs.shape)
    ax.hist(rvs[0], bins=binning)

    fig2, ax2 = plt.subplots()
    ax2.hist(rvs_sum, bins=binning*2)

    plt.show()


if __name__ == '__main__':
    main()
