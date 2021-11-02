#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 01 1:15 PM 2021
Created in PyCharm
Created as QGP_Scripts/anti_cluster_multi

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def main():
    means = [0.7, 3.5, 2.4]
    sds = [0.2, 0.3, 0.5]
    amps = 1 / (np.sqrt(2 * np.pi) * np.asarray(sds))
    aclust_amp = 0.4
    dists = [norm(mean, sd) for mean, sd in zip(means, sds)]
    x = np.linspace(0, 2 * np.pi, 1000)
    base = 1
    base_int = base * (x[-1] - x[0])
    p_clust = np.ones(1000)
    p_aclust = np.ones(1000)

    fig_ind_clust, ax_ind_clust = plt.subplots()
    fig_ind_aclust, ax_ind_aclust = plt.subplots()
    for dist, amp, sd in zip(dists, amps, sds):
        # clust_pdf = base + clust_amp(aclust_amp, sd, base_int) * dist.pdf(x) / amp
        clust_pdf = base + aclust_amp * dist.pdf(x) / amp
        # clust_pdf /= np.sum(clust_pdf) * (x[-1] - x[0]) / len(x)
        aclust_pdf = base - aclust_amp * dist.pdf(x) / amp
        # aclust_pdf /= np.sum(aclust_pdf) * (x[-1] - x[0]) / len(x)
        ax_ind_clust.plot(x, clust_pdf)
        ax_ind_aclust.plot(x, aclust_pdf)
        ax_ind_clust.set_ylim(bottom=0)
        ax_ind_aclust.set_ylim(bottom=0)
        p_clust *= clust_pdf
        p_aclust *= aclust_pdf

    fig_clust, ax_clust = plt.subplots()
    fig_aclust, ax_aclust = plt.subplots()
    ax_clust.plot(x, p_clust)
    ax_aclust.plot(x, p_aclust)
    ax_clust.set_ylim(bottom=0)
    ax_aclust.set_ylim(bottom=0)

    plt.show()
    print('donzo')


def clust_amp(aclust_amp, sigma, base_int):
    return 1 / (1 / aclust_amp - 2 * sigma * np.sqrt(2 * np.pi) / base_int)


if __name__ == '__main__':
    main()
