#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 28 12:27 PM 2021
Created in PyCharm
Created as QGP_Scripts/cbwc_sim.py

@author: Dylan Neff, dylan
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
from DistStats import DistStats


def main():
    ref3_vals = range(490, 700)
    ref3_events = np.asarray([int(1e7*np.exp(-0.02 * ref3)) for ref3 in ref3_vals])
    binom_n = np.asarray([20 for ref3 in ref3_vals])
    binom_p = 0.4

    plot_input_dists(ref3_vals, ref3_events, binom_n)
    hists = simulate(ref3_vals, ref3_events, binom_n, binom_p)
    plot_dists(hists)
    moments = calc_moments(hists)
    plot_moments(moments, dict(zip(ref3_vals, binom_n)), binom_p)

    print('donzo')


def plot_input_dists(ref3_vals, ref3_events, binom_n):
    fig1, ax1 = plt.subplots()
    ax1.plot(ref3_vals, ref3_events)

    fig2, ax2 = plt.subplots()
    ax2.plot(ref3_vals, binom_n)

    plt.show()


def simulate(ref3_vals, ref3_events, binom_n, binom_p):
    hists = {}
    for ref3, events, n in zip(ref3_vals, ref3_events, binom_n):
        ref3_binom = binom(n, binom_p)
        hist, bin_edges = np.histogram(ref3_binom.rvs(size=events), np.arange(-0.5, n + 1.5, 1))
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        hists.update({ref3: hist})

    return hists


def plot_dists(hists):
    n_plots_y = int(len(hists)**0.5 + 0.5)
    n_plots_x = int(len(hists) / n_plots_y + 0.5)
    fig, axs = plt.subplots(n_plots_y, n_plots_x, sharey=True, sharex=True)

    ref3s = list(hists.keys())
    i = 0
    for y in range(n_plots_y):
        for x in range(n_plots_x):
            if i < len(ref3s):
                norm = sum(hists[ref3s[i]].values())
                axs[y, x].bar(hists[ref3s[i]].keys(), np.asarray(list(hists[ref3s[i]].values())) / norm, 1.0,
                              align='center')
                axs[y, x].text(0, 0.25, norm)
            i += 1

    fig.canvas.manager.window.showMaximized()
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()


def calc_moments(hists):
    ref3s = []
    sds = []
    skews = []
    kurts = []
    for ref3, hist in hists.items():
        stats = DistStats(hist)
        kurts.append(stats.get_kurtosis().val)
        skews.append(stats.get_skewness().val)
        sds.append(stats.get_sd().val)
        ref3s.append(ref3)

    return ref3s, sds, skews, kurts


def plot_moments(moments, ns, p):
    ref3s, sds, skews, kurts = moments

    q = 1 - p
    binom_sd = lambda n: np.sqrt(n * p * q)
    binom_skew = lambda n: (q - p) / np.sqrt(n * p * q)
    binom_kurt = lambda n: (1 - 6 * p * q) / (n * p * q)

    fig1, ax1 = plt.subplots()
    ax1.scatter(ref3s, sds, label='Simulation')
    ax1.plot(ref3s, [binom_sd(ns[ref]) for ref in ref3s], 'r', label='Binomial Expectation')
    ax1.legend()
    ax1.set_ylabel('Standard Deviation')

    fig2, ax2 = plt.subplots()
    ax2.scatter(ref3s, skews, label='Simulation')
    ax2.plot(ref3s, [binom_skew(ns[ref]) for ref in ref3s], 'r', label='Binomial Expectation')
    ax2.legend()
    ax2.set_ylabel('Skewness')

    fig3, ax3 = plt.subplots()
    ax3.scatter(ref3s, kurts, label='Simulation')
    ax3.plot(ref3s, [binom_kurt(ns[ref]) for ref in ref3s], 'r', label='Binomial Expectation')
    ax3.legend()
    ax3.set_ylabel('Kurtosis')

    plt.show()


if __name__ == '__main__':
    main()
