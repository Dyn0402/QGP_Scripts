#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 08 6:44 PM 2021
Created in PyCharm
Created as QGP_Scripts/sub_event_poc

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.stats import norm
from multiprocessing import Pool

from DistStats import DistStats
from sub_event_resample_algorithm import get_resamples


def main():
    """
    Simulate binomials and test resampling/bootstrap against known answer
    :return:
    """

    # get_resamples(np.array([0.5, 1.5]), 2.09, 3)
    # return
    n_tracks = 15
    # n_samples = np.arange(3, 100, 2)
    n_sample = 3
    n_events = 100
    bin_width = np.deg2rad(120)
    bootstraps = 250
    experiments = 10

    n = n_tracks
    p = bin_width / (2 * np.pi)
    q = 1 - p
    binom_ans = {'mean': n * p, 'standard deviation': np.sqrt(n * p * q), 'skewness': (q - p) / np.sqrt(n * p * q),
                 'non-excess kurtosis': (1 - 6 * p * q) / (n * p * q) + 3}

    pois = poisson(1)

    stats = {'mean': getattr(DistStats, 'get_mean'), 'standard deviation': getattr(DistStats, 'get_sd'),
             'skewness': getattr(DistStats, 'get_skewness'), 'kurtosis': getattr(DistStats, 'get_kurtosis'),
             'non-excess kurtosis': getattr(DistStats, 'get_non_excess_kurtosis')}

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']

    stats_list = {stat: [] for stat in stats_plt}
    stats_err_list = {stat: [] for stat in stats_plt}
    for experiment in range(experiments):
        print(f'Experiment {experiment}')
        # data = {bin_count: 0 for bin_count in range(n_tracks + 1)}
        # data_bs = [{bin_count: 0 for bin_count in range(n_tracks + 1)} for i in range(bootstraps)]
        data = np.zeros(n_tracks + 1, dtype=int)
        data_bs = np.zeros((bootstraps, n_tracks + 1), dtype=int)
        for i in range(n_events):
            event = np.sort(gen_event(n_tracks))
            hist = get_resamples(event, bin_width, n_sample)
            data += hist
            for bootstrap in data_bs:
                for x in range(pois.rvs()):
                    bootstrap += hist

        data_stats = DistStats(data)
        data_bs_stats = [DistStats(bs) for bs in data_bs]
        for stat in stats_plt:
            stats_list[stat].append(stats[stat](data_stats).val)
            bs_list = []
            for i in range(len(data_bs_stats)):
                bs_list.append(stats[stat](data_bs_stats[i]).val)
            stats_err_list[stat].append(np.std(bs_list))

    for stat in stats_plt:
        fig_scat, ax_scat = plt.subplots()
        ax_scat.errorbar(range(len(stats_list[stat])), stats_list[stat], yerr=stats_err_list[stat],
                         marker='o', ls='none')
        ax_scat.axhline(binom_ans[stat], ls='--', color='black')
        ax_scat.set_xlabel('Exerpiment #')
        ax_scat.set_title(stat)
        ax_scat.grid()

        sigmas = []
        for i in range(len(stats_list[stat])):
            sigmas.append((stats_list[stat][i] - binom_ans[stat]) / stats_err_list[stat][i])

        fig_hist, ax_hist = plt.subplots()
        ax_hist.hist(sigmas, density=True)
        x = np.linspace(min(sigmas), max(sigmas), 1000)
        y = norm.pdf(x)
        ax_hist.plot(x, y, color='red', alpha=0.7)
        ax_hist.set_xlabel('Sigmas from True')
        ax_hist.set_title(stat)

    plt.show()

    plot_vs_samples()
    print('donzo')


def plot_vs_samples():
    pass


def gen_event(n_tracks):
    return np.random.random(n_tracks) * 2 * np.pi


def gen_experiment(n_tracks, n_events):
    return np.random.random((n_events, n_tracks))


def run_experiment(n_tracks, n_events, samples):
    pass


if __name__ == '__main__':
    main()
