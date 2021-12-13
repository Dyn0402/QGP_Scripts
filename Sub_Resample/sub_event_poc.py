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

from DistStats import DistStats
from sub_event_resample_algorithm import get_resamples


def main():
    """
    Simulate binomials and test resampling/bootstrap against known answer
    :return:
    """
    n_tracks = 15
    n_samples = np.arange(3, 100, 2)
    n_events = 100
    bin_width = np.deg2rad(120)
    bootstraps = 100
    experiments = 10

    n = n_tracks
    p = bin_width / (2 * np.pi)
    q = 1 - p
    binom_ans = {'mean': n * p, 'standard deviation': np.sqrt(n * p * q),
                 'non-excess kurtosis': (1 - 6 * p * q) / (n * p * q) + 3}

    pois = poisson(1)

    stats = {'mean': getattr(DistStats, 'get_mean'), 'standard deviation': getattr(DistStats, 'get_sd'),
             'skewness': getattr(DistStats, 'get_skewness'), 'kurtosis': getattr(DistStats, 'get_kurtosis'),
             'non-excess kurtosis': getattr(DistStats, 'get_non_excess_kurtosis')}

    stats_plt = ['non-excess kurtosis']

    ensemble = []

    stats_list = {stat: [] for stat in stats_plt}
    stats_err_list = {stat: [] for stat in stats_plt}
    for n_sample in n_samples:
        print(f'{n_sample} samples')
        data = {bin_count: 0 for bin_count in range(n_tracks + 1)}
        data_bs = [{bin_count: 0 for bin_count in range(n_tracks + 1)} for i in range(bootstraps)]
        for experiment in range(experiments):
            for i in range(n_events):
                event = list(sorted(gen_event(n_tracks)))
                hist = get_resamples(event, bin_width, n_sample)
                for count in hist:
                    data[count] += 1
                for bootstrap in data_bs:
                    for x in range(pois.rvs()):
                        for count in hist:
                            bootstrap[count] += 1

            data_stats = DistStats(data)
            data_bs_stats = [DistStats(bs) for bs in data_bs]
            for stat in stats_plt:
                stats_list[stat].append(stats[stat](data_stats).val)
                bs_list = []
                for i in range(len(data_bs_stats)):
                    bs_list.append(stats[stat](data_bs_stats[i]).val)
                stats_err_list[stat].append(np.std(bs_list))

    plt.errorbar(n_samples, stats_list[stats_plt[0]], yerr=stats_err_list[stats_plt[0]], marker='o',
                 ls='none')
    plt.axhline(binom_ans[stats_plt[0]])
    plt.show()

    plot_vs_samples()
    print('donzo')


def plot_vs_samples():
    pass


def gen_event(n_tracks):
    return np.random.random(n_tracks) * 2 * np.pi


if __name__ == '__main__':
    main()
