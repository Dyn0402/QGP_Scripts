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
import seaborn as sns
from scipy.stats import poisson
from scipy.stats import norm
from scipy.stats import binom
from multiprocessing import Pool

from DistStats import DistStats
from Measure import Measure
from sub_event_resample_algorithm import get_resamples


def main():
    """
    Simulate binomials and test resampling/bootstrap against known answer
    :return:
    """
    threads = 1
    n_tracks = 15
    # n_samples = np.arange(3, 100, 2)
    n_sample = 1440
    n_events = 100
    bin_width = np.deg2rad(120)
    bootstraps = 250
    experiments = 1
    plot_out_dir = '/home/dylan/Research/Results/Resample_POC/'

    n = n_tracks
    p = bin_width / (2 * np.pi)
    q = 1 - p
    binom_ans = {'mean': n * p, 'standard deviation': np.sqrt(n * p * q), 'skewness': (q - p) / np.sqrt(n * p * q),
                 'non-excess kurtosis': (1 - 6 * p * q) / (n * p * q) + 3}

    # pois = poisson(1)

    stat_methods = {'mean': getattr(DistStats, 'get_mean'), 'standard deviation': getattr(DistStats, 'get_sd'),
                    'skewness': getattr(DistStats, 'get_skewness'), 'kurtosis': getattr(DistStats, 'get_kurtosis'),
                    'non-excess kurtosis': getattr(DistStats, 'get_non_excess_kurtosis'),
                    'k2': getattr(DistStats, 'get_k_stat(2)')}

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis', 'k2']

    with Pool(threads) as pool:
        exp_stats = pool.starmap(run_experiment,
                                 [(n_tracks, n_events, bin_width, n_sample, bootstraps, stat_methods, stats_plt,
                                   binom_ans, n_exp, True, plot_out_dir) for n_exp in range(experiments)])

    n_exps = []
    stats_list = {stat: [] for stat in stats_plt}
    stats_err_list = {stat: [] for stat in stats_plt}

    for n_exp, stat_vals, stat_errs in sorted(exp_stats):
        print(n_exp)
        n_exps.append(n_exp)
        for stat in stats_plt:
            stats_list[stat].append(stat_vals[stat])
            stats_err_list[stat].append(stat_errs[stat])

    # for n_exp in range(experiments):
    #     print(f'Experiment {n_exp}')
    #     nexp, stat_vals, stat_errs = run_experiment(n_tracks, n_events, bin_width, n_sample, bootstraps, pois, stats,
    #                                           stats_plt, binom_ans, n_exp, True, plot_out_dir)
    #     for stat in stats_plt:
    #         stats_list[stat].append(stat_vals[stat])
    #         stats_err_list[stat].append(stat_errs[stat])

    for stat in stats_plt:
        fig_scat, ax_scat = plt.subplots()
        ax_scat.errorbar(range(len(stats_list[stat])), stats_list[stat], yerr=stats_err_list[stat],
                         marker='o', ls='none')
        ax_scat.axhline(binom_ans[stat], ls='--', color='black')
        ax_scat.set_xlabel('Exerpiment #')
        ax_scat.set_title(stat)
        ax_scat.grid()

        # sigmas = []
        # for i in range(len(stats_list[stat])):
        #     sigmas.append((stats_list[stat][i] - binom_ans[stat]) / stats_err_list[stat][i])

        # fig_hist, ax_hist = plt.subplots()
        # ax_hist.hist(sigmas, density=True)
        # x = np.linspace(min(sigmas), max(sigmas), 1000)
        # y = norm.pdf(x)
        # ax_hist.plot(x, y, color='red', alpha=0.7)
        # ax_hist.set_xlabel('Sigmas from True')
        # ax_hist.set_title(stat)

    plt.show()

    plot_vs_samples()
    print('donzo')


def plot_vs_samples():
    pass


def gen_event(n_tracks):
    return np.random.random(n_tracks) * 2 * np.pi


def gen_experiment(n_events, n_tracks, rng=np.random.default_rng()):
    return np.sort(rng.random((n_events, n_tracks))) * 2 * np.pi


def run_experiment(n_tracks, n_events, bin_width, samples, bootstraps, stat_methods,
                   stats_plt, binom_ans=None, n_exp=None, plot=False, out_dir=''):
    print(f'Experiment {n_exp}')
    rng = np.random.default_rng()
    experiment = gen_experiment(n_events, n_tracks, rng)
    data, data_bs = bin_experiment(experiment, n_tracks, bin_width, samples, bootstraps, rng)

    if plot:
        x = range(len(data))
        sns.displot(x=x, weights=data, discrete=True, kde=False, label='Simulation')
        plt.scatter(x, np.sum(data) * binom.pmf(x, n_tracks, bin_width / (2 * np.pi)), color='red', label='Binomial',
                    marker='_', zorder=4)
        plt.title(f'Experiment #{n_exp} Distribution')
        plt.xlabel('Particles in Bin')
        plt.legend()
        plt.savefig(f'{out_dir}Experiment_{n_exp}_dist.png', bbox_inches='tight')
        plt.close()

    data_stats = DistStats(data)
    data_bs_stats = [DistStats(bs) for bs in data_bs]
    stat_vals = {}
    stat_errs = {}
    for stat in stats_plt:
        stat_vals.update({stat: stat_methods[stat](data_stats).val})
        bs_list = [stat_methods[stat](bs_dist).val for bs_dist in data_bs_stats]
        stat_errs.update({stat: np.std(bs_list)})
        if plot:
            sns.displot(bs_list, kde=True, rug=True, label='Bootstrap Estimates')
            plt.title(f'Experiment #{n_exp}')
            plt.axvline(stat_vals[stat], color='red', ls='--', label='Experiment Estimate')
            plt.axvline(binom_ans[stat], color='green', ls='--', label='Binomial True')
            plt.annotate(f'{Measure(stat_vals[stat], stat_errs[stat])}', xy=(0.02, 0.98), xycoords='axes fraction',
                         bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3), verticalalignment='top')
            plt.xlabel(stat)
            plt.legend(loc='upper right')
            plt.savefig(f'{out_dir}Experiment_{n_exp}_{stat}.png', bbox_inches='tight')
            plt.close()

    return n_exp, stat_vals, stat_errs


def bin_experiment(experiment, n_tracks, bin_width, samples, bootstraps, rng):
    data = np.zeros(n_tracks + 1, dtype=int)
    data_bs = np.zeros((bootstraps, n_tracks + 1), dtype=int)
    for event in experiment:
        hist = get_resamples(event, bin_width, samples)
        data += hist
        for bootstrap in data_bs:
            for x in range(rng.poisson(1)):
                bootstrap += hist

    return data, data_bs


if __name__ == '__main__':
    main()
