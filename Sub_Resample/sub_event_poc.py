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
from scipy.stats import sem
from multiprocessing import Pool

from DistStats import DistStats
from Measure import Measure
from sub_event_resample_algorithm import get_resamples
from pickle_methods import *


def main():
    """
    Simulate binomials and test resampling/bootstrap against known answer
    :return:
    """
    threads = 16
    n_tracks = 15
    n_sample = 1
    n_events = 100
    bin_width = np.deg2rad(120)
    bootstraps = 250
    experiments = 1000
    plot_out_dir = '/home/dylan/Research/Results/Resample_POC/'

    stats = define_stats(n_tracks, bin_width)

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']

    with Pool(threads) as pool:
        exp_stats = pool.starmap(run_experiment,
                                 [(n_tracks, n_events, bin_width, n_sample, bootstraps, stats, stats_plt,
                                   n_exp, True, plot_out_dir) for n_exp in range(experiments)])

    n_exps = []
    stats_list = {stat: [] for stat in stats_plt}
    stats_err_list = {stat: [] for stat in stats_plt}
    stats_err_delta_list = {stat: [] for stat in stats_plt}

    for n_exp, stat_vals, stat_errs, stat_errs_delta in exp_stats:
        n_exps.append(n_exp)
        for stat in stats_plt:
            stats_list[stat].append(stat_vals[stat])
            stats_err_list[stat].append(stat_errs[stat])
            stats_err_delta_list[stat].append(stat_errs_delta[stat])

    for stat in stats_plt:
        plot_exp_scatter(n_exps, stats_list[stat], stats_err_list[stat], stats_err_delta_list[stat],
                         stats[stat]['true'], stat, plot_out_dir)
        plot_bs_vs_delta(n_exps, stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        plot_bs_vs_delta_hist(stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        # plot_exp_deviation(stats_list[stat], stats[stat]['true'], stat, plot_out_dir)
        # plot_exp_sigmas(stats_list[stat], stats_err_list[stat], stats[stat]['true'], stat, plot_out_dir)

    plt.show()

    plot_vs_samples()
    print('donzo')


def plot_vs_samples():
    pass


def gen_event(n_tracks):
    return np.random.random(n_tracks) * 2 * np.pi


def gen_experiment(n_events, n_tracks, rng=np.random.default_rng()):
    return np.sort(rng.random((n_events, n_tracks))) * 2 * np.pi


def run_experiment(n_tracks, n_events, bin_width, samples, bootstraps, stats,
                   stats_plt, n_exp=None, plot=False, out_dir=''):
    print(f'Experiment {n_exp}')
    rng = np.random.default_rng()
    experiment = gen_experiment(n_events, n_tracks, rng)
    data, data_bs = bin_experiment(experiment, n_tracks, bin_width, samples, bootstraps, rng)

    if plot:
        plot_dist(data, n_tracks, bin_width, n_exp, out_dir)

    data_stats = DistStats(data)
    data_bs_stats = [DistStats(bs) for bs in data_bs]
    stat_vals = {}
    stat_errs = {}
    stat_errs_delta = {}
    for stat in stats_plt:
        meas = stats[stat]['meth'](data_stats)
        stat_vals.update({stat: meas.val})
        stat_errs_delta.update({stat: meas.err})
        bs_list = [stats[stat]['meth'](bs_dist).val for bs_dist in data_bs_stats]
        stat_errs.update({stat: np.std(bs_list)})
        if plot:
            plot_bootstraps(bs_list, stat, stat_vals[stat], stat_errs[stat], stats[stat]['true'], n_exp, out_dir)

    return n_exp, stat_vals, stat_errs, stat_errs_delta


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


def plot_dist(data, n_tracks, bin_width, n_exp, out_dir):
    x = range(len(data))
    sns.displot(x=x, weights=data, discrete=True, kde=False, label='Simulation')
    plt.scatter(x, np.sum(data) * binom.pmf(x, n_tracks, bin_width / (2 * np.pi)), color='red', label='Binomial',
                marker='_', zorder=4)
    plt.title(f'Experiment #{n_exp} Distribution')
    plt.xlabel('Particles in Bin')
    plt.legend()
    plt.savefig(f'{out_dir}Experiment_{n_exp}_dist.png', bbox_inches='tight')
    plt.close()


def plot_bootstraps(bs_list, stat, exp_val, exp_err, binom_val, n_exp, out_dir):
    sns.displot(bs_list, kde=True, rug=True, label='Bootstrap Estimates')
    plt.title(f'Experiment #{n_exp}')
    plt.axvline(exp_val, color='red', ls='--', label='Experiment Estimate')
    plt.axvline(binom_val, color='green', ls='--', label='Binomial True')
    plt.annotate(f'{Measure(exp_val, exp_err)}', xy=(0.02, 0.98), xycoords='axes fraction',
                 bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3), verticalalignment='top')
    plt.xlabel(stat)
    plt.legend(loc='upper right')
    plt.savefig(f'{out_dir}Experiment_{n_exp}_{stat}.png', bbox_inches='tight')
    plt.close()


def plot_exp_scatter(n_exps, stat_vals, stat_errs, stat_errs_delta, binom_val, stat, out_dir):
    fig, ax = plt.subplots()
    ax.grid()
    ax.errorbar(n_exps, stat_vals, yerr=stat_errs, marker='o', ls='none', label='Bootstrap Errors')
    ax.errorbar(n_exps, stat_vals, yerr=stat_errs_delta, marker='', ls='none', color='green', elinewidth=3,
                alpha=0.5, label='Delta Theorem Errors')
    ax.axhline(binom_val, ls='--', color='black')
    ax.set_xlabel('Exerpiment #')
    ax.set_title(stat)
    ax.legend()
    fig.savefig(f'{out_dir}Scatter_{stat}.png', bbox_inches='tight')


def plot_bs_vs_delta(n_exps, bs_errs, delta_errs, stat, out_dir):
    fig, ax = plt.subplots()
    ax.set_title(stat)
    ax.grid()
    ax.scatter(n_exps, np.array(bs_errs) - np.array(delta_errs))
    ax.axhline(0, color='black')
    ax.set_ylabel('Bootstrap Error - Delta Error')
    ax.set_xlabel('Experiment #')
    fig.savefig(f'{out_dir}Bs_minus_Delta_Scatter_{stat}.png', bbox_inches='tight')


def plot_bs_vs_delta_hist(bs_errs, delta_errs, stat, out_dir):
    # fig = plt.figure()
    sns.displot(x=(np.array(bs_errs) - np.array(delta_errs)), kde=True, rug=True)
    plt.title(stat)
    plt.axvline(0, color='black')
    plt.xlabel('Bootstrap Error - Delta Error')
    plt.savefig(f'{out_dir}Bs_minus_Delta_Hist_{stat}.png', bbox_inches='tight')


def plot_exp_deviation(stat_vals, binom_val, stat, out_dir):
    deviations = []
    for i in range(len(stat_vals)):
        deviations.append(stat_vals[i] - binom_val)

    fig, ax = plt.subplots()
    ax.hist(deviations)
    mean = np.mean(deviations)
    mean_err = sem(deviations)
    ax.axvspan(mean - mean_err, mean + mean_err, color='green', alpha=0.8, label='mean')
    ax.axvline(0, ls='-', color='gray', alpha=1.0)
    ax.set_xlabel('Deviation from True')
    ax.set_title(stat)
    ax.legend()
    fig.savefig(f'{out_dir}Deviation_{stat}.png', bbox_inches='tight')


def plot_exp_sigmas(stat_vals, stat_errs, binom_val, stat, out_dir):
    sigmas = []
    for i in range(len(stat_vals)):
        sigmas.append((stat_vals[i] - binom_val) / stat_errs[i])

    fig, hist = plt.subplots()
    hist.hist(sigmas, density=True)
    mean = np.mean(sigmas)
    mean_err = sem(sigmas)
    hist.axvspan(mean - mean_err, mean + mean_err, color='green', alpha=0.8, label='mean')
    hist.axvline(0, ls='-', color='gray', alpha=1.0)
    x = np.linspace(min(sigmas), max(sigmas), 1000)
    y = norm.pdf(x)
    hist.plot(x, y, color='red', alpha=0.7)
    hist.set_xlabel('Sigmas from True')
    hist.set_title(stat)
    hist.legend()
    fig.savefig(f'{out_dir}Sigmas_{stat}.png', bbox_inches='tight')


def define_stats(n_tracks, bin_width):
    n = n_tracks
    p = bin_width / (2 * np.pi)
    q = 1 - p
    stats = {'mean': {'meth': get_mean_meas, 'true': n * p},
             'standard deviation': {'meth': get_sd_meas, 'true': np.sqrt(n * p * q)},
             'skewness': {'meth': get_skewness_meas, 'true': (q - p) / np.sqrt(n * p * q)},
             'kurtosis': {'meth': get_kurtosis_meas, 'true': (1 - 6 * p * q) / (n * p * q)},
             'non-excess kurtosis': {'meth': get_nekurtosis_meas, 'true': (1 - 6 * p * q) / (n * p * q) + 3},
             'c1': {'meth': get_c1_meas, 'true': n * p},
             'c2': {'meth': get_c2_meas, 'true': p * q * n},
             'c3': {'meth': get_c3_meas, 'true': p * q * (1 - 2 * p) * n},
             'c4': {'meth': get_c4_meas, 'true': p * q * (6 * p ** 2 - 6 * p + 1) * n},
             'c5': {'meth': get_c5_meas, 'true': p * q * (1 - 2 * p) * (12 * p ** 2 - 12 * p + 1) * n},
             'c6': {'meth': get_c6_meas, 'true': p * q * (120 * p ** 4 - 240 * p ** 3 + 150 * p ** 2 - 30 * p + 1) * n},
             'k1': {'meth': get_k1_meas, 'true': n * p},
             'k2': {'meth': get_k2_meas, 'true': p * q * n},
             'k3': {'meth': get_k3_meas, 'true': p * q * (1 - 2 * p) * n},
             'k4': {'meth': get_k4_meas, 'true': p * q * (6 * p ** 2 - 6 * p + 1) * n},
             'k5': {'meth': get_k5_meas, 'true': p * q * (1 - 2 * p) * (12 * p ** 2 - 12 * p + 1) * n},
             'k6': {'meth': get_k6_meas, 'true': p * q * (120 * p ** 4 - 240 * p ** 3 + 150 * p ** 2 - 30 * p + 1) * n},
             'c4/c2': {'meth': get_c4_div_c2_meas, 'true': 6 * p ** 2 - 6 * p + 1},
             'k4/k2': {'meth': get_k4_div_k2_meas, 'true': 6 * p ** 2 - 6 * p + 1},
             'c6/c2': {'meth': get_c6_div_c2_meas, 'true': 120 * p ** 4 - 240 * p ** 3 + 150 * p ** 2 - 30 * p + 1},
             'k6/k2': {'meth': get_k6_div_k2_meas, 'true': 120 * p ** 4 - 240 * p ** 3 + 150 * p ** 2 - 30 * p + 1},
             'c4/c2 - k4/k2': {'meth': get_c4_div_c2_sub_k4_div_k2_meas, 'true': 0},
             'c6/c2 - k6/k2': {'meth': get_c6_div_c2_sub_k6_div_k2_meas, 'true': 0},
             }

    return stats


if __name__ == '__main__':
    main()
