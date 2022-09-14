#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 14 12:55 PM 2022
Created in PyCharm
Created as QGP_Scripts/poc_functions.py

@author: Dylan Neff, Dylan
"""

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binom

from Measure import Measure
from sub_event_resample_algorithm import get_resamples
from pickle_methods import *


def gen_event(n_tracks):
    return np.random.random(n_tracks) * 2 * np.pi


def gen_experiment(n_events, n_tracks, rng=np.random.default_rng()):
    return np.sort(rng.random((n_events, n_tracks))) * 2 * np.pi


def run_experiment(n_tracks, n_events, bin_width, samples, bootstraps, stats,
                   stats_plt, seed, n_exp=None, plot=False, out_dir=''):
    rng = np.random.default_rng(seed)
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
            plot_bootstraps(bs_list, stat, stat_vals[stat], stat_errs[stat], stats[stat]['true'], n_exp, n_events,
                            out_dir)

    return n_exp, stat_vals, stat_errs, stat_errs_delta


def bin_experiment(experiment, n_tracks, bin_width, samples, bootstraps, rng):
    data = np.zeros(n_tracks + 1, dtype=int)
    data_bs = np.zeros((bootstraps, n_tracks + 1), dtype=int)
    for event in experiment:
        # event = rotate_event(event, rng.random() * 2 * np.pi)  # This doesn't matter with no phi dependence
        hist = get_resamples(event, bin_width, samples)
        data += hist
        for bootstrap in data_bs:
            for x in range(rng.poisson(1)):
                bootstrap += hist

    return data, data_bs


def run_experiment_no_bs(n_tracks, n_events, bin_width, samples, stats,
                         stats_plt, seed, n_exp=None):
    rng = np.random.default_rng(seed)
    experiment = gen_experiment(n_events, n_tracks, rng)
    data = bin_experiment_no_bs(experiment, n_tracks, bin_width, samples, rng)

    data_stats = DistStats(data)
    stat_vals = {}
    stat_errs_delta = {}
    for stat in stats_plt:
        meas = stats[stat]['meth'](data_stats)
        stat_vals.update({stat: meas.val})
        stat_errs_delta.update({stat: meas.err})

    return n_exp, samples, n_events, bin_width, n_tracks, stat_vals, stat_errs_delta


def bin_experiment_no_bs(experiment, n_tracks, bin_width, samples, rng):
    data = np.zeros(n_tracks + 1, dtype=int)
    for event in experiment:
        # event = rotate_event(event, rng.random() * 2 * np.pi)  # This doesn't matter with no phi dependence
        hist = get_resamples(event, bin_width, samples)
        data += hist

    return data


def rotate_event(event, rotate_angle):
    """
    Rotate all tracks (float values [0,2pi)) in event by rotate_angle. Ensure output range remains [0, 2pi)
    :param event:
    :param rotate_angle:
    :return:
    """
    event = event + rotate_angle
    while np.any(event >= 2 * np.pi):
        event = np.where(event >= 2 * np.pi, event - 2 * np.pi, event)

    return np.sort(event)


def get_bs_sem(vals, num_bs=250):
    bs_means = []
    for bs_index in range(num_bs):
        bs_vals = np.random.choice(vals, size=vals.size)
        bs_means.append(np.mean(bs_vals))

    return np.std(np.array(bs_means))


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


def write_info_file(plot_out_dir, threads, n_tracks, n_sample, n_events, bin_widths, bootstraps, experiments,
                    stats_plt):
    with open(plot_out_dir + 'info.txt', 'w') as file:
        file.write(f'threads: {threads}\n')
        file.write(f'n_tracks: {n_tracks}\n')
        file.write(f'n_sample: {n_sample}\n')
        file.write(f'n_events: {n_events}\n')
        file.write(f'bin_widths: {bin_widths}\n')
        file.write(f'bootstraps: {bootstraps}\n')
        file.write(f'experiments: {experiments}\n')
        file.write('stats: ')
        for stat in stats_plt:
            file.write(f'{stat}, ')


def plot_dist(data, n_tracks, bin_width, n_exp, out_dir=None):
    x = range(len(data))
    sns.displot(x=x, weights=data, discrete=True, kde=False, label='Simulation')
    plt.scatter(x, np.sum(data) * binom.pmf(x, n_tracks, bin_width / (2 * np.pi)), color='red', label='Binomial',
                marker='_', zorder=4)
    plt.title(f'Experiment #{n_exp} Distribution')
    plt.xlabel('Particles in Bin')
    plt.legend()
    if out_dir is not None:
        plt.savefig(f'{out_dir}Experiment_{n_exp}_dist.png', bbox_inches='tight')
        plt.close()


def plot_bootstraps(bs_list, stat, exp_val, exp_err, binom_val, n_exp, n_events, out_dir):
    sns.displot(bs_list, kde=True, rug=True, label='Bootstrap Estimates')
    plt.title(f'Experiment #{n_exp}')
    plt.axvline(exp_val, color='red', ls='--', label='Experiment Estimate')
    plt.axvspan(exp_val - exp_err, exp_val + exp_err, color='red', alpha=0.3)
    plt.axvline(binom_val, color='green', ls='--', label='Binomial True')
    binom_val_dec_match = str(Measure(binom_val, exp_err)).split(' ')[0]
    plt.annotate(f'Estimate {Measure(exp_val, exp_err)}\nBinom True {binom_val_dec_match}\n{n_events} Events',
                 xy=(0.02, 0.98), xycoords='axes fraction', bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
                 verticalalignment='top')
    plt.xlabel(stat)
    plt.legend(loc='upper right')
    plt.savefig(f'{out_dir}Experiment_{n_exp}_{stat}.png', bbox_inches='tight')
    plt.close()
