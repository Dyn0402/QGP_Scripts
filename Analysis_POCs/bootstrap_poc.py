#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 14 12:53 PM 2022
Created in PyCharm
Created as QGP_Scripts/bootstrap_poc.py

@author: Dylan Neff, Dylan
"""

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from scipy.stats import sem
import os
from multiprocessing import Pool
import tqdm
import istarmap

from pickle_methods import *
from poc_functions import *


def main():
    bootstrap_validation()
    print('donzo')


def bootstrap_validation():
    """
    Simulate binomials and test resampling bootstrap uncertainties against known answer
    :return:
    """
    seed = 13434
    threads = 15
    n_tracks = 15
    n_sample = 1
    n_events = 10000
    bin_width = np.deg2rad(120)
    bootstraps = 2000
    experiments = 100
    alg = 4
    # plot_out_dir = '/home/dylan/Research/Results/Resample_POC/nsample1440_nevent10000/'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Bootstrap_Validation/'
    plot_out_name = 'nsample1_nevent10k_bw120_ntrack15_nexp100_nbs2k/'
    plot_out_dir = plot_out_base + plot_out_name
    plot_exps_out_dir = f'{plot_out_dir}Experiments/'
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    try:
        os.mkdir(plot_exps_out_dir)
    except FileExistsError:
        pass
    show_plot = False

    stats = define_stats(n_tracks, bin_width)

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']

    write_info_file(plot_out_dir, threads, n_tracks, n_sample, n_events, bin_width, bootstraps, experiments, stats_plt)

    seeds = np.random.SeedSequence(seed).spawn(experiments)
    jobs = [(n_tracks, n_events, bin_width, n_sample, bootstraps, stats, stats_plt, seed,
             n_exp, True, plot_exps_out_dir, alg) for n_exp, seed in enumerate(seeds)]

    exp_stats = []
    with Pool(threads) as pool:
        for exp_stat in tqdm.tqdm(pool.istarmap(run_experiment, jobs), total=len(jobs)):
            exp_stats.append(exp_stat)

    plt.clf()
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
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir, stats_err_list[stat])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir,
                         stat_errs_delta=stats_err_delta_list[stat])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir, stats_err_list[stat],
                         stats_err_delta_list[stat])
        plot_bs_vs_delta(n_exps, stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        plot_bs_vs_delta_hist(stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        plot_err_vs_val_diff(stats_list[stat], stats[stat]['true'], stat, stats_err_list[stat],
                             stats_err_delta_list[stat], plot_out_dir)
        plot_exp_deviation(stats_list[stat], stats[stat]['true'], stat, plot_out_dir)
        plot_exp_sigmas(stats_list[stat], stats_err_list[stat], stats[stat]['true'], stat, plot_out_dir)

    plot_exp_scatter_stats(stats_plt, n_exps, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list)
    plot_exp_sigmas_stats(stats_plt, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list)
    plot_err_vs_val_diff_stats(stats_plt, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list)

    if show_plot:
        plt.show()


def bootstrap_vs_nsamples():
    """
    Simulate binomials and test resampling bootstrap uncertainties against known answer
    :return:
    """
    seed = 13434
    threads = 15
    n_tracks = 15
    n_sample = 1
    n_events = 250
    bin_width = np.deg2rad(120)
    bootstraps = 250
    experiments = 100
    alg = 4
    # plot_out_dir = '/home/dylan/Research/Results/Resample_POC/nsample1440_nevent10000/'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Bootstrap_Validation/'
    plot_out_name = 'nsample1_nevent250_bw120_ntrack15_nexp100/'
    plot_out_dir = plot_out_base + plot_out_name
    plot_exps_out_dir = f'{plot_out_dir}Experiments/'
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    try:
        os.mkdir(plot_exps_out_dir)
    except FileExistsError:
        pass
    show_plot = True

    stats = define_stats(n_tracks, bin_width)

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']

    write_info_file(plot_out_dir, threads, n_tracks, n_sample, n_events, bin_width, bootstraps, experiments, stats_plt)

    seeds = np.random.SeedSequence(seed).spawn(experiments)
    jobs = [(n_tracks, n_events, bin_width, n_sample, bootstraps, stats, stats_plt, seed,
             n_exp, True, plot_exps_out_dir, alg) for n_exp, seed in enumerate(seeds)]

    exp_stats = []
    with Pool(threads) as pool:
        for exp_stat in tqdm.tqdm(pool.istarmap(run_experiment, jobs), total=len(jobs)):
            exp_stats.append(exp_stat)

    plt.clf()
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
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir, stats_err_list[stat])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir,
                         stat_errs_delta=stats_err_delta_list[stat])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, plot_out_dir, stats_err_list[stat],
                         stats_err_delta_list[stat])
        plot_bs_vs_delta(n_exps, stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        plot_bs_vs_delta_hist(stats_err_list[stat], stats_err_delta_list[stat], stat, plot_out_dir)
        plot_exp_deviation(stats_list[stat], stats[stat]['true'], stat, plot_out_dir)
        plot_exp_sigmas(stats_list[stat], stats_err_list[stat], stats[stat]['true'], stat, plot_out_dir)

    plot_exp_scatter_stats(stats_plt, n_exps, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list)
    plot_exp_sigmas_stats(stats_plt, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list)

    if show_plot:
        plt.show()


def plot_exp_scatter(n_exps, stat_vals, binom_val, stat, out_dir=None, stat_errs_bs=None, stat_errs_delta=None):
    if out_dir is not None:
        fig, ax = plt.subplots()
    else:
        ax = plt.gca()
    ax.grid()
    ax.axhline(binom_val, ls='--', color='black', label='Binomial True Value')
    title = f'Scatter_{stat}'
    if stat_errs_bs is not None:
        ax.errorbar(n_exps, stat_vals, yerr=stat_errs_bs, marker='', ls='none', label='Bootstrap Errors')
        title += '_bserr'
    if stat_errs_delta is not None:
        ax.errorbar(n_exps, stat_vals, yerr=stat_errs_delta, marker='', ls='none', color='green', elinewidth=3,
                    alpha=0.5, label='Delta Theorem Errors')
        title += '_dterr'
    ax.scatter(n_exps, stat_vals, marker='o')

    if out_dir is not None:
        ax.set_xlabel('Exerpiment #')
        ax.set_title(stat.capitalize())
        ax.legend()
        fig.savefig(f'{out_dir}{title}.png', bbox_inches='tight')


def plot_err_vs_val_diff(stat_vals, binom_val, stat, bs_errs=None, delta_errs=None, out_dir=None):
    if out_dir is not None:
        fig, ax = plt.subplots()
    else:
        ax = plt.gca()
    ax.set_title(stat.title())
    ax.grid()
    ax.axvline(0, ls='--', color='black')
    if bs_errs is not None:
        ax.scatter(stat_vals - binom_val, bs_errs, color='blue', alpha=0.8, label='Bootstrap')
    if delta_errs is not None:
        ax.scatter(stat_vals - binom_val, delta_errs, color='green', alpha=0.8, label='Delta Theorem')
    ax.set_xlabel('Experiment Value - Binomial Value')
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim() + binom_val)
    ax2.tick_params(direction='in')
    ax2.set_xlabel('Experiment Value')
    if out_dir is not None:
        ax.set_ylabel('Error Value')
        ax.legend()
        fig.savefig(f'{out_dir}Err_vs_Val_Scatter_{stat}.png', bbox_inches='tight')


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
    ax.axvline(mean, color='green', ls=':', label='mean')
    ax.axvspan(mean - mean_err, mean + mean_err, color='green', alpha=0.6)
    ax.axvline(0, ls='--', color='black')
    ax.set_xlabel('Deviation from True')
    ax.set_title(stat)
    ax.legend()
    fig.savefig(f'{out_dir}Deviation_{stat}.png', bbox_inches='tight')


def plot_exp_sigmas(stat_vals, stat_errs, binom_val, stat, out_dir=None):
    sigmas = []
    for i in range(len(stat_vals)):
        sigmas.append((stat_vals[i] - binom_val) / stat_errs[i])

    if out_dir is not None:
        fig, hist = plt.subplots()
    else:
        hist = plt.gca()
    bin_edges = hist.hist(sigmas)[1]
    bin_width = bin_edges[1] - bin_edges[0]
    mean = np.mean(sigmas)
    mean_err = sem(sigmas)
    hist.axvline(mean, color='green', ls=':', label='mean')
    hist.axvspan(mean - mean_err, mean + mean_err, color='green', alpha=0.6)
    hist.axvline(0, ls='--', color='black')
    x = np.linspace(min(sigmas), max(sigmas), 1000)
    y = norm.pdf(x) * len(sigmas) * bin_width
    hist.plot(x, y, color='red', alpha=0.7, label='Standard Normal')

    if out_dir is not None:
        hist.set_xlabel('Sigmas from True')
        hist.set_ylabel('Number of Experiments')
        hist.set_title(stat.capitalize())
        hist.legend()
        fig.savefig(f'{out_dir}Sigmas_{stat}.png', bbox_inches='tight')


def plot_exp_scatter_stats(stats_plt, n_exps, stats_list, stats, plot_out_dir, stats_err_bs_list, stats_err_dt_list):
    fig_dt, axs_dt = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(10, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_dt[stat_index])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, stat_errs_delta=stats_err_dt_list[stat])
        axs_dt[stat_index].set_ylabel(stat.capitalize())
    axs_dt[-1].legend()
    axs_dt[-1].set_xlabel('Experiment #')
    fig_dt.tight_layout()
    fig_dt.savefig(f'{plot_out_dir}Scatter_Stats_dterr.png', bbox_inches='tight')

    fig_bs, axs_bs = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(10, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_bs[stat_index])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, stat_errs_bs=stats_err_bs_list[stat])
        axs_bs[stat_index].set_ylabel(stat.capitalize())
    axs_bs[-1].legend()
    axs_bs[-1].set_xlabel('Experiment #')
    fig_bs.tight_layout()
    fig_bs.savefig(f'{plot_out_dir}Scatter_Stats_bserr.png', bbox_inches='tight')

    fig_bsdt, axs_bsdt = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(10, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_bsdt[stat_index])
        plot_exp_scatter(n_exps, stats_list[stat], stats[stat]['true'], stat, stat_errs_bs=stats_err_bs_list[stat],
                         stat_errs_delta=stats_err_dt_list[stat])
        axs_bsdt[stat_index].set_ylabel(stat.capitalize())
    axs_bsdt[-1].legend()
    axs_bsdt[-1].set_xlabel('Experiment #')
    fig_bsdt.tight_layout()
    fig_bsdt.savefig(f'{plot_out_dir}Scatter_Stats_bsdterr.png', bbox_inches='tight')


def plot_exp_sigmas_stats(stats_plt, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list):
    fig_dt, axs_dt = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_dt[stat_index])
        plot_exp_sigmas(stats_list[stat], stats_err_delta_list[stat], stats[stat]['true'], stat)
        axs_dt[stat_index].set_ylabel('Number of Experiments')
        plt.annotate(stat.capitalize(), xy=(0.02, 0.98), xycoords='axes fraction', verticalalignment='top')
        # bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
    axs_dt[0].legend()
    axs_dt[-1].set_xlabel('Sigmas from True')
    fig_dt.tight_layout()
    fig_dt.savefig(f'{plot_out_dir}Sigmas_Stats_dterr.png', bbox_inches='tight')

    fig_bs, axs_bs = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_bs[stat_index])
        plot_exp_sigmas(stats_list[stat], stats_err_list[stat], stats[stat]['true'], stat)
        axs_bs[stat_index].set_ylabel('Number of Experiments')
        plt.annotate(stat.capitalize(), xy=(0.02, 0.98), xycoords='axes fraction', verticalalignment='top')
    axs_bs[0].legend()
    axs_bs[-1].set_xlabel('Sigmas from True')
    fig_bs.tight_layout()
    fig_bs.savefig(f'{plot_out_dir}Sigmas_Stats_bserr.png', bbox_inches='tight')


def plot_err_vs_val_diff_stats(stats_plt, stats_list, stats, plot_out_dir, stats_err_list, stats_err_delta_list):
    fig_bsdt, axs_bsdt = plt.subplots(1, len(stats_plt), figsize=(12, 6))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_bsdt[stat_index])
        plot_err_vs_val_diff(stats_list[stat], stats[stat]['true'], stat, stats_err_list[stat],
                             stats_err_delta_list[stat])
        # plt.annotate(stat.capitalize(), xy=(0.02, 0.98), xycoords='axes fraction', verticalalignment='top')
        # bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
    axs_bsdt[-1].legend()
    axs_bsdt[0].set_ylabel('Error Value')
    fig_bsdt.tight_layout()
    fig_bsdt.savefig(f'{plot_out_dir}Err_vs_Val_Scatter_Stats_dtbs.png', bbox_inches='tight')

    fig_bs, axs_bs = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    for stat_index, stat in enumerate(stats_plt):
        plt.sca(axs_bs[stat_index])
        plot_err_vs_val_diff(stats_list[stat], stats[stat]['true'], stat, stats_err_list[stat])
        # plt.annotate(stat.capitalize(), xy=(0.02, 0.98), xycoords='axes fraction', verticalalignment='top')
    axs_bs[-1].legend()
    axs_bs[0].set_ylabel('Error Value')
    fig_bs.tight_layout()
    fig_bs.savefig(f'{plot_out_dir}Err_vs_Val_Scatter_Stats_bs.png', bbox_inches='tight')


if __name__ == '__main__':
    main()
