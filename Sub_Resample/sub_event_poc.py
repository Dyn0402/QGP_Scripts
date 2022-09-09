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
from matplotlib.cm import get_cmap
import pandas as pd
import seaborn as sns
from scipy.stats import poisson
from scipy.stats import norm
from scipy.stats import binom
from scipy.stats import sem
import os
from itertools import product
from multiprocessing import Pool
import tqdm

from DistStats import DistStats
from Measure import Measure
from sub_event_resample_algorithm import get_resamples
from pickle_methods import *
import istarmap


def main():
    resample_validation()
    # resample_with_nsamples()
    # correlated_dists()
    # multiple_exps_dist()
    # bootstrap_validation()
    print('donzo')


def resample_validation():
    """
    Simulate binomials and test resampling values against known answer
    :return:
    """
    seed = 1432
    threads = 15
    n_tracks = [15]
    # n_samples = [1, 3, 1440]
    n_samples = [1, 2, 3]
    # n_samples = np.array([1, 2, 3, 4, 5, 6, 7] + list(np.arange(10, 5000, 50)))
    # n_samples = np.arange(1, 20, 1)
    n_events = np.arange(100, 2000, 5)
    # n_events = [250]
    # bin_widths = np.deg2rad([60, 120, 240, 300])
    bin_widths = np.deg2rad([120])
    experiments = 1000
    # plot_out_dir = '/home/dylan/Research/Results/Resample_POC/nsample1440_nevent10000/'
    # plot_out_base = 'F:/Research/Resample_POC/Resample_Validation/'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Resample_Validation/'
    # plot_out_base = 'C:/Users/Dyn04/Desktop/Resample_POC/Resample_Validation/'
    # plot_out_name = 'vs_nsamplesto5k_bws_ntrack15_nevent250/'
    # plot_out_name = 'vs_nsamplesto20_isobws_ntrack15_nevent250/'
    plot_out_name = 'vs_neventsto2k_bw120_ntrack15_nsamples123/'
    plot_out_dir = plot_out_base + plot_out_name
    plot_sds = True
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    show_plot = True

    stats = {bin_width: {n_track: define_stats(n_track, bin_width) for n_track in n_tracks}
             for bin_width in bin_widths}

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']
    # stats_plt = ['skewness', 'non-excess kurtosis']
    # stats_plt = ['standard deviation']

    write_info_file(plot_out_dir, threads, n_tracks, n_samples, n_events, bin_widths, 'resample_validiation',
                    experiments, stats_plt)

    num_indep_exps = experiments * len(n_events) * len(n_samples) * len(n_tracks) * len(bin_widths)
    seeds = iter(np.random.SeedSequence(seed).spawn(num_indep_exps))
    plot_data = []

    jobs = [(n_track, n_event, bin_width, n_sample, stats[bin_width][n_track], stats_plt, next(seeds), n_exp)
            for n_exp in range(experiments) for n_event in n_events for n_sample in n_samples
            for bin_width in bin_widths for n_track in n_tracks]

    with Pool(threads) as pool:
        for exp_stat in tqdm.tqdm(pool.istarmap(run_experiment_no_bs, jobs), total=len(jobs)):
            n_exp, n_samples_exp, n_events_exp, bin_width, n_track, stat_vals, stat_errs_delta = exp_stat
            for stat, val in stat_vals.items():
                plot_data.append({'stat': stat, 'n_exp': n_exp, 'n_samples': n_samples_exp, 'n_events': n_events_exp,
                                  'n_tracks': n_track, 'bin_width': bin_width, 'val': val,
                                  'delta_err': stat_errs_delta[stat],
                                  })

    plot_data = pd.DataFrame(plot_data)

    plot_vs_indep_var(plot_data, stats_plt, stats, 'n_events', plot_out_dir, plot_sds)
    # plot_vs_indep_var(plot_data, stats_plt, stats, 'n_samples', plot_out_dir, plot_sds)
    # plot_vs_indep_var(plot_data, stats_plt, stats, 'bin_width', plot_out_dir, plot_sds)

    if show_plot:
        plt.show()


def plot_vs_indep_var(plot_data, stats_plt, stats, indep_var, plot_out_dir, plot_sd=True):
    for stat in stats_plt:
        color = iter(get_cmap('Set1').colors)
        # color_binom = iter(get_cmap('tab20b').colors)
        stat_df = plot_data[plot_data['stat'] == stat]
        fig, ax = plt.subplots()
        fig_del, ax_del = plt.subplots()
        fig_del_norm, ax_del_norm = plt.subplots()
        ax.grid()
        ax_del.grid()
        ax_del_norm.grid()
        ax_del.axhline(0, color='black')
        ax_del_norm.axhline(0, color='black')
        ax_del_norm.axhline(1, color='black')

        indep_vals = pd.unique(plot_data[indep_var])  # Assume here a square lattice

        var_string_consts = {
            'n_tracks': {'x-label': 'Number of Tracks', 'legend': ' tracks'},
            'bin_width': {'x-label': 'Azimuthal Partition Width', 'legend': '° width'},
            'n_events': {'x-label': 'Number of Events', 'legend': ' events'},
            'n_samples': {'x-label': 'Number of Samples', 'legend': ' samples'}}

        set_vars = list(var_string_consts.keys())

        set_vars.remove(indep_var)
        uniques = {var: pd.unique(stat_df[var]) for var in set_vars}
        set_combos = product(*uniques.values())

        binom_val_label = False
        if indep_var == 'bin_width':
            for n_track in uniques['n_tracks']:
                binom_vals = []
                for bin_width in indep_vals:
                    binom_vals.append(stats[bin_width][n_track][stat]['true'])
                if not binom_val_label:
                    ax.plot(indep_vals, binom_vals, ls='--', color='black', label='Binomial Value')
                    binom_val_label = True
                else:
                    ax.plot(indep_vals, binom_vals, ls='--', color='black')
        elif indep_var == 'n_tracks':
            for bin_width in uniques['bin_width']:
                binom_vals = []
                for n_track in indep_vals:
                    binom_vals.append(stats[bin_width][n_track][stat]['true'])
                    if not binom_val_label:
                        ax.plot(indep_vals, binom_vals, ls='--', color='black', label='Binomial Value')
                        binom_val_label = True
                    else:
                        ax.plot(indep_vals, binom_vals, ls='--', color='black')
        else:
            for bin_width in uniques['bin_width']:
                for n_track in uniques['n_tracks']:
                    # label = f'{n_track} track, {round(np.rad2deg(bin_width))}° Binomial Value'
                    # c = next(color_binom)
                    if not binom_val_label:
                        ax.axhline(stats[bin_width][n_track][stat]['true'], ls='--', color='black',
                                   label='Binomial Value')
                        binom_val_label = True
                    else:
                        ax.axhline(stats[bin_width][n_track][stat]['true'], ls='--', color='black')

        for set_var_vals in set_combos:
            set_var_vals = dict(zip(set_vars, set_var_vals))
            set_df = stat_df
            for var, var_val in set_var_vals.items():
                set_df = set_df[set_df[var] == var_val]

            c = next(color)
            means, sds, sems, deltas, delta_sems, delta_sds = [], [], [], [], [], []
            for indep_val in indep_vals:
                vals = set_df[set_df[indep_var] == indep_val]['val']
                means.append(np.mean(vals))
                sds.append(np.std(vals))
                sems.append(sds[-1] / np.sqrt(vals.size))
                # delts = np.power(vals - stats[stat]['true'], 2)
                if indep_var == 'bin_width':
                    binom_val = stats[indep_val][set_var_vals['n_tracks']][stat]['true']
                elif indep_var == 'n_tracks':
                    binom_val = stats[set_var_vals['bin_width']][indep_val][stat]['true']
                else:
                    binom_val = stats[set_var_vals['bin_width']][set_var_vals['n_tracks']][stat]['true']
                delts = np.abs(vals - binom_val)
                # fig_deltas, ax_deltas = plt.subplots()
                # ax_deltas.set_title(f'{indep_var} = {indep_val} {set_var_vals}')
                # fig_deltas.canvas.manager.set_window_title(f'{indep_var} = {indep_val} {set_var_vals}')
                # sns.histplot(delts)
                # delta_sems.append(get_bs_sem(delts))  # Looks to be about equivalent to s/root(n)
                deltas.append(np.mean(delts))
                delta_sds.append(np.std(delts))
                delta_sems.append(np.std(delts) / np.sqrt(vals.size))
                # print(f'{indep_val}: {delta_sds[-1]}')
            means, sds, sems, deltas, delta_sems = (np.array(x) for x in (means, sds, sems, deltas, delta_sems))

            label = []
            for var in set_vars:
                if len(uniques[var]) > 1:
                    val = set_var_vals[var] if var != 'bin_width' else round(np.rad2deg(set_var_vals[var]))
                    label.append(f'{val}{var_string_consts[var]["legend"]}')
            label = ', '.join(label)
            ax.plot(indep_vals, means, label=label, color=c)
            ax.fill_between(indep_vals, means - sems, means + sems, color=c, alpha=0.6)

            ax_del.plot(indep_vals, deltas, label=label, color=c)
            ax_del.fill_between(indep_vals, deltas - delta_sems, deltas + delta_sems, color=c, alpha=0.5)

            del_max = max(deltas)
            ax_del_norm.plot(indep_vals, deltas / del_max, label=label, color=c)
            ax_del_norm.fill_between(indep_vals, (deltas - delta_sems) / del_max, (deltas + delta_sems) / del_max,
                                     color=c, alpha=0.5)

            if plot_sd:
                ax.fill_between(indep_vals, means - sds, means + sds, color=c, alpha=0.2)
                ax_del.fill_between(indep_vals, deltas - delta_sds, deltas + delta_sds, color=c, alpha=0.2)
                ax_del_norm.fill_between(indep_vals, (deltas - delta_sds) / del_max, (deltas + delta_sds) / del_max,
                                         color=c, alpha=0.2)
        ax.set_xlabel(var_string_consts[indep_var]['x-label'])
        # ax.set_ylabel('')
        ax.legend()
        ax_del.set_xlabel(var_string_consts[indep_var]['x-label'])
        ax_del.legend()
        ax_del_norm.set_xlabel(var_string_consts[indep_var]['x-label'])
        ax_del_norm.legend()

        title = f'{stat.title()} vs {var_string_consts[indep_var]["x-label"]}'
        title_del = f'{stat.title()} Deviations vs {var_string_consts[indep_var]["x-label"]}'
        title_del_norm = f'{stat.title()} Normalized Deviations vs {var_string_consts[indep_var]["x-label"]}'
        title_extra = [f'{vals[0]}{var_string_consts[var]["legend"]}' if var != 'bin_width'
                       else f'{round(np.rad2deg(vals[0]))}{var_string_consts[var]["legend"]}' for var, vals
                       in uniques.items() if len(vals) == 1]
        title_extra = ', '.join(title_extra)
        if title_extra != '':
            title += f'\n{title_extra}'
            title_del += f'\n{title_extra}'
            title_del_norm += f'\n{title_extra}'

        ax.set_title(title)
        ax_del.set_title(title_del)
        ax_del_norm.set_title(title_del_norm)

        fig_names = {
            fig: f'{stat}_vs_{indep_var.replace("_", "")}',
            fig_del: f'{stat}_absdev_vs_{indep_var.replace("_", "")}',
            fig_del_norm: f'{stat}_absdev_norm_vs_{indep_var.replace("_", "")}',
        }

        for fig_obj, fig_name in fig_names.items():
            fig_obj.tight_layout()
            fig_obj.canvas.manager.set_window_title(fig_name)
            fig_obj.savefig(f'{plot_out_dir}{fig_name}.png', bbox_inches='tight')


def bootstrap_validation():
    """
    Simulate binomials and test resampling bootstrap uncertainties against known answer
    :return:
    """
    seed = 13434
    threads = 15
    n_tracks = 15
    n_sample = 1440
    n_events = 400
    bin_width = np.deg2rad(120)
    bootstraps = 250
    experiments = 100
    # plot_out_dir = '/home/dylan/Research/Results/Resample_POC/nsample1440_nevent10000/'
    plot_out_base = 'F:/Research/Resample_POC/Bootstrap_Validation/'
    plot_out_name = 'nsample1440_nevent100_new/'
    plot_out_dir = plot_out_base + plot_out_name
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    show_plot = False

    stats = define_stats(n_tracks, bin_width)

    stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']

    write_info_file(plot_out_dir, threads, n_tracks, n_sample, n_events, bin_width, bootstraps, experiments, stats_plt)

    seeds = np.random.SeedSequence(seed).spawn(experiments)
    jobs = [(n_tracks, n_events, bin_width, n_sample, bootstraps, stats, stats_plt, seed,
             n_exp, True, plot_out_dir) for n_exp, seed in enumerate(seeds)]

    exp_stats = []
    with Pool(threads) as pool:
        for exp_stat in tqdm.tqdm(pool.istarmap(run_experiment, jobs), total=len(jobs)):
            exp_stats.append(exp_stat)

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
        plot_exp_deviation(stats_list[stat], stats[stat]['true'], stat, plot_out_dir)
        plot_exp_sigmas(stats_list[stat], stats_err_list[stat], stats[stat]['true'], stat, plot_out_dir)

    if show_plot:
        plt.show()


def correlated_dists():
    threads = 14
    seed_0 = 1432
    n_tracks = 15
    n_samples = 1
    n_events_sim = np.arange(100, 1e5 + 1, 1000, dtype=int)
    n_events_dist_plot = np.array([1e2, 1e3, 1e5], dtype=int)
    bin_width = np.deg2rad(120)
    stat_plt = 'standard deviation'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Visualizations/'
    plot_out_name = 'test3/'
    plot_out_dir = plot_out_base + plot_out_name
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    show_plot = True

    stats = define_stats(n_tracks, bin_width)

    n_events_sim = np.array(sorted([*set(list(n_events_sim) + list(n_events_dist_plot))]))

    # Same string of number for all experiments! I want it like this here to show filling of one experiment but be aware
    seeds_repeat = [seed_0 for n_event in n_events_sim]
    # Independent and uncorrlated seeds for each experiment
    seeds_indep = iter(np.random.SeedSequence(seed_0).spawn(len(n_events_sim)))

    for seeds, exp_name in zip([seeds_repeat, seeds_indep], ['Single_Exp', 'Independent_Exps']):
        print(f'Starting {exp_name}')
        fig, axs = plt.subplots(len(n_events_dist_plot), 1, sharex=True)
        event_axes = dict(zip(n_events_dist_plot, axs))
        plt.subplots_adjust(hspace=0)

        fig_stat, ax_stat = plt.subplots()
        ax_stat.set_xlabel('Number of Events')
        ax_stat.set_ylabel(stat_plt)

        # fig_dev, ax_dev = plt.subplots()
        # ax_dev.set_xlabel('Number of Events')
        # ax_dev.set_ylabel(f'{stat_plt} deviation from binomial')

        fig_dev_abs, ax_dev_abs = plt.subplots()
        ax_dev_abs.set_xlabel('Number of Events')
        ax_dev_abs.set_ylabel(r'$|\sigma_{data} - \sigma_{binomial}|$')

        jobs = []
        for n_event, seed in zip(n_events_sim, seeds):
            jobs.append((n_tracks, n_event, bin_width, n_samples, stats, [stat_plt], seed))

            if n_event in n_events_dist_plot:
                rng = np.random.default_rng(seed)
                experiment = gen_experiment(n_event, n_tracks, rng)
                hist = bin_experiment_no_bs(experiment, n_tracks, bin_width, n_samples)
                ax = event_axes[n_event]
                x = range(len(hist))
                sns.histplot(x=x, weights=hist, discrete=True, kde=False, label='Simulation', ax=ax)
                scatter = ax.scatter(x, np.sum(hist) * binom.pmf(x, n_tracks, bin_width / (2 * np.pi)), color='red',
                                     marker='_', zorder=4)
                if ax == axs[0]:
                    scatter.set_label('Binomial')
                    ax.legend()
                ax.text(0.7, 0.1, f'{n_event} Events', fontsize=12, transform=ax.transAxes)
                if ax == axs[-1]:
                    ax.set_xlabel('Particles in Bin')

        stat_vals, stat_errs, stat_meases = [], [], []
        with Pool(threads) as pool:
            for n_exp, samples, n_events, bin_width, n_track, stats_vals, stats_errs_delta in \
                    tqdm.tqdm(pool.istarmap(run_experiment_no_bs, jobs), total=len(jobs)):
                stat_meases.append(Measure(stats_vals[stat_plt], stats_errs_delta[stat_plt]))
                stat_vals.append(stats_vals[stat_plt])
                stat_errs.append(stats_errs_delta[stat_plt])

        devs = np.power(np.array(stat_vals) - stats[stat_plt]['true'], 2)
        devs_err = np.power(np.array(stat_meases) - stats[stat_plt]['true'], 2)
        devs_err = [x.err for x in devs_err]
        devs_abs = np.abs(np.array(stat_vals) - stats[stat_plt]['true'])

        # ax_dev.axhline(0, color='black')
        ax_dev_abs.axhline(0, color='black')

        # ax_dev.errorbar(n_events_sim, devs, yerr=devs_err, marker='o', ls='none')
        ax_dev_abs.errorbar(n_events_sim, devs_abs, yerr=stat_errs, marker='o', ls='none')

        ax_stat.axhline(stats[stat_plt]['true'], ls='--', color='red', label='Binomial')
        ax_stat.errorbar(n_events_sim, stat_vals, yerr=stat_errs, ls='none', marker='o', label='Simulation')
        ax_stat.legend()

        fig_names = {
            fig: f'dists_vs_binom_with_nevents_{exp_name}',
            fig_stat: f'sd_vs_events_{exp_name}',
            # fig_dev: f'sd_dev2_vs_events_{exp_name}',
            fig_dev_abs: f'sd_devabs_vs_events_{exp_name}',
        }

        for fig_obj, fig_name in fig_names.items():
            fig_obj.canvas.manager.set_window_title(fig_name)
            fig_obj.tight_layout()
            if plot_out_dir is not None:
                fig_obj.savefig(f'{plot_out_dir}{fig_name}.png')

    if show_plot:
        plt.show()


def multiple_exps_dist():
    threads = 15
    seed_0 = 1432
    n_tracks = 15
    n_samples = 1
    n_events = np.array([100, 1e3, 1e5], dtype=int)
    n_experiments = 300
    bin_width = np.deg2rad(120)
    stat_plt = 'standard deviation'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Visualizations/'
    plot_out_name = 'test4/'
    plot_out_dir = plot_out_base + plot_out_name
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    show_plot = True

    stats = define_stats(n_tracks, bin_width)

    # Independent and uncorrlated seeds for each experiment
    seeds = iter(np.random.SeedSequence(seed_0).spawn(n_experiments * len(n_events)))

    # print(f'Starting {exp_name}')
    fig, axs = plt.subplots(len(n_events), 1, sharex=True)
    event_axes = dict(zip(n_events, axs))
    plt.subplots_adjust(hspace=0)

    fig_abs, axs_abs = plt.subplots(len(n_events), 1, sharex=True)
    event_axes_abs = dict(zip(n_events, axs_abs))
    plt.subplots_adjust(hspace=0)

    jobs = [(n_tracks, n_event, bin_width, n_samples, stats, [stat_plt], next(seeds), n_exp)
            for n_exp in range(n_experiments) for n_event in n_events]

    stat_vals = {n_event: [] for n_event in n_events}
    with Pool(threads) as pool:
        for n_exp, n_sample, n_event, bin_width, n_track, stats_vals, stats_errs_delta in \
                tqdm.tqdm(pool.istarmap(run_experiment_no_bs, jobs), total=len(jobs)):
            stat_vals[n_event].append(stats_vals[stat_plt])

    for n_event in n_events:
        ax = event_axes[n_event]
        sns.histplot(x=stat_vals[n_event], kde=False, label='Simulation', ax=ax)
        mean = np.mean(stat_vals[n_event])
        sd = np.std(stat_vals[n_event])
        semean = sd / np.sqrt(n_experiments)
        ax.axvline(stats[stat_plt]['true'], ls='--', color='red', label='Binomial')
        ax.axvline(mean, color='black', label='mean', ls=':')
        ax.axvspan(mean - sd, mean + sd, color='gray', alpha=0.1, label='Standard Deviation')
        ax.axvspan(mean - semean, mean + semean, color='gray', alpha=0.6, label='Error on the Mean')
        if ax == axs[-1]:
            ax.legend(loc='upper left')
            ax.set_xlabel('Standard Deviation of Distribution')
        ax.text(0.75, 0.35, f'{n_event} Events', fontsize=12, transform=ax.transAxes)

        ax = event_axes_abs[n_event]
        abs_vals = abs(np.array(stat_vals[n_event]) - stats[stat_plt]['true'])
        sns.histplot(x=abs_vals, kde=False, label='Simulation', ax=ax)
        mean = np.mean(abs_vals)
        sd = np.std(abs_vals)
        semean = sd / np.sqrt(n_experiments)
        ax.axvline(0, color='black')
        ax.axvline(mean, color='black', label='mean', ls=':')
        ax.axvspan(mean - sd, mean + sd, color='gray', alpha=0.1, label='Standard Deviation')
        ax.axvspan(mean - semean, mean + semean, color='gray', alpha=0.6, label='Error on the Mean')
        if ax == axs_abs[-1]:
            ax.legend(loc='upper right')
            ax.set_xlabel(r'$|\sigma - \sigma_{binom}|$')
            ax.text(0.3, 0.35, f'{n_event} Events', fontsize=12, transform=ax.transAxes)
        else:
            ax.text(0.75, 0.35, f'{n_event} Events', fontsize=12, transform=ax.transAxes)

    fig.tight_layout()
    fig_abs.tight_layout()

    fig.savefig(f'{plot_out_dir}nexperiment_sd_dists.png')
    fig_abs.savefig(f'{plot_out_dir}nexperiment_absdev_dists.png')

    if show_plot:
        plt.show()


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
            plot_bootstraps(bs_list, stat, stat_vals[stat], stat_errs[stat], stats[stat]['true'], n_exp, out_dir)

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


if __name__ == '__main__':
    main()
