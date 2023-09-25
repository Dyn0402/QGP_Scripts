#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 21 10:47 AM 2022
Created in PyCharm
Created as QGP_Scripts/visualizations.py

@author: Dylan Neff, Dylan
"""

import os
from multiprocessing import Pool

import pandas as pd
import tqdm
import istarmap

from poc_functions import *


def main():
    # correlated_dists()
    multiple_exps_dist()
    print('donzo')


def correlated_dists():
    threads = 14
    seed_0 = 1432
    n_tracks = 15
    n_samples = 1
    n_events_sim = np.arange(100, 1e5 + 1, 1000, dtype=int)
    n_events_dist_plot = np.array([1e2, 1e3, 1e5], dtype=int)
    bin_width = np.deg2rad(120)
    stat_plt = 'variance'
    # plot_out_base = 'E:/Transfer/Research/Resample_POC/Visualizations/'
    plot_out_base = 'F:/Research/Resample_POC/Visualizations/'
    plot_out_name = 'test3/'
    plot_out_dir = plot_out_base + plot_out_name
    # plt.rcParams['figure.figsize'] = (8, 4)
    plt.rcParams['figure.dpi'] = 144
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
        fig, axs = plt.subplots(len(n_events_dist_plot), 1, sharex=True, figsize=(8, 5))
        event_axes = dict(zip(n_events_dist_plot, axs))
        plt.subplots_adjust(hspace=0)

        fig_stat, ax_stat = plt.subplots(figsize=(8, 3))
        ax_stat.set_xlabel('Number of Events')
        ax_stat.set_ylabel(stat_plt)

        # fig_dev, ax_dev = plt.subplots()
        # ax_dev.set_xlabel('Number of Events')
        # ax_dev.set_ylabel(f'{stat_plt} deviation from binomial')

        fig_dev_abs, ax_dev_abs = plt.subplots(figsize=(8, 3))
        ax_dev_abs.set_xlabel('Number of Events')
        ax_dev_abs.set_ylabel(r'$|\sigma^2_{data} - \sigma^2_{binomial}|$')

        jobs = []
        for n_event, seed in zip(n_events_sim, seeds):
            jobs.append((n_tracks, n_event, bin_width, n_samples, stats, [stat_plt], seed))

            if n_event in n_events_dist_plot:
                rng = np.random.default_rng(seed)
                experiment = gen_experiment(n_event, n_tracks, rng)
                hist = bin_experiment_no_bs(experiment, n_tracks, bin_width, n_samples, rng)
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
            for n_exp, samples, n_events, bin_width, n_track, stats_vals, stats_errs_delta, alg in \
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
            fig_stat: f'var_vs_events_{exp_name}',
            # fig_dev: f'var_dev2_vs_events_{exp_name}',
            fig_dev_abs: f'var_devabs_vs_events_{exp_name}',
        }

        for fig_obj, fig_name in fig_names.items():
            fig_obj.canvas.manager.set_window_title(fig_name)
            fig_obj.tight_layout()
            if plot_out_dir is not None:
                fig_obj.savefig(f'{plot_out_dir}{fig_name}.png')
                fig_obj.savefig(f'{plot_out_dir}{fig_name}.pdf')

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
    stat_plt = 'variance'
    # plot_out_base = 'E:/Transfer/Research/Resample_POC/Visualizations/'
    plot_out_base = 'F:/Research/Resample_POC/Visualizations/'
    plot_out_name = 'test4/'
    plot_out_dir = plot_out_base + plot_out_name
    plt.rcParams['figure.figsize'] = (8, 5)
    plt.rcParams['figure.dpi'] = 144
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
        for n_exp, n_sample, n_event, bin_width, n_track, stats_vals, stats_errs_delta, alg in \
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
            ax.set_xlabel('Variance of Distribution')
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
            ax.set_xlabel(r'$|\sigma^2 - \sigma^2_{binom}|$')
            ax.text(0.3, 0.35, f'{n_event} Events', fontsize=12, transform=ax.transAxes)
        else:
            ax.text(0.75, 0.35, f'{n_event} Events', fontsize=12, transform=ax.transAxes)

    fig.tight_layout()
    fig_abs.tight_layout()

    fig.savefig(f'{plot_out_dir}nexperiment_var_dists.png')
    fig.savefig(f'{plot_out_dir}nexperiment_var_dists.pdf')
    fig_abs.savefig(f'{plot_out_dir}nexperiment_absdev_dists.png')
    fig_abs.savefig(f'{plot_out_dir}nexperiment_absdev_dists.pdf')

    if show_plot:
        plt.show()


def single_event_nsample_convergance():
    # Mostly deprecated for sub_event_resample_algorithm animate_nsamples_resamples3
    threads = 15
    seed_0 = 1432
    n_tracks = 15
    n_samples = np.arange(4, 5000, 1)
    n_events = 1
    n_experiments = 2
    bin_width = np.deg2rad(120)
    stats_plt = ['standard deviation']
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
    seeds = iter(np.random.SeedSequence(seed_0).spawn(n_experiments))

    fig, axs = plt.subplots(len(stats_plt), 1, sharex=True)
    axs[-1].set_xlabel('Number of Samples per Event')

    jobs = [(n_tracks, n_events, bin_width, n_samples, stats, stats_plt, next(seeds), n_exp)
            for n_exp in range(n_experiments)]

    plot_data = []
    with Pool(threads) as pool:
        for n_exp, n_sample, n_event, bin_width, n_track, stats_vals, stats_errs_delta, alg in \
                tqdm.tqdm(pool.istarmap(run_experiment_no_bs, jobs), total=len(jobs)):
            for stat, val in stats_vals.items():
                plot_data.append({'n_exp': n_exp, 'n_sample': n_sample, 'stat': stat, 'val': val})
    plot_data = pd.DataFrame(plot_data)

    for stat_index, stat in enumerate(stats_plt):
        stat_df = plot_data[plot_data['stat'] == stat]
        axs[stat_index].set_ylabel(stat.capitalize())
        n_exps_unq = pd.unique(stat_df['n_exp'])
        for n_exp in n_exps_unq:
            exp_df = stat_df['stat' == stat]
            axs[stat_index].scatter(exp_df['n_sample'], exp_df['val'], label=f'Experiment #{n_exp}')

    if n_experiments > 1:
        axs[0].legend()

    fig.tight_layout()
    fig.savefig(f'{plot_out_dir}single_event_nsample_convergence.png')
    fig.savefig(f'{plot_out_dir}single_event_nsample_convergence.pdf')

    if show_plot:
        plt.show()


if __name__ == '__main__':
    main()
