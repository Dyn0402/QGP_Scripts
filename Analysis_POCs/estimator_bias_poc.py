#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 14 12:58 PM 2022
Created in PyCharm
Created as QGP_Scripts/estimator_bias_poc.py

@author: Dylan Neff, Dylan
"""

import pandas as pd
import os
from multiprocessing import Pool
import tqdm
import istarmap

from poc_functions import *


def main():
    print('donzo')


def estimator_bias():
    """
        Simulate binomials and test resampling values against known answer
        :return:
        """
    seed = 1432
    threads = 15
    n_tracks = [15]
    # n_samples = [1, 1440]
    n_samples = [1, 1440]
    # n_samples = np.array([1, 2, 3, 4, 5, 6, 7] + list(np.arange(10, 5000, 50)))
    # n_samples = np.arange(1, 31, 1)
    # n_events = np.arange(100, 2000, 5)
    n_events = np.arange(4, 31, 1)
    # n_events = [250]
    # bin_widths = np.deg2rad([60, 120, 240, 300])
    bin_widths = np.deg2rad([120])
    experiments = 1000
    # plot_out_dir = '/home/dylan/Research/Results/Resample_POC/nsample1440_nevent10000/'
    # plot_out_base = 'F:/Research/Resample_POC/Resample_Validation/'
    plot_out_base = 'E:/Transfer/Research/Resample_POC/Resample_Validation/'
    # plot_out_base = 'C:/Users/Dyn04/Desktop/Resample_POC/Resample_Validation/'
    # plot_out_name = 'vs_nsamplesto5k_bws_ntrack15_nevent250/'
    # plot_out_name = 'vs_nsamplesto20_bw180ish_ntrack15_nevent250/'
    # plot_out_name = 'vs_nsamplesto20_bw60and120_ntrack15and30_nevent2000/'
    # plot_out_name = 'vs_neventsto2k_bw60and120_ntrack15and30_nsamples1440/'
    plot_out_name = 'kstats_vs_nevento30_nsamples1and1440_bw120_ntrack15/'
    plot_out_dir = plot_out_base + plot_out_name
    plot_sds = False
    try:
        os.mkdir(plot_out_dir)
    except FileExistsError:
        pass
    show_plot = True

    stats = {bin_width: {n_track: define_stats(n_track, bin_width) for n_track in n_tracks}
             for bin_width in bin_widths}

    stats_plt = ['c4', 'k4', 'c4/c2', 'k4/k2']
    # stats_plt = ['standard deviation', 'skewness', 'non-excess kurtosis']
    # stats_plt = ['skewness', 'non-excess kurtosis']
    # stats_plt = ['standard deviation']

    write_info_file(plot_out_dir, threads, n_tracks, n_samples, n_events, bin_widths, 'estimator_bias',
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

    pars = {'n_events': len(n_events), 'n_samples': len(n_samples), 'bin_width': len(bin_widths)}
    indep_var = max(pars, key=pars.get)  # Use the longest parameter list as the independent variable
    plot_vs_indep_var(plot_data, stats_plt, stats, indep_var, plot_out_dir, plot_sds)
    # plot_vs_indep_var(plot_data, stats_plt, stats, 'n_events', plot_out_dir, plot_sds)
    # plot_vs_indep_var(plot_data, stats_plt, stats, 'n_samples', plot_out_dir, plot_sds)
    # plot_vs_indep_var(plot_data, stats_plt, stats, 'bin_width', plot_out_dir, plot_sds)

    if show_plot:
        plt.show()


def plot(plot_data, stats_plt, stats, indep_var, plot_out_dir, plot_sd=True):
    stat_combo_fig, stat_combo_ax = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    stat_combo_fig_del, stat_combo_ax_del = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    stat_combo_fig_del_norm, stat_combo_ax_del_norm = plt.subplots(len(stats_plt), 1, sharex=True, figsize=(8, 8))
    for stat_index, stat in enumerate(stats_plt):
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
        stat_combo_ax[stat_index].grid()
        stat_combo_ax_del[stat_index].grid()
        stat_combo_ax_del_norm[stat_index].grid()
        stat_combo_ax_del[stat_index].axhline(0, color='black')
        stat_combo_ax_del_norm[stat_index].axhline(0, color='black')

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
        num_set_combos = np.product([len(unique_vals) for unique_vals in uniques.values()])

        binom_val_label = False
        if indep_var == 'bin_width':
            for n_track in uniques['n_tracks']:
                binom_vals = []
                for bin_width in indep_vals:
                    binom_vals.append(stats[bin_width][n_track][stat]['true'])
                if not binom_val_label:
                    ax.plot(indep_vals, binom_vals, ls='--', color='black', label='Binomial Value')
                    stat_combo_ax[stat_index].plot(indep_vals, binom_vals, ls='--', color='black',
                                                   label='Binomial Value')
                    binom_val_label = True
                else:
                    ax.plot(indep_vals, binom_vals, ls='--', color='black')
                    stat_combo_ax[stat_index].plot(indep_vals, binom_vals, ls='--', color='black')
        elif indep_var == 'n_tracks':
            for bin_width in uniques['bin_width']:
                binom_vals = []
                for n_track in indep_vals:
                    binom_vals.append(stats[bin_width][n_track][stat]['true'])
                    if not binom_val_label:
                        ax.plot(indep_vals, binom_vals, ls='--', color='black', label='Binomial Value')
                        stat_combo_ax[stat_index].plot(indep_vals, binom_vals, ls='--', color='black',
                                                       label='Binomial Value')
                        binom_val_label = True
                    else:
                        ax.plot(indep_vals, binom_vals, ls='--', color='black')
                        stat_combo_ax[stat_index].plot(indep_vals, binom_vals, ls='--', color='black')
        else:
            for bin_width in uniques['bin_width']:
                for n_track in uniques['n_tracks']:
                    # label = f'{n_track} track, {round(np.rad2deg(bin_width))}° Binomial Value'
                    # c = next(color_binom)
                    if not binom_val_label:
                        ax.axhline(stats[bin_width][n_track][stat]['true'], ls='--', color='black',
                                   label='Binomial Value')
                        stat_combo_ax[stat_index].axhline(stats[bin_width][n_track][stat]['true'], ls='--',
                                                          color='black', label='Binomial Value')
                        binom_val_label = True
                    else:
                        ax.axhline(stats[bin_width][n_track][stat]['true'], ls='--', color='black')
                        stat_combo_ax[stat_index].axhline(stats[bin_width][n_track][stat]['true'], ls='--',
                                                          color='black')

        for set_var_vals in set_combos:
            set_var_vals = dict(zip(set_vars, set_var_vals))
            set_df = stat_df
            for var, var_val in set_var_vals.items():
                set_df = set_df[set_df[var] == var_val]

            c = next(color)
            means, sds, sems, deltas, delta_sems, delta_sds, delta_norms, delta_norm_sems, delta_norm_sds = \
                [], [], [], [], [], [], [], [], []
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
                delta_norms.append(deltas[-1] / binom_val)
                delta_sds.append(np.std(delts))
                delta_norm_sds.append(np.std(delts / binom_val))
                delta_sems.append(np.std(delts) / np.sqrt(vals.size))
                delta_norm_sems.append(np.std(delts / binom_val) / np.sqrt(vals.size))
                # print(f'{indep_val}: {delta_sds[-1]}')
            means, sds, sems, deltas, delta_sems, delta_sds, delta_norms, delta_norm_sems, delta_norm_sds = \
                (np.array(x) for x in (means, sds, sems, deltas, delta_sems, delta_sds, delta_norms, delta_norm_sems,
                                       delta_norm_sds))

            label = []
            for var in set_vars:
                if len(uniques[var]) > 1:
                    val = set_var_vals[var] if var != 'bin_width' else round(np.rad2deg(set_var_vals[var]))
                    label.append(f'{val}{var_string_consts[var]["legend"]}')
            label = ', '.join(label)

            for ax_ in [ax, stat_combo_ax[stat_index]]:
                ax_.plot(indep_vals, means, label=label, color=c)
                ax_.fill_between(indep_vals, means - sems, means + sems, color=c, alpha=0.6)
                if plot_sd:
                    ax_.fill_between(indep_vals, means - sds, means + sds, color=c, alpha=0.2)

            for ax_del_ in [ax_del, stat_combo_ax_del[stat_index]]:
                ax_del_.plot(indep_vals, deltas, label=label, color=c)
                ax_del_.fill_between(indep_vals, deltas - delta_sems, deltas + delta_sems, color=c, alpha=0.5)
                if plot_sd:
                    ax_del.fill_between(indep_vals, deltas - delta_sds, deltas + delta_sds, color=c, alpha=0.2)

            # del_max = max(deltas)
            for ax_del_norm_ in [ax_del_norm, stat_combo_ax_del_norm[stat_index]]:
                ax_del_norm_.plot(indep_vals, delta_norms, label=label, color=c)
                ax_del_norm_.fill_between(indep_vals, delta_norms - delta_norm_sems, delta_norms + delta_norm_sems,
                                          color=c, alpha=0.5)
                if plot_sd:
                    ax_del_norm.fill_between(indep_vals, delta_norms - delta_norm_sds, delta_norms + delta_norm_sds,
                                             color=c, alpha=0.2)

        ax.set_xlabel(var_string_consts[indep_var]['x-label'])
        ax_del.set_xlabel(var_string_consts[indep_var]['x-label'])
        ax_del_norm.set_xlabel(var_string_consts[indep_var]['x-label'])
        ax.legend()
        if num_set_combos > 1:
            ax_del.legend()
            ax_del_norm.legend()

        for ax_ in [stat_combo_ax[stat_index], stat_combo_ax_del[stat_index], stat_combo_ax_del_norm[stat_index]]:
            ax_.set_ylabel(stat)

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

        stat_str = stat.replace('/', '_')
        fig_names = {
            fig: f'{stat_str}_vs_{indep_var.replace("_", "")}',
            fig_del: f'{stat_str}_absdev_vs_{indep_var.replace("_", "")}',
            fig_del_norm: f'{stat_str}_absdev_norm_vs_{indep_var.replace("_", "")}',
        }

        for fig_obj, fig_name in fig_names.items():
            fig_obj.tight_layout()
            fig_obj.canvas.manager.set_window_title(fig_name)
            fig_obj.savefig(f'{plot_out_dir}{fig_name}.png', bbox_inches='tight')

    stat_combo_ax[-1].set_xlabel(var_string_consts[indep_var]['x-label'])
    stat_combo_ax_del[-1].set_xlabel(var_string_consts[indep_var]['x-label'])
    stat_combo_ax_del_norm[-1].set_xlabel(var_string_consts[indep_var]['x-label'])
    stat_combo_ax[0].legend()
    if num_set_combos > 1:
        stat_combo_ax_del[0].legend()
        stat_combo_ax_del_norm[0].legend()

    title = f'Stats vs {var_string_consts[indep_var]["x-label"]}'
    title_del = f'Stats Deviations vs {var_string_consts[indep_var]["x-label"]}'
    title_del_norm = f'Stats Normalized Deviations vs {var_string_consts[indep_var]["x-label"]}'
    title_extra = [f'{vals[0]}{var_string_consts[var]["legend"]}' if var != 'bin_width'
                   else f'{round(np.rad2deg(vals[0]))}{var_string_consts[var]["legend"]}' for var, vals
                   in uniques.items() if len(vals) == 1]
    title_extra = ', '.join(title_extra)
    if title_extra != '':
        title += f'\n{title_extra}'
        title_del += f'\n{title_extra}'
        title_del_norm += f'\n{title_extra}'

    stat_combo_fig.suptitle(title)
    stat_combo_fig_del.suptitle(title_del)
    stat_combo_fig_del_norm.suptitle(title_del_norm)

    fig_combo_names = {
        stat_combo_fig: f'stats_vs_{indep_var.replace("_", "")}',
        stat_combo_fig_del: f'stats_absdev_vs_{indep_var.replace("_", "")}',
        stat_combo_fig_del_norm: f'stats_absdev_norm_vs_{indep_var.replace("_", "")}',
    }

    for fig_obj, fig_name in fig_combo_names.items():
        fig_obj.tight_layout()
        # fig_obj.subplots_adjust(hspace=0)
        fig_obj.canvas.manager.set_window_title(fig_name)
        fig_obj.savefig(f'{plot_out_dir}{fig_name}.png', bbox_inches='tight')


if __name__ == '__main__':
    main()
