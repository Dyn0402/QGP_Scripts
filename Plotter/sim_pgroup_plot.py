#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 09 3:55 PM 2021
Created in PyCharm
Created as QGP_Scripts/sim_pgroup_plot

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from SysReader import SysReader


def main():
    pd.options.mode.chained_assignment = None
    sys_path = '/home/dylan/Research/Results/Azimuth_Analysis/sys_vals_11-1-21_anticlust_sim.txt'
    # sys_path = 'C:\\Users\\Dyn04\\Downloads\\sys_vals_11-1-21_anticlust_sim.txt'
    df = SysReader(sys_path).values

    df = add_pgrp_spread(df)

    # plot1(df)
    # plot2(df)
    plot3(df)
    # plot4(df)

    print('donzo')


def plot1(df):
    # data_set = ['raw', 'mix', 'pull_raw', 'pull_mix', 'divide', 'pull_divide']
    data_set = ['raw', 'mix', 'divide', 'pull_divide']
    dist_plts = ['single10']  # , 'poisson10']
    div_plt = [60, 90, 120, 240, 300]
    s_plts = [0.5]  # 0.002, 0.5
    stats = ['standard_deviation', 'skewness', 'non_excess_kurtosis']  # ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    for s_plt in s_plts:
        for dset in data_set:
            for dist_plt in dist_plts:
                plot_vs_pgroup_3sub(df, dset, div_plt, s_plt, dist_plt, stats)
                if 'divide' in data_set:
                    # plot_nsigma1_vs_pgroup(df, dset, div_plt, s_plt, dist_plt, stats)
                    plot_nsigma1_vs_pgroup_3sub(df, dset, div_plt, s_plt, dist_plt, stats)

    plt.show()


def plot2(df):
    data_set = ['raw', 'mix', 'divide']
    dist_plts = ['single48']
    div_plt = [60, 90, 120, 240, 300]
    s_plts = [0.002]
    stats = ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    for s_plt in s_plts:
        for dset in data_set:
            for dist_plt in dist_plts:
                plot_vs_pgroup(df, dset, div_plt, s_plt, dist_plt, stats)

    plt.show()


def plot3(df):
    data_set = ['raw', 'mix', 'divide', 'pull_divide']
    dist_plt = 'single'
    div_plt = [60, 90, 120, 240, 300]
    s_plts = [0.5]
    pgroup = 0.1
    stats = ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    print(np.unique(df['pgroup']))

    for s_plt in s_plts:
        for dset in data_set:
            plot_vs_nevent(df, dset, div_plt, s_plt, dist_plt, pgroup, stats)

    plt.show()


def plot4(df):
    data_set = ['raw', 'mix', 'divide', 'pull_divide']
    dist_plt = 'poisson'
    div_plt = [60, 90, 120, 240, 300]
    s_plts = [0.002]
    pgroup = 2
    stats = ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    print(np.unique(df['pgroup']))

    for s_plt in s_plts:
        for dset in data_set:
            plot_vs_nevent(df, dset, div_plt, s_plt, dist_plt, pgroup, stats)

    plt.show()


def add_pgrp_spread(df):
    pgroups = []
    ss = []
    dists = []
    nevent = []
    for set_name in df['set']:
        try:
            dist, amp, s = set_name.split('_')[1:]
        except ValueError:
            try:
                dist, anticlust, amp, s = set_name.split('_')[1:]
            except ValueError:
                print("Can't read name format")
        dists.append(dist)
        nevent.append(int(dist.strip('single').strip('poisson')))
        pgroups.append(float('0.' + amp.strip('pgroup')))
        ss.append(float('0.' + s.strip('spread')))
    df['dist'] = dists
    df['pgroup'] = pgroups
    df['spread'] = ss
    df['nevent'] = nevent

    return df


def plot_vs_pgroup(df, data_set, div_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(2, 2, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['pgroup'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['pgroup'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)

    for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        y_lab = ax[ax_num].set_ylabel(stat)
        if ax_num[1] == 1:
            ax[ax_num].yaxis.tick_right()
            ax[ax_num].yaxis.set_label_position('right')
            y_lab.set_rotation(-90)
            ax[ax_num].yaxis.labelpad = 13
        if ax_num[0] == 1:
            ax[ax_num].set_xlabel('Probability of Grouping (%)')
        ax[ax_num].grid()
    ax[(0, 0)].legend()
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0, left=0.1, right=0.9, top=0.9, bottom=0.1)


def plot_vs_pgroup_3sub(df, data_set, div_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(3, 1, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [0, 1, 2]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['pgroup'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['pgroup'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)
            if 'divide' in data_set:
                ax[ax_num].axhline(1, color='black', ls='--', alpha=0.6, zorder=5)

    for stat, ax_num in zip(stats, [0, 1, 2]):
        ax[ax_num].set_ylabel(stat)
        if ax_num == 2:
            ax[ax_num].set_xlabel('Probability of Grouping (%)')
        ax[ax_num].grid()
    ax[0].legend(ncol=2)
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, right=0.99, top=0.9, bottom=0.1)


def plot_vs_nevent(df, data_set, div_plt, s_plt, dist_plt, pgroup, stats):
    fig, ax = plt.subplots(2, 2, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['pgroup'] == pgroup) & (df['dist'].str.contains(dist_plt))]
        # df_div['nevent'] = df_div['dist'].str.strip(dist_plt).astype('int32')
        # df_div.assign(nevent=df_div['dist'].str.strip(dist_plt).astype('int32').values)

        for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['nevent'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['nevent'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['nevent'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)

    for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        y_lab = ax[ax_num].set_ylabel(stat)
        if ax_num[1] == 1:
            ax[ax_num].yaxis.tick_right()
            ax[ax_num].yaxis.set_label_position('right')
            y_lab.set_rotation(-90)
            ax[ax_num].yaxis.labelpad = 13
        if ax_num[0] == 1:
            ax[ax_num].set_xlabel('Mean Number of Protons per Event')
        ax[ax_num].grid()
    ax[(0, 0)].legend()
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread {pgroup}% pgroup')
    plt.subplots_adjust(wspace=0.05, hspace=0, left=0.1, right=0.9, top=0.9, bottom=0.1)


def plot_nsigma1_vs_pgroup(df, data_set, div_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(2, 2, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['pgroup'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].plot(df_stat['pgroup'], (df_stat['val'] - 1) / df_stat['sys_err'], marker='.',
                               label=f'{div}°', alpha=0.7)

    for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        y_lab = ax[ax_num].set_ylabel(stat)
        if ax_num[1] == 1:
            ax[ax_num].yaxis.tick_right()
            ax[ax_num].yaxis.set_label_position('right')
            y_lab.set_rotation(-90)
            ax[ax_num].yaxis.labelpad = 13
        if ax_num[0] == 1:
            ax[ax_num].set_xlabel('Probability of Grouping (%)')
        ax[ax_num].grid()
    ax[(0, 0)].legend()
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0, left=0.1, right=0.9, top=0.9, bottom=0.1)


def plot_nsigma1_vs_pgroup_3sub(df, data_set, div_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(3, 1, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [0, 1, 2]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['pgroup'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].plot(df_stat['pgroup'], (df_stat['val'] - 1) / df_stat['sys_err'], marker='.',
                               label=f'{div}°', alpha=0.7)
            if 'divide' in data_set:
                ax[ax_num].axhline(1, color='black', ls='--', alpha=0.6, zorder=5)

    for stat, ax_num in zip(stats, [0, 1, 2]):
        ax[ax_num].set_ylabel(stat)
        if ax_num == 2:
            ax[ax_num].set_xlabel('Probability of Grouping (%)')
        ax[ax_num].grid()
    ax[0].legend(ncol=2)
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, right=0.99, top=0.9, bottom=0.1)


def plot_nsigma1_vs_nevent(df, data_set, div_plt, s_plt, dist_plt, pgroup, stats):
    fig, ax = plt.subplots(2, 2, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['pgroup'] == pgroup) & (df['dist'].str.contains(dist_plt))]
        # df_div['nevent'] = df_div['dist'].str.strip(dist_plt).astype('int32')
        # df_div.assign(nevent=df_div['dist'].str.strip(dist_plt).astype('int32').values)

        for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['nevent'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].plot(df_stat['nevent'], (df_stat['val'] - 1) / df_stat['sys_err'], marker='.',
                               label=f'{div}°', alpha=0.7)

    for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        y_lab = ax[ax_num].set_ylabel(stat)
        if ax_num[1] == 1:
            ax[ax_num].yaxis.tick_right()
            ax[ax_num].yaxis.set_label_position('right')
            y_lab.set_rotation(-90)
            ax[ax_num].yaxis.labelpad = 13
        if ax_num[0] == 1:
            ax[ax_num].set_xlabel('Mean Number of Protons per Event')
        ax[ax_num].grid()
    ax[(0, 0)].legend()
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread {pgroup}% pgroup')
    plt.subplots_adjust(wspace=0.05, hspace=0, left=0.1, right=0.9, top=0.9, bottom=0.1)


if __name__ == '__main__':
    main()
