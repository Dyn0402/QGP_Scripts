#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 10 5:13 PM 2021
Created in PyCharm
Created as QGP_Scripts/sim_amp_spread_plot

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from SysReader import SysReader

from sim_pgroup_plot import add_pgrp_spread


def main():
    pd.options.mode.chained_assignment = None
    sys_path = '/home/dylan/Research/Results/Azimuth_Analysis/sys_vals_11-10-21_clustmulti_sim.txt'
    # sys_path = 'C:\\Users\\Dyn04\\Downloads\\sys_vals_11-1-21_anticlust_sim.txt'
    df = SysReader(sys_path).values

    df = add_pgrp_spread(df)

    print(df.head())
    print(df.columns)

    plot1(df)
    plot2(df)
    plot3(df)

    plt.show()

    print('donzo')


def plot1(df):
    # data_set = ['raw', 'mix', 'pull_raw', 'pull_mix', 'divide', 'pull_divide']
    data_set = ['divide']
    dist_plts = ['single10']  # , 'poisson10']
    div_plt = [60, 90, 120, 240, 300]
    s_plts = [0.3]
    stats = ['standard_deviation', 'skewness', 'non_excess_kurtosis']

    for s_plt in s_plts:
        for dset in data_set:
            for dist_plt in dist_plts:
                plot_vs_amp_3sub(df, dset, div_plt, s_plt, dist_plt, stats)

    # plt.show()


def plot2(df):
    data_set = ['divide']
    dist_plts = ['single10']  # , 'poisson10']
    div_plt = [60, 90, 120, 240, 300]
    amp_plts = [0.1]
    stats = ['standard_deviation', 'skewness', 'non_excess_kurtosis']

    for amp_plt in amp_plts:
        for dset in data_set:
            for dist_plt in dist_plts:
                plot_vs_spread_3sub(df, dset, div_plt, amp_plt, dist_plt, stats)

    # plt.show()


def plot3(df):
    data_set = ['raw', 'mix', 'divide', 'pull_divide']
    dist_plts = ['single']  # , 'poisson10']
    div_plt = [60, 90, 120, 240, 300]
    amp_plts = [0.05]
    s_plts = [0.1]
    stats = ['standard_deviation', 'skewness', 'non_excess_kurtosis']

    for amp_plt in amp_plts:
        for s_plt in s_plts:
            for dset in data_set:
                for dist_plt in dist_plts:
                    plot_vs_nevent_3sub(df, dset, div_plt, amp_plt, s_plt, dist_plt, stats)


def plot_vs_amp_3sub(df, data_set, div_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(3, 1, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['spread'] == s_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [0, 1, 2]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['amp'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['amp'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['amp'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)
            if 'divide' in data_set:
                ax[ax_num].axhline(1, color='black', ls='--', alpha=0.6, zorder=5)

    for stat, ax_num in zip(stats, [0, 1, 2]):
        ax[ax_num].set_ylabel(stat)
        if ax_num == 2:
            ax[ax_num].set_xlabel('Gaussian Kernel Amplitude')
        ax[ax_num].grid()
    ax[0].legend(ncol=2)
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, right=0.99, top=0.9, bottom=0.1)


def plot_vs_spread_3sub(df, data_set, div_plt, amp_plt, dist_plt, stats):
    fig, ax = plt.subplots(3, 1, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['amp'] == amp_plt) &
                    (df['dist'] == dist_plt)]

        for stat, ax_num in zip(stats, [0, 1, 2]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['spread'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['spread'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['spread'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)
            if 'divide' in data_set:
                ax[ax_num].axhline(1, color='black', ls='--', alpha=0.6, zorder=5)

    for stat, ax_num in zip(stats, [0, 1, 2]):
        ax[ax_num].set_ylabel(stat)
        if ax_num == 2:
            ax[ax_num].set_xlabel('Gaussian Kernel Sigma')
        ax[ax_num].grid()
    ax[0].legend(ncol=2)
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {amp_plt} amp')
    plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, right=0.99, top=0.9, bottom=0.1)


def plot_vs_nevent_3sub(df, data_set, div_plt, amp_plt, s_plt, dist_plt, stats):
    fig, ax = plt.subplots(3, 1, sharex=True)

    for div in div_plt:
        df_div = df[(df['data_type'] == data_set) & (df['div'] == div) & (df['amp'] == amp_plt) &
                    (df['spread'] == s_plt) & (df['dist'].str.contains(dist_plt))]

        for stat, ax_num in zip(stats, [0, 1, 2]):
            df_stat = df_div[df_div['stat'] == stat]
            df_stat = df_stat.sort_values(by=['nevent'])
            # ax[ax_num].errorbar(df_stat['pgroup'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
            #                     label=f'{div}°')
            ax[ax_num].fill_between(df_stat['nevent'], df_stat['val'] + df_stat['stat_err'],
                                    df_stat['val'] - df_stat['stat_err'], label=f'{div}°', alpha=0.7)
            ax[ax_num].errorbar(df_stat['nevent'], df_stat['val'], yerr=df_stat['sys_err'], marker='.', ls='none',
                                elinewidth=4, alpha=0.7)
            if 'divide' in data_set:
                ax[ax_num].axhline(1, color='black', ls='--', alpha=0.6, zorder=5)

    for stat, ax_num in zip(stats, [0, 1, 2]):
        ax[ax_num].set_ylabel(stat)
        if ax_num == 2:
            ax[ax_num].set_xlabel('Number of Protons per Event')
        ax[ax_num].grid()
    ax[0].legend(ncol=2)
    fig.suptitle(f'{data_set} {dist_plt} {div_plt}° Divisions {amp_plt} amp {s_plt} spread')
    plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.1, right=0.99, top=0.9, bottom=0.1)


if __name__ == '__main__':
    main()
