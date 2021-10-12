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

from SysReader import SysReader


def main():
    # data_set = ['raw', 'mix', 'pull_raw', 'pull_mix', 'divide', 'pull_divide']
    # data_set = ['divide', 'pull_divide']
    data_set = ['raw', 'mix', 'divide']
    dist_plts = ['single10']  # , 'poisson10']
    energy_plt = 62
    div_plt = [60, 90, 120, 240]
    cent_plt = 8
    s_plts = [0.002, 0.5]
    stats = ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    sys_path = '/home/dylan/Research/Results/Azimuth_Analysis/sys_vals_10-11-21_sim.txt'
    df = SysReader(sys_path).values

    ps = []
    ss = []
    dists = []
    for set_name in df['set']:
        dist, p, s = set_name.split('_')[1:]
        dists.append(dist)
        ps.append(float('0.' + p.strip('pgroup')) * 100)
        ss.append(float('0.' + s.strip('spread')))
    df['dist'] = dists
    df['pgroup'] = ps
    df['spread'] = ss

    for s_plt in s_plts:
        for dset in data_set:
            for dist_plt in dist_plts:
                plot(df, dset, div_plt, s_plt, dist_plt, stats)

    plt.show()

    print('donzo')


def plot(df, data_set, div_plt, s_plt, dist_plt, stats):
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

    # plt.show()


if __name__ == '__main__':
    main()
