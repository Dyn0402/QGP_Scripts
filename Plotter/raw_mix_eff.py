#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 10 10:46 PM 2021
Created in PyCharm
Created as QGP_Scripts/raw_mix_eff

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from SysReader import SysReader


def main():
    # data_types = ['raw', 'mix', 'pull_raw', 'pull_mix', 'divide', 'pull_divide']
    data_types = ['pull_raw', 'pull_mix']
    # set_pre = 'Sim_pgroup'
    # dist_plt = 'single10'
    data_sets_plt = ['AMPT', 'BES1']
    energy_plts = [7, 19, 39, 62]
    div_plts = [120]
    cent_plt = 8
    stats = ['mean', 'standard_deviation', 'skewness', 'non_excess_kurtosis']

    sys_path = '/home/dylan/Research/Results/Azimuth_Analysis/sys_vals_9-22-21_current.txt'
    df = SysReader(sys_path).values
    print(df.head())
    print(np.unique(df['set']))

    data_sets = []
    effs = []
    for set_name in df['set']:
        split = set_name.split('_')
        data_sets.append(split[0])
        if len(split) == 2 and 'Eff' in split[1]:
            ef = float('0.' + split[1].strip('Eff')) * 100
            if split[0] == 'BES1':
                ef = (1 - 0.68 * (1 - ef / 100)) * 100
            effs.append(ef)
        elif set_name == 'BES1':
            effs.append(32)
        else:
            effs.append(0)
    df['data_set'] = data_sets
    df['eff'] = effs

    for energy_plt in energy_plts:
        for div_plt in div_plts:
            for data_type in data_types:
                plot(df, cent_plt, div_plt, energy_plt, data_type, data_sets_plt, stats)

    plt.show()

    print('donzo')


def plot(df, cent, div, energy, data_type, data_sets, stats):
    fig, ax = plt.subplots(2, 2, sharex=True)
    c = {'BES1': 'r', 'AMPT': 'b'}

    for data_set in data_sets:
        df_ds = df[(df['data_set'] == data_set) & (df['div'] == div) & (df['data_type'] == data_type) &
                   (df['energy'] == energy) & (~df['set'].str.contains('rapid')) & (df['cent'] == cent)]

        for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
            df_stat = df_ds[df_ds['stat'] == stat]
            df_stat = df_stat.sort_values(by=['eff'])
            ax[ax_num].errorbar(df_stat['eff'], df_stat['val'], yerr=df_stat['stat_err'], marker='o', ls='none',
                                color=c[data_set], label=f'{data_set}')
            ax[ax_num].errorbar(df_stat['eff'], df_stat['val'], yerr=df_stat['sys_err'], marker='', ls='none',
                                color=c[data_set], elinewidth=4, alpha=0.4)

    for stat, ax_num in zip(stats, [(0, 0), (0, 1), (1, 0), (1, 1)]):
        y_lab = ax[ax_num].set_ylabel(stat)
        if ax_num[1] == 1:
            ax[ax_num].yaxis.tick_right()
            ax[ax_num].yaxis.set_label_position('right')
            y_lab.set_rotation(-90)
            ax[ax_num].yaxis.labelpad = 13
        if ax_num[0] == 1:
            ax[ax_num].set_xlabel('Efficiency (%)')
        ax[ax_num].grid()
    ax[(0, 0)].legend()
    fig.suptitle(f'{energy}GeV {data_type} Centrality {cent} {div}° Divisions')
    plt.subplots_adjust(wspace=0.05, hspace=0, left=0.1, right=0.9, top=0.9, bottom=0.1)


if __name__ == '__main__':
    main()
