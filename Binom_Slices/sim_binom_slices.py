#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 15 1:30 PM 2021
Created in PyCharm
Created as QGP_Scripts/sim_binom_slices

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.optimize import curve_fit as cf

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats


def main():
    base_path = '/home/dylan/Research/'
    # base_path = 'D:/Transfer/Research/'
    df_path = '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # div = 120
    read_data = False  # Read data from text files instead of reading dataframe
    divs = [60, 72, 89, 90, 120, 240, 270, 288, 300]
    cents = [8]
    sim_pars = [('amp05', 'spread01'), ('amp1', 'spread01'), ('amp05', 'spread03'), ('amp1', 'spread02'),
                ('amp1', 'spread03'), ('amp15', 'spread01'), ('amp15', 'spread02'), ('amp15', 'spread03'),
                ('amp05', 'spread02'),
                ('amp01', 'spread1'), ('amp02', 'spread1'), ('amp03', 'spread1'), ('amp04', 'spread1'),
                ('amp01', 'spread15'), ('amp02', 'spread15'), ('amp03', 'spread15'), ('amp04', 'spread15'),
                ('amp01', 'spread08'), ('amp02', 'spread08'), ('amp03', 'spread08'), ('amp04', 'spread08'),
                ('amp01', 'spread2'), ('amp02', 'spread2'), ('amp03', 'spread2'), ('amp04', 'spread2')]
    stats = {'mean': getattr(DistStats, 'get_mean'), 'standard deviation': getattr(DistStats, 'get_sd'),
             'skewness': getattr(DistStats, 'get_skewness'), 'kurtosis': getattr(DistStats, 'get_kurtosis')}
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.05), 'skewness': (0, 1.25), 'kurtosis': (-3, 2)}
    stats_plot = ['standard deviation']  # , 'skewness', 'kurtosis']
    div_plt = 60
    total_protons_plt = 32
    cent_plt = 8
    data_type_plt = 'divide'
    data_sets_plt = ['sim_amp02_spread15']

    if read_data:
        df = get_dist_stats(base_path, divs, cents, sim_pars, stats)
        df.to_csv(df_path)
    else:
        df = pd.read_csv(df_path)

    df = df.dropna()

    # plot_vs_protons(df, stats_plot, div_plt, cent_plt, data_type_plt, y_ranges)
    # plot_vs_divs(df, stats_plot, total_protons_plt, cent_plt, data_type_plt, y_ranges)
    for div_plt in [240]:
        fit_vs_protons(df, stats_plot, div_plt, cent_plt, data_type_plt, y_ranges, data_sets_plt)
    # fit_vs_divs(df, stats_plot, total_protons_plt, cent_plt, data_type_plt, y_ranges, data_sets_plt)
    plt.show()
    print('donzo')


def get_dist_stats(base_path, divs, cents, sim_pars, stats):
    sim_energy = 62
    ampt_energy = 7

    ampt_set = 'default/Ampt_rapid05_n1ratios_'

    sim_dfs = []
    ampt_dfs = []
    for div in divs:
        for cent in cents:
            for pars in sim_pars:
                print(f'Get sim {pars}, {div} div, {cent} cent')
                sim_dist = get_sim_dists(base_path + 'Data_Sim/', range(11), div, cent, sim_energy,
                                         exact_keys=['anticlmulti', *pars], contain_keys=['single'])
                sim_dist_mix = get_sim_dists(base_path + 'Data_Sim_Mix/', range(11), div, cent, sim_energy,
                                             exact_keys=['anticlmulti', *pars], contain_keys=['single'])
                stats_dict = get_stats(sim_dist, sim_dist_mix, stats)
                pars_dict = {'div': div, 'cent': cent, 'energy': sim_energy, 'amp': pars[0], 'spread': pars[1],
                             'name': f'sim_{pars[0]}_{pars[1]}'}
                sim_dfs.append(pd.DataFrame([{**entry, **pars_dict} for entry in stats_dict]))

            print(f'Get ampt, {div} div, {cent} cent')
            ampt_dist = get_ampt_dist(base_path + 'Data_Ampt/' + ampt_set, div, cent, ampt_energy, range(60))
            ampt_dist_mix = get_ampt_dist(base_path + 'Data_Ampt_Mix/' + ampt_set, div, cent, ampt_energy, range(60))
            stats_dict = get_stats(ampt_dist, ampt_dist_mix, stats)
            pars_dict = {'div': div, 'cent': cent, 'energy': ampt_energy, 'amp': 'ampt', 'spread': 'ampt',
                         'name': 'ampt'}
            ampt_dfs.append(pd.DataFrame([{**entry, **pars_dict} for entry in stats_dict]))

    df = pd.concat(sim_dfs, ignore_index=True)
    df = df.append(pd.concat(ampt_dfs, ignore_index=True))

    return df


def get_stats(raw_dists, mix_dists, stat_names):
    data_types = ['raw', 'mix', 'divide']
    stat_meases = {stat_name: {data_type: {} for data_type in data_types} for stat_name in stat_names}
    for set_num in mix_dists:
        for total_protons in sorted(mix_dists[set_num]):
            if total_protons in raw_dists[set_num]:
                raw_stats = DistStats(raw_dists[set_num][total_protons])
                mix_stats = DistStats(mix_dists[set_num][total_protons])
                for stat_name, stat_meth in stat_names.items():
                    raw_meas = stat_meth(raw_stats)
                    mix_meas = stat_meth(mix_stats)
                    divide_meas = raw_meas / mix_meas
                    for data_type, meas in zip(data_types, [raw_meas, mix_meas, divide_meas]):
                        # print(stat_meth(raw_stats), stat_meth(mix_stats))
                        # stat_meas = stat_meth(raw_stats) / stat_meth(mix_stats)
                        if total_protons in stat_meases[stat_name][data_type]:
                            stat_meases[stat_name][data_type][total_protons].append(meas)
                        else:
                            stat_meases[stat_name][data_type].update({total_protons: [meas]})

    data = []
    num_sets = len(mix_dists)
    for stat in stat_meases:
        for data_type in data_types:
            for total_protons in sorted(stat_meases[stat][data_type]):
                meases = stat_meases[stat][data_type][total_protons]
                if len(meases) == num_sets:  # Only save if all sets are present
                    med = np.median(meases)
                    data.append({'stat': stat, 'data_type': data_type, 'total_protons': total_protons,
                                 'val': med.val, 'err': med.err, 'sys': np.std(meases).val})

    return data


def get_ampt_dist(path, div, cent, energy, set_nums):
    dists = {}
    for set_num in set_nums:
        txt_path = f'{path}{set_num}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
        dists[set_num] = AzimuthBinData(div, txt_path).get_dist()

    return dists


def get_sim_dists(path, set_nums, div, cent, energy, exact_keys=[], contain_keys=[]):
    dists = {set_num: {} for set_num in set_nums}
    for dir_name in os.listdir(path):
        keep_dir = True
        for key in exact_keys:
            if key not in dir_name.split('_'):
                keep_dir = False
                break
        for key in contain_keys:
            if key not in dir_name:
                keep_dir = False
                break
        if keep_dir:
            for set_name in os.listdir(path + dir_name):
                set_num = int(set_name.split('_')[-1])
                if set_num in set_nums:
                    txt_path = f'{path}{dir_name}/' \
                               f'{set_name}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
                    data = AzimuthBinData(div, txt_path)
                    new_dist = data.get_dist()
                    for total_protons in new_dist:
                        dists[set_num][total_protons] = new_dist[total_protons]

    return dists


def line(x, a, b):
    return a * x + b


def quad_180(x, a, c):
    return a * (x - 180) ** 2 + c


def plot_vs_protons(df, stat_names, div, cent, data_type, y_ranges):
    # div = 60
    # cent = 8
    # data_type = 'divide'
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(stat)
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])

        for data_set in np.unique(df['name']):
            df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['div'] == div) &
                        (df['cent'] == cent) & (df['stat'] == stat)]
            df_set = df_set.sort_values(by=['total_protons'])
            if 'ampt' in data_set:
                ax.errorbar(df_set['total_protons'], df_set['val'], df_set['err'], label=data_set, marker='o', ls='',
                            color='blue')
                ax.errorbar(df_set['total_protons'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                            color='blue', alpha=0.3)
            else:
                ax.fill_between(df_set['total_protons'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                                label=data_set)
        ax.legend()
        fig.tight_layout()

    # plt.show()


def plot_vs_divs(df, stat_names, total_protons, cent, data_type, y_ranges):
    # cent = 8
    # data_type = 'divide'
    # total_protons = 30
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(stat)
        ax.set_xlabel('Division Width')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])

        for data_set in np.unique(df['name']):
            df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) &
                        (df['total_protons'] == total_protons) & (df['cent'] == cent) & (df['stat'] == stat)]
            df_set = df_set.sort_values(by=['div'])
            if 'ampt' in data_set:
                ax.errorbar(df_set['div'], df_set['val'], df_set['err'], label=data_set, marker='o', ls='',
                            color='blue')
                ax.errorbar(df_set['div'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                            color='blue', alpha=0.3)
            else:
                ax.fill_between(df_set['div'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                                label=data_set)
        ax.legend()
        fig.tight_layout()

    # plt.show()


def fit_vs_protons(df, stat_names, div, cent, data_type, y_ranges, data_set_plt):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {div}° Divisions')
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(np.unique(df['name'])))))

        for data_set in np.unique(df['name']):
            c = next(color)
            df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['div'] == div) &
                        (df['cent'] == cent) & (df['stat'] == stat)]
            df_set = df_set.sort_values(by=['total_protons'])
            if 'ampt' in data_set:
                ax.errorbar(df_set['total_protons'], df_set['val'], df_set['err'], label=data_set, marker='o', ls='',
                            color=c, alpha=0.8)
                ax.errorbar(df_set['total_protons'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                            color=c, alpha=0.2)
            elif data_set in data_set_plt:
                ax.fill_between(df_set['total_protons'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                                label=data_set, color=c, alpha=0.8)
            if len(df_set['val']) > 1:
                # print(f'df_set[total_protons]: {df_set["total_protons"]}, df_set[val]: {df_set["val"]}')
                # Why are some linear fits failing?
                popt, pcov = cf(line, df_set['total_protons'], df_set['val'], sigma=df_set['err'], absolute_sigma=True)
                if 'ampt' in data_set or data_set in data_set_plt:
                    print(df_set['total_protons'])
                    print(df_set['val'])
                    print(df_set['err'])
                    ax.plot(df_set['total_protons'], line(df_set['total_protons'], *popt), ls='--', color=c)
        ax.legend()
        fig.tight_layout()


def fit_vs_divs(df, stat_names, total_protons, cent, data_type, y_ranges, data_set_plt):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {total_protons} protons')
        ax.set_xlabel('Division Width')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])

        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(np.unique(df['name'])))))
        fit_pars = []
        for data_set in np.unique(df['name']):
            c = next(color)
            df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) &
                        (df['total_protons'] == total_protons) & (df['cent'] == cent) & (df['stat'] == stat)]
            df_set = df_set.sort_values(by=['div'])
            if 'ampt' in data_set:
                ax.errorbar(df_set['div'], df_set['val'], df_set['err'], label=data_set, marker='o', ls='',
                            color=c, alpha=0.8)
                ax.errorbar(df_set['div'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                            color=c, alpha=0.3)
            elif data_set in data_set_plt:
                ax.fill_between(df_set['div'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                                label=data_set, color=c, alpha=0.8)
            if len(df_set['val'] > 1):
                popt, pcov = cf(quad_180, df_set['div'], df_set['val'], sigma=df_set['err'], absolute_sigma=True)
                perr = np.sqrt(np.diag(pcov))
                fit_pars.append({'data_set': data_set, 'curvature': popt[0], 'baseline': popt[1], 'color': c,
                                 'spread': df_set['spread'].iloc[0], 'amp': df_set['amp'].iloc[0],
                                 'curve_err': perr[0], 'base_err': perr[1]})
                x = np.linspace(0, 360, 100)
                if 'ampt' in data_set or data_set in data_set_plt:
                    ax.plot(x, quad_180(x, *popt), ls='--', color=c)
        ax.legend()
        fig.tight_layout()

        df_pars = pd.DataFrame(fit_pars)
        fig2, ax2 = plt.subplots()
        ax2.set_title(stat)
        ax2.set_xlabel('baseline')
        ax2.set_ylabel('curvature')
        spreads = np.unique(df_pars['spread'])
        print(spreads)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads))))
        for spread in spreads:
            c = next(color)
            df_spread = df_pars[df_pars['spread'] == spread]
            if len(df_spread['curvature']) > 1:
                popt, pcov = cf(line, df_spread['baseline'], df_spread['curvature'], sigma=df_spread['curve_err'],
                                absolute_sigma=True)
                x_plt = np.linspace(y_ranges[stat][0], 1, 3)
                ax2.plot(x_plt, line(x_plt, *popt), ls='--', alpha=0.7, color=c)
            ax2.errorbar(df_spread['baseline'], df_spread['curvature'], xerr=df_spread['base_err'], marker='o',
                         yerr=df_spread['curve_err'], color=c, label=spread, ls='none', alpha=0.6)
        ax2.legend()
        ax2.grid()
        fig2.tight_layout()


if __name__ == '__main__':
    main()
