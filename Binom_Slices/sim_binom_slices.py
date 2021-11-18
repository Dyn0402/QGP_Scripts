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

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats


def main():
    base_path = '/home/dylan/Research/'
    # div = 120
    divs = [60, 120, 240, 300]
    cents = [8]
    sim_pars = [('amp05', 'spread1'), ('amp1', 'spread1'), ('amp05', 'spread3'), ('amp1', 'spread2'),
                ('amp1', 'spread3'), ('amp15', 'spread1'), ('amp15', 'spread2'), ('amp15', 'spread3'),
                ('amp05', 'spread2')]
    stats = {'mean': getattr(DistStats, 'get_mean'), 'standard deviation': getattr(DistStats, 'get_sd'),
             'skewness': getattr(DistStats, 'get_skewness'), 'kurtosis': getattr(DistStats, 'get_kurtosis')}
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.05), 'skewness': (0, 1.25), 'kurtosis': (-3, 2)}
    stats_plot = ['standard deviation', 'skewness', 'kurtosis']

    # div_stats = {}
    # for div in divs:
    #     div_stats[div] = get_dist_stats(base_path, div, cent, sim_pars, stats)

    df_sim, df_ampt = get_dist_stats(base_path, divs, cents, sim_pars, stats)

    # print(div_stats[120])

    # dist_stats = get_dist_stats(base_path, div, cent, sim_pars, stats)
    #
    plot(df_sim, df_ampt, stats_plot, y_ranges)
    # plot_vs_divs(div_stats, stats_plot, y_ranges)
    print('donzo')


def get_dist_stats(base_path, divs, cents, sim_pars, stats):
    # dist_stats = {}

    sim_energy = 62
    ampt_energy = 62

    ampt_set = 'default/Ampt_rapid05_n1ratios_'

    df_sim = pd.DataFrame()
    df_ampt = pd.DataFrame()

    for div in divs:
        for cent in cents:
            for pars in sim_pars:
                sim_dist = get_sim_dists(base_path + 'Data_Sim/', range(11), div, cent, sim_energy,
                                         exact_keys=['anticlmulti', *pars], contain_keys=['single'])
                sim_dist_mix = get_sim_dists(base_path + 'Data_Sim_Mix/', range(11), div, cent, sim_energy,
                                             exact_keys=['anticlmulti', *pars], contain_keys=['single'])
                stats_dict = get_stats(sim_dist, sim_dist_mix, stats)
                pars_dict = {'div': div, 'cent': cent, 'energy': sim_energy, 'amp': sim_pars[0], 'spread': sim_pars[1]}
                new_df = pd.DataFrame([entry.update(pars_dict) for entry in stats_dict])
                df_sim.append(new_df)
                # dist_stats.update({f'sim_{pars[0]}_{pars[1]}': get_stats(sim_dist, sim_dist_mix, stats)})

            ampt_dist = get_ampt_dist(base_path + 'Data_Ampt/' + ampt_set, div, cent, ampt_energy, range(60))
            ampt_dist_mix = get_ampt_dist(base_path + 'Data_Ampt_Mix/' + ampt_set, div, cent, ampt_energy, range(60))
            stats_dict = get_stats(ampt_dist, ampt_dist_mix, stats)
            pars_dict = {'div': div, 'cent': cent, 'energy': ampt_energy, 'amp': sim_pars[0], 'spread': sim_pars[1]}
            new_df = pd.DataFrame([entry.update(pars_dict) for entry in stats_dict])
            df_ampt.append(new_df)
            # dist_stats.update({'ampt': get_stats(ampt_dist, ampt_dist_mix, stats)})

    return df_sim, df_ampt


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

    # stat_vals = {stat_name: [] for stat_name in stat_names}
    # stat_errs = {stat_name: [] for stat_name in stat_names}
    # stat_syss = {stat_name: [] for stat_name in stat_names}
    data = []
    # total_protons_list = []
    num_sets = len(mix_dists)
    for stat in stat_meases:
        for data_type in data_types:
            for total_protons in sorted(stat_meases[stat][data_type]):
                meases = stat_meases[stat][total_protons]
                if len(meases) == num_sets:  # Only save if all sets are present
                    # if stat == list(stat_names.keys())[0]:  # Only append to list for one stat (first one in key list)
                    #     total_protons_list.append(total_protons)
                    med = np.median(meases)
                    # stat_vals[stat].append(med.val)
                    # stat_errs[stat].append(med.err)
                    # stat_syss[stat].append(np.std(meases).val)
                    data.append({'stat': stat, 'data_type': data_type, 'total_protons': total_protons,
                                 'val': med.val, 'err': med.err, 'sys': np.std(meases).val})
                    # print(sorted(meases))
                    # print(np.median(meases))
                    # print(np.std(meases))

    return data
    # return total_protons_list, stat_vals, stat_errs, stat_syss


def plot(df_sim, df_ampt, stat_names, y_ranges):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(stat)
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])

        for data_set_name, data_set in data_sets.items():
            total_protons, stat_vals, stat_errs, stat_syss = data_set
            vals = np.asarray(stat_vals[stat])
            errs = np.asarray(stat_errs[stat])
            syss = np.asarray(stat_syss[stat])
            if 'ampt' in data_set_name:
                ax.errorbar(total_protons, vals, errs, label=data_set_name, marker='o', ls='', color='blue')
                ax.errorbar(total_protons, vals, syss, marker='', ls='', elinewidth=3, color='blue', alpha=0.3)
            else:
                ax.fill_between(total_protons, vals - errs, vals + errs, label=data_set_name)
        ax.legend()
        fig.tight_layout()

    plt.show()


def plot_vs_divs(div_stats, stat_names, y_ranges):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(stat)
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--')
        ax.set_ylim(y_ranges[stat])
        for data_set_name, data_set in data_sets.items():
            total_protons, stat_vals, stat_errs, stat_syss = data_set
            vals = np.asarray(stat_vals[stat])
            errs = np.asarray(stat_errs[stat])
            syss = np.asarray(stat_syss[stat])
            if 'ampt' in data_set_name:
                ax.errorbar(total_protons, vals, errs, label=data_set_name, marker='o', ls='', color='blue')
                ax.errorbar(total_protons, vals, syss, marker='', ls='', elinewidth=3, color='blue', alpha=0.3)
            else:
                ax.fill_between(total_protons, vals - errs, vals + errs, label=data_set_name)
        ax.legend()
        fig.tight_layout()

    plt.show()


def plot1(dists):
    fig, ax = plt.subplots()
    ax.set_xlabel('Total Protons')
    ax.set_ylabel('Slice Mean')
    ax.grid()
    for dist in dists:
        for raw_mix in dist:
            total_proton_list = []
            # stats = {'mean': [], 'standard deviation': [], 'skewness': [], 'kurtosis': []}
            means = []
            for total_protons in dist:
                stats = DistStats(dist[total_protons])
                total_proton_list.append(total_protons)
                means.append(stats.get_mean().val)
            ax.scatter(total_proton_list, means)
    plt.show()


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


if __name__ == '__main__':
    main()
