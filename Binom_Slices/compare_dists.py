#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 25 2:44 PM 2022
Created in PyCharm
Created as QGP_Scripts/compare_dists.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import warnings
from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from AzimuthBinData import AzimuthBinData as Abd
from Bootstrap_Az_Bin import BootstrapAzBin as Babd

from calc_binom_slices import find_sim_sets, get_name_amp_spread


def main():
    # comp_test()
    # bs_test()
    sim_diff_comp()
    # chi2_test()
    # chi2_test_all()
    # sum_chi2()
    print('donzo')


def sum_chi2():
    weights_path = 'D:/Research/Data_Ampt/default_resample/Ampt_rapid05_resample_norotate_0/62GeV/' \
                   'ratios_divisions_60_centrality_8_local.txt'
    total_proton_dist = Abd(path=weights_path).get_total_particle_dist()
    print(total_proton_dist)
    tp_sum = float(np.sum(np.array(list(total_proton_dist.values()), dtype=np.longlong)))
    total_proton_weights = {total_protons: counts / tp_sum for total_protons, counts in total_proton_dist.items()}
    print(total_proton_weights)
    chi2_indiv_path = 'D:/Research/Results/Azimuth_Analysis/chi2_all_dist_test.csv'
    chi2_sum_out_path = 'D:/Research/Results/Azimuth_Analysis/chi2_sum_dist_test.csv'
    df = pd.read_csv(chi2_indiv_path)
    sums_df = []
    for energy in pd.unique(df['energy']):
        df_energy = df[df['energy'] == energy]
        for data_set in pd.unique(df_energy['data_set']):
            df_data_set = df_energy[df_energy['data_set'] == data_set]
            for sim_set in pd.unique(df_data_set['sim_name']):
                df_sim_set = df_data_set[df_data_set['sim_name'] == sim_set]
                chi2_sum = 0
                for total_protons in pd.unique(df_sim_set['total_protons']):
                    df_tp_set = df_sim_set[df_sim_set['total_protons'] == total_protons]
                    chi2_sum += total_proton_weights[total_protons] * np.sum(df_tp_set['chi2_sum'])
                amp, spread = get_name_amp_spread(sim_set)
                sums_df.append({'data_set': data_set, 'sim_name': sim_set, 'amp': amp, 'spread': spread,
                                'chi2_sum': chi2_sum, 'energy': energy})

    pd.DataFrame(sums_df).to_csv(chi2_sum_out_path, index=False)


def chi2_test_all():
    base_path = 'D:/Research/'
    chi2_out_path = 'D:/Research/Results/Azimuth_Analysis/chi2_all_dist_test2.csv'
    energy = 62
    cent = 8
    divs = [60]  #, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    threads = 10

    data_sets = [
        (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
         energy, cent, divs, 'Data_Ampt', 'Data_Ampt_Mix'),
    ]

    df_sim_sets = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'])

    sim_sets = []

    for amp in np.unique(df_sim_sets['amp']):
        if amp not in ['0']:
            continue
        df_amp = df_sim_sets[df_sim_sets['amp'] == amp]
        for spread in np.unique(df_amp['spread']):
            if spread not in ['001']:
                continue
            sim_sets.append((base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                             f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                             energy, cent, divs, 'Data_Sim', 'Data_Sim_Mix'))
    chi2_df = []
    for data_set in data_sets:
        div_diffs_data = get_set_all(data_set)
        jobs = [(sim_set, div_diffs_data, divs, energy) for sim_set in sim_sets]
        with Pool(threads) as pool:
            for df_subset in tqdm.tqdm(pool.istarmap(run_sim_set, jobs), total=len(jobs)):
                chi2_df.extend(df_subset)

    # chi2_df = []
    # for data_set in data_sets:
    #     div_diffs_data = get_set_all(data_set)
    #     for sim_set in sim_sets:
    #         print(sim_set[1])
    #         div_diffs_sim = get_set_all(sim_set)
    #         for div in divs:
    #             for total_protons, (diff_def_data, diff_sds_data) in div_diffs_data[div].items():
    #                 diff_def_sim, diff_sds_sim = div_diffs_sim[div][total_protons]
    #
    #                 data_sim_diff = diff_def_data - diff_def_sim
    #                 data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)
    #                 chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
    #                 chi2_sum = np.sum(chi2)
    #                 amp, spread = get_name_amp_spread(sim_set[1], as_type='string')
    #                 amp_f, spread_f = get_name_amp_spread(sim_set[1], as_type='float')
    #                 chi2_df.append({'sim_name': f'sim_amp{amp}_spread{spread}', 'chi2_sum': chi2_sum, 'amp': amp_f,
    #                                 'spread': spread_f, 'data_set': 'ampt_resample_def', 'energy': energy,
    #                                 'total_protons': total_protons, 'divs': div})

    chi2_df = pd.DataFrame(chi2_df)
    chi2_df.to_csv(chi2_out_path, index=False)


def run_sim_set(sim_set, div_diffs_data, divs, energy):
    # print(sim_set[1])
    div_diffs_sim = get_set_all(sim_set)
    chi2_df = []
    for div in divs:
        for total_protons, (diff_def_data, diff_sds_data) in div_diffs_data[div].items():
            diff_def_sim, diff_sds_sim = div_diffs_sim[div][total_protons]

            data_sim_diff = diff_def_data - diff_def_sim
            data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)
            chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
            chi2_sum = np.nansum(chi2)
            num_points = np.sum(~np.isnan(chi2))
            amp, spread = get_name_amp_spread(sim_set[1], as_type='string')
            amp_f, spread_f = get_name_amp_spread(sim_set[1], as_type='float')
            chi2_df.append({'sim_name': f'sim_amp{amp}_spread{spread}', 'chi2_sum': chi2_sum, 'num_points': num_points,
                            'amp': amp_f, 'spread': spread_f, 'data_set': 'ampt_resample_def', 'energy': energy,
                            'total_protons': total_protons, 'divs': div})

    return chi2_df


def chi2_test():
    base_path = 'D:/Research/'
    chi2_out_path = 'D:/Research/Results/Azimuth_Analysis/chi2_dist_test.csv'
    energy = 62
    cent = 8
    div = 60

    total_protons = 4

    data_sets = [
        (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
         energy, cent, div, total_protons, 'Data_Ampt', 'Data_Ampt_Mix'),
    ]

    df_sim_sets = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'])

    sim_sets = []

    for amp in np.unique(df_sim_sets['amp']):
        df_amp = df_sim_sets[df_sim_sets['amp'] == amp]
        for spread in np.unique(df_amp['spread']):
            sim_sets.append((base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                             f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                             energy, cent, div, total_protons, 'Data_Sim', 'Data_Sim_Mix'))

    chi2_df = []
    for data_set in data_sets:
        diff_def_data, diff_sds_data = get_set(data_set)
        for sim_set in sim_sets:
            print(sim_set[1])
            diff_def_sim, diff_sds_sim = get_set(sim_set)
            # print('diff_def_data: ', diff_def_data)
            # print('diff_def_sim: ', diff_def_sim)
            # print('diff_sds_data: ', diff_sds_data)
            # print('diff_sds_sim: ', diff_sds_sim)
            data_sim_diff = diff_def_data - diff_def_sim
            # print('data_sim_diff: ', data_sim_diff)
            data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)
            # print('data_sim_diff_err2: ', data_sim_diff_err2)
            chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
            # print('chi2: ', chi2)
            chi2_sum = np.sum(chi2)
            # print('chi2_sum: ', chi2_sum)
            # chi2 = np.power(diff_def_data - diff_def_sim, 2) / (np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2))
            amp, spread = get_name_amp_spread(sim_set[1], as_type='string')
            amp_f, spread_f = get_name_amp_spread(sim_set[1], as_type='float')
            chi2_df.append({'sim_name': f'sim_amp{amp}_spread{spread}', 'chi2_sum': chi2_sum, 'amp': amp_f,
                            'spread': spread_f, 'data_set': 'ampt_resample_def', 'energy': energy})

    chi2_df = pd.DataFrame(chi2_df)
    chi2_df.to_csv(chi2_out_path, index=False)


def sim_diff_comp():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 60

    total_protons = 21

    sim_amp, sim_spread = '015', '1'
    data_sets = [
        (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
         energy, cent, div, total_protons, 'Data_Ampt', 'Data_Ampt_Mix'),
        (base_path, f'flat80_anticlmulti_spread{sim_spread}_amp{sim_amp}_resample',
         f'Sim_spread{sim_spread}_amp{sim_amp}_flat80_anticlmulti_norotate_resample_0',
         energy, cent, div, total_protons, 'Data_Sim', 'Data_Sim_Mix'),
    ]

    for data_set in data_sets:
        base_path, set_group, set_name, energy, cent, div, total_protons, raw_folder, mix_folder = data_set
        file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
        path_sufx = f'{set_group}/{set_name}/{energy}GeV/{file_name}'
        raw_tp_dist = get_norm_dist(f'{base_path}{raw_folder}/{path_sufx}', total_protons)
        mix_tp_dist = get_norm_dist(f'{base_path}{mix_folder}/{path_sufx}', total_protons)

        name = set_name[:set_name.find('_', set_name.find('_', set_name.find('_') + 1) + 1)]

        fig_comp, ax_comp = plt.subplots()
        ax_comp.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist), marker='o', alpha=0.7, label='Raw')
        ax_comp.plot(range(mix_tp_dist.size), mix_tp_dist / np.sum(mix_tp_dist), marker='o', alpha=0.7, label='Mix')
        ax_comp.set_title(f'{name} {total_protons} Protons')
        ax_comp.axhline(0, ls='--', color='black', zorder=0)
        ax_comp.set_xlabel('Protons in Bin')
        ax_comp.legend()

    fig_diff, ax_diff = plt.subplots()
    x = np.arange(total_protons + 1)

    for data_set in data_sets:
        set_name = data_set[2]
        name = set_name[:set_name.find('_', set_name.find('_', set_name.find('_') + 1) + 1)]
        diff_def, diff_sds = get_set(data_set)
        ax_diff.fill_between(x, diff_def + diff_sds, diff_def - diff_sds, alpha=0.5, label=name)
        ax_diff.plot(x, diff_def)
    ax_diff.axhline(0, ls='--', color='black')
    ax_diff.set_xlabel('Protons in Bin')
    ax_diff.set_ylabel('Raw - Mix / Avg')
    ax_diff.text(0.05, 0.3, f'{total_protons} Protons', fontsize='large', transform=ax_diff.transAxes)
    ax_diff.legend()
    fig_diff.tight_layout()

    plt.show()


def get_set(set_pars):
    base_path, set_group, set_name, energy, cent, div, total_protons, raw_folder, mix_folder = set_pars
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy}GeV/{file_name}'
    diff_def, diff_sds = get_diff(f'{base_path}{raw_folder}/{path_sufx}', f'{base_path}{mix_folder}/{path_sufx}',
                                  total_protons)

    return diff_def, diff_sds


def get_set_all(set_pars):
    base_path, set_group, set_name, energy, cent, divs, raw_folder, mix_folder = set_pars
    div_data = {}
    for div in divs:
        file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
        path_sufx = f'{set_group}/{set_name}/{energy}GeV/{file_name}'
        total_proton_data = get_diff_all(f'{base_path}{raw_folder}/{path_sufx}', f'{base_path}{mix_folder}/{path_sufx}')
        div_data.update({div: total_proton_data})

    return div_data


def bs_test():
    base_path = 'D:/Research/'
    set_group = 'default_resample'
    set_name = 'Ampt_rapid05_resample_norotate_0'
    energy = 62
    cent = 8
    div = 120

    total_protons = 20

    raw_folder = 'Data_Ampt'
    mix_folder = 'Data_Ampt_Mix'

    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy}GeV/{file_name}'

    diff_def, diff_sds = get_diff(f'{base_path}{raw_folder}/{path_sufx}', f'{base_path}{mix_folder}/{path_sufx}',
                                  total_protons)

    fig_diff, ax_diff = plt.subplots()
    x = np.arange(total_protons + 1)
    ax_diff.fill_between(x, diff_def + diff_sds, diff_def - diff_sds, alpha=0.5)
    ax_diff.plot(x, diff_def)
    # ax_diff.plot(range(len(diff)), diff, marker='o')
    ax_diff.axhline(0, ls='--', color='black')
    ax_diff.set_xlabel('Protons in Bin')
    ax_diff.set_ylabel('Raw - Mix / Avg')

    plt.show()


def get_diff(raw_path, mix_path, total_protons):
    raw_def, raw_bss = get_norm_dists(raw_path, total_protons)
    mix_def, mix_bss = get_norm_dists(mix_path, total_protons)

    # with warnings.catch_warnings():  # Ignore divide by zero warnings when raw=mix=0. From 0 padding
    #     warnings.simplefilter('ignore')
    # ---------------------------------------------------------------------------------------
    diff_def = (raw_def - mix_def)  # / (2 * (raw_def + mix_def))
    # ---------------------------------------------------------------------------------------
    raw_bs_tile, mix_bs_repeat = np.tile(raw_bss, (len(raw_bss), 1)), np.repeat(mix_bss, len(mix_bss), axis=0)
    # ---------------------------------------------------------------------------------------
    bs_diffs = (raw_bs_tile - mix_bs_repeat)  # / (2 * (raw_bs_tile + mix_bs_repeat))
    # ---------------------------------------------------------------------------------------
    diff_sds = np.std(bs_diffs, axis=0)
    diff_def = np.where((diff_sds == 0) | (np.isnan(diff_sds)), float('nan'), diff_def)

    return diff_def, diff_sds


def get_total_proton_list(raw_def, mix_def, raw_bss, mix_bss, min_bs=100):
    raw_def = set(raw_def.keys())
    mix_def = set(mix_def.keys())
    raw_bss_num = {total_protons for total_protons, bss in raw_bss.items() if len(bss) >= min_bs}
    mix_bss_num = {total_protons for total_protons, bss in mix_bss.items() if len(bss) >= min_bs}

    return raw_def.intersection(mix_def, raw_bss_num, mix_bss_num)


def get_diff_all(raw_path, mix_path):
    raw_def_tp, raw_bss_tp = get_norm_dists_all(raw_path)
    mix_def_tp, mix_bss_tp = get_norm_dists_all(mix_path)

    total_protons = get_total_proton_list(raw_def_tp, mix_def_tp, raw_bss_tp, mix_bss_tp)

    total_proton_data = {}
    # with warnings.catch_warnings():  # Ignore divide by zero warnings when raw=mix=0. From 0 padding
    #     warnings.simplefilter('ignore')
    for total_proton in total_protons:
        raw_def, mix_def = raw_def_tp[total_proton], mix_def_tp[total_proton]
        raw_bss, mix_bss = raw_bss_tp[total_proton], mix_bss_tp[total_proton]

        # ---------------------------------------------------------------------------------------
        diff_def = (raw_def - mix_def)  # / (2 * (raw_def + mix_def))
        # ---------------------------------------------------------------------------------------
        raw_bs_tile, mix_bs_repeat = np.tile(raw_bss, (len(mix_bss), 1)), np.repeat(mix_bss, len(raw_bss), axis=0)
        # ---------------------------------------------------------------------------------------
        bs_diffs = (raw_bs_tile - mix_bs_repeat)  # / (2 * (raw_bs_tile + mix_bs_repeat))
        # ---------------------------------------------------------------------------------------
        diff_sd = np.std(bs_diffs, axis=0)
        diff_def = np.where((diff_sd == 0) | (np.isnan(diff_sd)), float('nan'), diff_def)
        total_proton_data.update({total_proton: (diff_def, diff_sd)})

    return total_proton_data


def comp_test():
    base_path = 'D:/Research/'
    set_group = 'default_resample'
    set_name = 'Ampt_rapid05_resample_norotate_0'
    energy = 62
    cent = 8
    div = 120

    total_protons = 20

    raw_folder = 'Data_Ampt'
    mix_folder = 'Data_Ampt_Mix'

    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy}GeV/{file_name}'

    raw_tp_dist = get_norm_dist(f'{base_path}{raw_folder}/{path_sufx}', total_protons)
    mix_tp_dist = get_norm_dist(f'{base_path}{mix_folder}/{path_sufx}', total_protons)
    # max_len = max(raw_dist.size, mix_dist.size)
    # raw_tp_dist = np.pad(raw_dist, (0, max_len - raw_dist.size))
    # mix_tp_dist = np.pad(mix_dist, (0, max_len - mix_dist.size))

    fig_comp, ax_comp = plt.subplots()
    ax_comp.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist), marker='o', alpha=0.7, label='Raw')
    ax_comp.plot(range(mix_tp_dist.size), mix_tp_dist / np.sum(mix_tp_dist), marker='o', alpha=0.7, label='Mix')
    ax_comp.axhline(0, ls='--', color='black', zorder=0)
    ax_comp.set_xlabel('Protons in Bin')
    ax_comp.legend()

    with warnings.catch_warnings():  # Ignore divide by zero warnings when raw=mix=0. From 0 padding
        warnings.simplefilter('ignore')
        diff = (raw_tp_dist - mix_tp_dist) / (2 * (raw_tp_dist + mix_tp_dist))
    fig_diff, ax_diff = plt.subplots()
    ax_diff.plot(range(len(diff)), diff, marker='o')
    ax_diff.axhline(0, ls='--', color='black')
    ax_diff.set_xlabel('Protons in Bin')
    ax_diff.set_ylabel('Raw - Mix / Avg')

    plt.show()


def get_norm_dist(path, total_protons):
    dist = np.array(Abd(path=path).get_dist()[total_protons])
    return np.pad(dist / np.sum(dist), (0, total_protons + 1 - dist.size))


def get_norm_dists(path, total_protons):
    abin_data = Babd(path=path)
    def_dist = np.array(abin_data.get_dist()[total_protons])
    def_dist_norm = np.pad(def_dist / np.sum(def_dist), (0, total_protons + 1 - def_dist.size))

    bs_dists_norm = []
    for bs_dist in abin_data.get_dist_bs():
        bs_dist = np.array(bs_dist[total_protons])
        bs_dist_norm = np.pad(bs_dist / np.sum(bs_dist), (0, total_protons + 1 - bs_dist.size))
        bs_dists_norm.append(bs_dist_norm)

    return def_dist_norm, bs_dists_norm


def get_norm_dists_all(path):
    abin_data = Babd(path=path)
    def_dist = abin_data.get_dist()
    def_dists_norm = {}
    for total_protons, dist in def_dist.items():
        dist = np.array(dist)
        def_dists_norm.update({total_protons: np.pad(dist / np.sum(dist), (0, total_protons + 1 - dist.size))})

    bs_dists_norm = {}
    for bs_dist in abin_data.get_dist_bs():
        for total_protons, dist in bs_dist.items():
            dist = np.array(dist)
            bs_dist_norm = np.pad(dist / np.sum(dist), (0, total_protons + 1 - dist.size))
            if total_protons in bs_dists_norm:
                bs_dists_norm[total_protons].append(bs_dist_norm)
            else:
                bs_dists_norm.update({total_protons: [bs_dist_norm]})

    return def_dists_norm, bs_dists_norm


if __name__ == '__main__':
    main()
