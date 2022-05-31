#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 25 2:44 PM 2022
Created in PyCharm
Created as QGP_Scripts/compare_dists.py

@author: Dylan Neff, Dylan
"""
import os.path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas.errors import EmptyDataError

import warnings
from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from AzimuthBinData import AzimuthBinData as Abd
from Bootstrap_Az_Bin import BootstrapAzBin as Babd

from calc_binom_slices import find_sim_sets, get_name_amp_spread


# import pprofile


def main():
    # comp_test()
    # bs_test()
    # sim_diff_comp()
    # chi2_test()
    # bes_info = ('default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0', 'Data', 'Data_Mix',
    #             'bes_resample')
    # chi2_test_all(bes_info, 'F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv')
    # sum_chi2('F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv',
    #          'F:/Research/Results/Azimuth_Analysis/bes_chi2_sum_dist_bs.csv',
    #          'F:/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
    #          )
    # ampt_old_info = ('default_resample', 'Ampt_rapid05_resample_norotate_0', 'Data_Ampt_Old', 'Data_Ampt_Old_Mix',
    #                  'ampt_resample_def')
    # chi2_test_all(ampt_old_info, 'F:/Research/Results/Azimuth_Analysis/ampt_old_chi2_all_dist_bs.csv')
    # sum_chi2('F:/Research/Results/Azimuth_Analysis/ampt_old_chi2_all_dist_bs.csv',
    #          'F:/Research/Results/Azimuth_Analysis/ampt_old_chi2_sum_dist_bs.csv',
    #          'F:/Research/Data_Ampt/default_resample/Ampt_rapid05_resample_norotate_0/')
    # ampt_types = ['new']
    # effs = ['3', '2', '1']
    # for ampt_type in ampt_types:
    #     type_flag = '_Old' if ampt_type == 'old' else ''
    #     for eff in effs:
    #         ampt_set_info = (f'Eff{eff}_resample', f'Ampt{type_flag}_rapid05_resample_norotate_Efficiency{eff}_0',
    #                          f'Data_Ampt{type_flag}', f'Data_Ampt{type_flag}_Mix',
    #                          f'ampt_{ampt_type}_resample_eff{eff}')
    #         chi2_test_all(ampt_set_info,
    #                       f'F:/Research/Results/Azimuth_Analysis/ampt_{ampt_type}_eff{eff}_chi2_all_dist_bs.csv')
    #         sum_chi2(f'F:/Research/Results/Azimuth_Analysis/ampt_{ampt_type}_eff{eff}_chi2_all_dist_bs.csv',
    #                  f'F:/Research/Results/Azimuth_Analysis/ampt_{ampt_type}_eff{eff}_chi2_sum_dist_bs.csv',
    #                  f'F:/Research/Data_Ampt/Eff{eff}_resample/Ampt_rapid05_resample_norotate_Efficiency{eff}_0/')
    spread_amps = (('08', '02'), ('08', '05'), ('05', '02'), ('1', '02'))
    # spread, amp = '08', '02'
    for spread, amp in spread_amps:
        sim_set_info = (f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                        f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                        'Data_Sim', 'Data_Sim_Mix', f'sim_sim_comp_s{spread}_a{amp}_resample')
        chi2_test_all(sim_set_info,
                      f'F:/Research/Results/Azimuth_Analysis/sim_sim_comp_s{spread}_a{amp}_chi2_all_dist_bs.csv',
                      energies=[62])
        sum_chi2(f'F:/Research/Results/Azimuth_Analysis/sim_sim_comp_s{spread}_a{amp}_chi2_all_dist_bs.csv',
                 f'F:/Research/Results/Azimuth_Analysis/sim_sim_comp_s{spread}_a{amp}_chi2_sum_dist_bs.csv',
                 f'F:/Research/Data_Sim/flat80_anticlmulti_spread{spread}_amp{amp}_resample/'
                 f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0/')
    print('donzo')


def sum_chi2(chi2_indiv_path='F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv',
             chi2_sum_out_path='F:/Research/Results/Azimuth_Analysis/bes_chi2_sum_dist_bs.csv',
             weights_path_pre='F:/Research/Data/default_resample/'
                              'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'):
    # weights_path_pre = 'F:/Research/Data_Ampt/default_resample/Ampt_rapid05_resample_norotate_0/'
    # weights_path_pre = '/home/dylan/Research/Data_Ampt/default_resample/Ampt_rapid05_resample_norotate_0/'
    # weights_path_pre = 'F:/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
    cent = 8
    div_weight = 60  # Proton distribution for all divs should be the same
    # chi2_indiv_path = 'F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv'
    # chi2_indiv_path = '/home/dylan/Research/Results/Azimuth_Analysis/chi2_all_dist_ampt_new.csv'
    # chi2_sum_out_path = 'F:/Research/Results/Azimuth_Analysis/bes_chi2_sum_dist_bs.csv'
    # chi2_sum_out_path = '/home/dylan/Research/Results/Azimuth_Analysis/chi2_sum_dist_ampt_new.csv'
    threads = 12

    print(f'Reading input csv {chi2_indiv_path}')
    df = pd.read_csv(chi2_indiv_path)

    if os.path.exists(chi2_sum_out_path):
        res = input(f'Output file exists {chi2_sum_out_path}\n'
                    f'Enter "a" to append to dataframe, "o" to overwrite, anything else to'
                    f' quit:\n')
        if res.lower() == 'a':
            pass
        elif res.lower() == 'o':
            os.remove(chi2_sum_out_path)

    for energy in pd.unique(df['energy']):
        df_energy = df[df['energy'] == energy]
        weights_path = f'{weights_path_pre}{energy}GeV/ratios_divisions_{div_weight}_centrality_{cent}_local.txt'
        total_proton_dist = Abd(path=weights_path).get_total_particle_dist()
        tp_sum = float(np.sum(np.array(list(total_proton_dist.values()), dtype=np.longlong)))
        total_proton_weights = {total_protons: counts / tp_sum for total_protons, counts in
                                total_proton_dist.items()}
        for data_set in pd.unique(df_energy['data_set']):
            print(f'{energy}GeV {data_set}')
            sums_df = []
            jobs = []
            df_data_set = df_energy[df_energy['data_set'] == data_set]
            for sim_set in pd.unique(df_data_set['sim_name']):
                df_sim_set = df_data_set[df_data_set['sim_name'] == sim_set]
                amp, spread = get_name_amp_spread(sim_set)
                df_entry = {'data_set': data_set, 'sim_name': sim_set, 'energy': energy, 'amp': amp, 'spread': spread}
                jobs.append((df_sim_set, total_proton_weights, df_entry))

            with Pool(threads) as pool:
                for df_entry in tqdm.tqdm(pool.istarmap(sum_chi2_set, jobs), total=len(jobs)):
                    sums_df.append(df_entry)

            sums_df = pd.DataFrame(sums_df)
            try:
                df_old = pd.read_csv(chi2_sum_out_path)
                sums_df = sums_df.append(df_old, ignore_index=True)
            except (FileNotFoundError, EmptyDataError) as e:
                pass
            print('Outputting to csv')
            sums_df.to_csv(chi2_sum_out_path, index=False)


def sum_chi2_set(df_sim_set, total_proton_weights, df_entry):
    chi2_sum, chi2_sum_num = 0, 0
    chi2_sum_bss, chi2_sum_num_bss = {}, {}
    for div in pd.unique((df_sim_set['divs'])):
        df_divs = df_sim_set[df_sim_set['divs'] == div]
        for total_protons in pd.unique(df_divs['total_protons']):
            df_tp_set = df_divs[df_divs['total_protons'] == total_protons]
            # chi2_sum += total_proton_weights[total_protons] * \
            #             np.sum(df_tp_set['chi2_sum'] / df_tp_set['num_points'])  # div weights 1
            # div weights are all 1
            if total_protons in total_proton_weights:
                tp_weight = total_proton_weights[total_protons]
            else:
                tp_weight = min(total_proton_weights.values())  # If not in original dist, give min non-zero weight
            chi2_sum += tp_weight * np.sum(df_tp_set['chi2_avg_def'])
            chi2_sum_num += 1
            bs_cols = [col for col in df_tp_set if 'chi2_avg_bs' in col]
            for col in bs_cols:
                if col not in chi2_sum_bss:
                    chi2_sum_bss.update({col: 0})
                    chi2_sum_num_bss.update({col: 0})
                chi2_sum_bss[col] += tp_weight * np.sum(df_tp_set[col])
                chi2_sum_num_bss[col] += 1
    df_entry.update({'chi2_avg': chi2_sum / chi2_sum_num})
    for col in chi2_sum_bss:
        df_entry.update({col: chi2_sum_bss[col] / chi2_sum_num_bss[col]})

    return df_entry


def chi2_test_all(data_set_info, chi2_out_path='F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv',
                  energies=None):
    # base_path = 'F:/Research/'
    base_path = 'F:/Research/'
    # base_path = '/home/dylan/Research/'
    # chi2_out_path = 'F:/Research/Results/Azimuth_Analysis/chi2_all_dist_ampt_new_bstest3.csv'
    # chi2_out_path = 'F:/Research/Results/Azimuth_Analysis/bes_chi2_all_dist_bs.csv'
    df_append = True
    # chi2_out_path = '/home/dylan/Research/Results/Azimuth_Analysis/chi2_all_dist_ampt_new.csv'
    # energy = 62
    energies = [7, 11, 19, 27, 39, 62] if energies is None else energies
    # energies = [11]
    sim_energy = 62
    cent = 8
    divs = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    # divs = [60]
    threads = 12

    if os.path.exists(chi2_out_path):
        res = input(f'File exists {chi2_out_path}\n'
                    f'Enter "a" to append, "o" to overwrite, anything else to exit:\n').lower()
        if res == 'a':
            df_append = True
        elif res == 'o':
            os.remove(chi2_out_path)
        else:
            print("Don't want to overwrite, exiting.")
            return

    print('Get sim set list')
    df_sim_sets = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'])
    sim_sets = []
    for amp in pd.unique(df_sim_sets['amp']):
        # if amp not in ['015']:
        #     continue
        df_amp = df_sim_sets[df_sim_sets['amp'] == amp]
        for spread in pd.unique(df_amp['spread']):
            # if spread not in ['1']:
            #     continue
            sim_sets.append((base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                             f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                             sim_energy, cent, divs, 'Data_Sim', 'Data_Sim_Mix'))

    # data_sets = [
    #     # (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
    #     #  energy, cent, divs, 'Data_Ampt', 'Data_Ampt_Mix'),
    # ]
    for energy in energies:
        print(f'Starting {energy}GeV')
        # data_sets.extend([(base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
        #                    energy, cent, divs, 'Data_Ampt', 'Data_Ampt_Mix', 'ampt_resample_def')
        #                   for energy in energies])
        # data_sets = [(base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
        #               energy, cent, divs, 'Data_Ampt', 'Data_Ampt_Mix', 'ampt_resample_def')]
        data_sets = [(base_path, data_set_info[0], data_set_info[1],
                      energy, cent, divs, data_set_info[2], data_set_info[3], data_set_info[4])]
        # data_sets.extend([(base_path, 'default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0',
        #                    energy, cent, divs, 'Data', 'Data_Mix', 'bes_resample') for energy in energies])

        chi2_df = []
        for data_set in data_sets:
            data_set, data_name = data_set[:-1], data_set[-1]
            data_energy = data_set[3]
            print(f'Get {data_name}, {data_energy}GeV distribution')
            div_diffs_data = get_set_all(data_set)
            # return
            jobs = [(sim_set, div_diffs_data, divs, data_energy, data_name) for sim_set in sim_sets]
            print(f'Starting {data_name}, {data_energy}GeV comparisons')
            with Pool(threads) as pool:
                for df_subset in tqdm.tqdm(pool.istarmap(run_sim_set, jobs), total=len(jobs)):
                    chi2_df.extend(df_subset)

        print('Making dataframe')
        chi2_df = pd.DataFrame(chi2_df)
        print('Saving dataframe')
        if df_append:
            try:
                df_old = pd.read_csv(chi2_out_path)
                chi2_df = chi2_df.append(df_old, ignore_index=True)
            except (FileNotFoundError, EmptyDataError) as e:
                pass
        chi2_df.to_csv(chi2_out_path, index=False)


def run_sim_set(sim_set, div_diffs_data, divs, energy, data_name):
    div_diffs_sim = get_set_all(sim_set)
    chi2_df = []
    for div in divs:
        for total_protons, (diff_def_data, diff_sds_data, diffs_bs_data) in div_diffs_data[div].items():
            diff_def_sim, diff_sds_sim, diffs_bs_sim = div_diffs_sim[div][total_protons]

            amp, spread = get_name_amp_spread(sim_set[1], as_type='string')
            amp_f, spread_f = get_name_amp_spread(sim_set[1], as_type='float')
            data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)  # Maybe do better?
            data_sim_diff = diff_def_data - diff_def_sim
            chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
            chi2_sum = np.nansum(chi2)
            num_points = np.sum(~np.isnan(chi2))
            chi2_avg = chi2_sum / num_points
            df_entry = {'sim_name': f'sim_amp{amp}_spread{spread}', 'data_set': data_name, 'energy': energy,
                        'divs': div, 'total_protons': total_protons, 'amp': amp_f, 'spread': spread_f,
                        'chi2_avg_def': chi2_avg}
            bs_range = min(len(diffs_bs_data), len(diffs_bs_sim))
            data_sim_bs_diffs = diffs_bs_data[:bs_range] - diffs_bs_sim[:bs_range]
            chi2_bss = np.power(data_sim_bs_diffs, 2) / data_sim_diff_err2
            chi2_sum_bss = np.nansum(chi2_bss, axis=1)
            num_points_bss = np.sum(~np.isnan(chi2_bss), axis=1)
            chi2_avg_bss = chi2_sum_bss / num_points_bss
            df_entry.update({f'chi2_avg_bs{bs_i}': chi2_avg_bs for bs_i, chi2_avg_bs in enumerate(chi2_avg_bss)})
            chi2_df.append(df_entry)

    return chi2_df


def chi2_test():
    base_path = 'F:/Research/'
    chi2_out_path = 'F:/Research/Results/Azimuth_Analysis/chi2_dist_test.csv'
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
    base_path = 'F:/Research/'
    energy = 62
    cent = 8
    divs = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    plot_raw = False

    total_protons = [21]

    sim_pars = [('0', '08'), ('0', '1'), ('0', '09'),
                ('005', '09')]  # ('0125', '08'), ('45', '35') , ('4', '1'), ('6', '45')

    for total_proton in total_protons:
        chi2s = {}
        fig_diff_divs, ax_diff_divs = plt.subplots(4, 3, sharex=True, figsize=(14, 8))
        ax_diff_divs = ax_diff_divs.flat
        div_index = 0
        for div in divs:
            data_sets = [
                # (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
                #  energy, cent, div, total_proton, 'Data_Ampt', 'Data_Ampt_Mix', 'ampt_new'),
                # (base_path, 'default_resample', 'Ampt_rapid05_resample_norotate_0',
                #  energy, cent, div, total_proton, 'Data_Ampt_Old', 'Data_Ampt_Old_Mix', 'ampt_old'),
                (base_path, 'default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0',
                 energy, cent, div, total_proton, 'Data', 'Data_Mix', 'bes'),
            ]

            sim_sets = []
            for sim_par in sim_pars:
                sim_sets.append((base_path, f'flat80_anticlmulti_spread{sim_par[1]}_amp{sim_par[0]}_resample',
                                 f'Sim_spread{sim_par[1]}_amp{sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
                                 62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
                                 f'Sim_spread{sim_par[1]}_amp{sim_par[0]}'))

            for data_set in sim_sets + data_sets:
                base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, \
                name = data_set
                file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
                path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
                raw_tp_dist = get_norm_dist(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
                mix_tp_dist = get_norm_dist(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

                # name = set_name[:set_name.find('_', set_name.find('_', set_name.find('_') + 1) + 1)]

                if plot_raw:
                    fig_comp, ax_comp = plt.subplots()
                    ax_comp.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist), marker='o', alpha=0.7,
                                 label='Raw')
                    ax_comp.plot(range(mix_tp_dist.size), mix_tp_dist / np.sum(mix_tp_dist), marker='o', alpha=0.7,
                                 label='Mix')
                    ax_comp.set_title(f'{name} {total_proton} Protons, {div} divs')
                    ax_comp.axhline(0, ls='--', color='black', zorder=0)
                    ax_comp.set_xlabel('Protons in Bin')
                    ax_comp.legend()
                    fig_comp.canvas.manager.set_window_title(f'{name} {total_proton} Protons, {div} divs')

            fig_diff, ax_diff = plt.subplots()
            x = np.arange(total_proton + 1)

            data_sets_diffs = {}
            for data_set in sim_sets + data_sets:
                # set_name = data_set[2]
                # name = set_name[:set_name.find('_', set_name.find('_', set_name.find('_') + 1) + 1)]
                data_set, name = data_set[:-1], data_set[-1]
                diff_def, diff_sds = get_set(data_set)
                data_sets_diffs.update({name: {'diff_def': diff_def, 'diff_sds': diff_sds}})
                ax_diff.fill_between(x, diff_def + diff_sds, diff_def - diff_sds, alpha=0.5, label=name)
                ax_diff.plot(x, diff_def)
                ax_diff_divs[div_index].fill_between(x, diff_def + diff_sds, diff_def - diff_sds, alpha=0.5, label=name)
                ax_diff_divs[div_index].plot(x, diff_def)
                ax_diff_divs[div_index].set_ylim(-0.005, 0.0082)
            ax_diff.axhline(0, ls='--', color='black')
            ax_diff.set_xlabel('Protons in Bin')
            ax_diff.set_ylabel('Raw - Mix')
            ax_diff.text(0.75, 0.1, f'{energy}GeV\n{total_proton} Protons\n{div}° Divisions', fontsize='large',
                         transform=ax_diff.transAxes)
            ax_diff_divs[div_index].axhline(0, ls='--', color='black')
            if div_index >= 9:
                ax_diff_divs[div_index].set_xlabel('Protons in Bin')
            ax_diff_divs[div_index].text(total_proton * 0.4, -0.004, f'{div}° Divisions', fontsize='large')
            ax_diff.legend()
            fig_diff.tight_layout()
            fig_diff.canvas.manager.set_window_title(f'{energy}GeV, {total_proton} Protons, {div}° Divisions')
            div_index += 1

            # chi2s.append([])
            for data_set in data_sets:
                for sim_set in sim_sets:
                    if data_set[-1] not in chi2s:
                        chi2s.update({data_set[-1]: {}})
                    if sim_set[-1] not in chi2s[data_set[-1]]:
                        chi2s[data_set[-1]].update({sim_set[-1]: []})
                    chi2 = get_chi2(data_sets_diffs[data_set[-1]]['diff_def'],
                                    data_sets_diffs[data_set[-1]]['diff_sds'],
                                    data_sets_diffs[sim_set[-1]]['diff_def'], data_sets_diffs[sim_set[-1]]['diff_sds'])
                    chi2s[data_set[-1]][sim_set[-1]].append(chi2)

            # chi2_ampt = get_chi2(data_sets_diffs[0]['diff_def'], data_sets_diffs[0]['diff_sds'],
            #                      data_sets_diffs[2]['diff_def'], data_sets_diffs[2]['diff_sds'])
            # chi2_bes = get_chi2(data_sets_diffs[1]['diff_def'], data_sets_diffs[1]['diff_sds'],
            #                     data_sets_diffs[2]['diff_def'], data_sets_diffs[2]['diff_sds'])
            #
            # print(f'energy: {energy}, cent: {cent}, div: {div}, total_protons: {total_proton}')
            # print(f'BES chi2: {chi2_bes}')
            # print(f'AMPT chi2: {chi2_ampt}')
        for data_set in data_sets:
            fig_chi_div, ax_chi_div = plt.subplots()
            for sim_set in sim_sets:
                sim_lab = sim_set[-1].replace('Sim_', '').replace('_flat80_anticlmulti_norotate_resample_0', '')
                ax_chi_div.plot(divs, chi2s[data_set[-1]][sim_set[-1]], marker='o', label=sim_lab)
            ax_chi_div.set_xlabel('Division Width')
            ax_chi_div.set_ylabel('Chi_Square')
            ax_chi_div.set_ylim(bottom=0)
            ax_chi_div.set_title(data_set[-1])
            ax_chi_div.legend()
            fig_chi_div.tight_layout()
            fig_chi_div.canvas.manager.set_window_title(data_set[-1])
        for data_set in sim_sets + data_sets:
            ax_diff_divs[-1].fill_between([], [], [], alpha=0.5, label=data_set[-1])
        ax_diff_divs[-1].legend()
        fig_diff_divs.canvas.manager.set_window_title(f'{energy}GeV, {total_proton} Protons, All Divisions')
        fig_diff_divs.tight_layout()
        fig_diff_divs.subplots_adjust(wspace=0.18, hspace=0.05)

    plt.show()


def get_chi2(diff_def_data, diff_sds_data, diff_def_sim, diff_sds_sim):
    data_sim_diff = diff_def_data - diff_def_sim
    data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)
    chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
    chi2_sum = np.nansum(chi2)

    return chi2_sum


def get_chi2_print(diff_def_data, diff_sds_data, diff_def_sim, diff_sds_sim):
    # print(f'Data_def: {diff_def_data}')
    # print(f'Sim_def: {diff_def_sim}')
    # print(f'Data_sds: {diff_sds_data}')
    # print(f'Sim_sds: {diff_sds_sim}')
    data_sim_diff = diff_def_data - diff_def_sim
    # print(f'Data - Sim: {data_sim_diff}')
    data_sim_diff_err2 = np.power(diff_sds_data, 2) + np.power(diff_sds_sim, 2)
    # print(f'Data - Sim Err: {data_sim_diff_err2}')
    chi2 = np.power(data_sim_diff, 2) / data_sim_diff_err2
    print(f'Chi2: {chi2}')
    chi2_sum = np.nansum(chi2)
    print(f'Chi2_sum: {chi2_sum}')

    return chi2_sum


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
    base_path = 'F:/Research/'
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
    raw_bss, mix_bss = np.array(raw_bss), np.array(mix_bss)

    # with warnings.catch_warnings():  # Ignore divide by zero warnings when raw=mix=0. From 0 padding
    #     warnings.simplefilter('ignore')
    # ---------------------------------------------------------------------------------------
    diff_def = (raw_def - mix_def)  # / (2 * (raw_def + mix_def))
    # ---------------------------------------------------------------------------------------
    bs_diffs = (raw_bss[..., None] - mix_bss.T).transpose(0, 2, 1).reshape(-1, len(raw_bss[0]))
    # raw_bs_tile, mix_bs_repeat = np.tile(raw_bss, (len(raw_bss), 1)), np.repeat(mix_bss, len(mix_bss), axis=0)
    # # ---------------------------------------------------------------------------------------
    # bs_diffs = (raw_bs_tile - mix_bs_repeat)  # / (2 * (raw_bs_tile + mix_bs_repeat))
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
    # prof = pprofile.Profile()
    # with prof():
    raw_def_tp, raw_bss_tp = get_norm_dists_all(raw_path)
    mix_def_tp, mix_bss_tp = get_norm_dists_all(mix_path)

    # print(raw_path)
    # print(mix_path)

    # print(f'raw shape: {np.array(raw_bss_tp[7]).shape} {len(raw_bss_tp[7])}')
    # print(f'mix shape: {np.array(mix_bss_tp[7]).shape} {len(mix_bss_tp[7])}')

    total_protons = get_total_proton_list(raw_def_tp, mix_def_tp, raw_bss_tp, mix_bss_tp)

    total_proton_data = {}
    # with warnings.catch_warnings():  # Ignore divide by zero warnings when raw=mix=0. From 0 padding
    #     warnings.simplefilter('ignore')
    for total_proton in total_protons:
        raw_def, mix_def = raw_def_tp[total_proton], mix_def_tp[total_proton]
        # raw_bss, mix_bss = raw_bss_tp[total_proton], mix_bss_tp[total_proton]
        raw_bss, mix_bss = np.array(raw_bss_tp[total_proton]), np.array(mix_bss_tp[total_proton])
        # print(f'total protons: {total_proton}')
        # print(f'raw shape: {np.array(raw_bss).shape} {len(raw_bss)}')
        # print(f'mix shape: {np.array(mix_bss).shape} {len(mix_bss)}')
        # print(raw_bss.shape[0])
        # print(f'mix shape: {np.array(mix_bss[:raw_bss.shape[0]]).shape} {len(mix_bss)}')

        # ---------------------------------------------------------------------------------------
        diff_def = (raw_def - mix_def)  # / (2 * (raw_def + mix_def))
        # ---------------------------------------------------------------------------------------
        # bs_diffs = np.array([raw_bs - mix_bs for raw_bs in raw_bss for mix_bs in mix_bss])  # 26s to 17s, slower
        # print(raw_bss)
        bs_diffs = (raw_bss[..., None] - mix_bss.T).transpose(0, 2, 1).reshape(-1, len(raw_bss[0]))

        # New test
        diffs_bs_slim = raw_bss - mix_bss[:raw_bss.shape[0]]  # Sometimes fewer raw than mix
        # New test
        # bs_diffs = np.broadcast_to(mix_bss.T, (raw_bss.shape[0], mix_bss.shape[1],
        #                                        mix_bss.shape[0])).transpose(2, 0, 1).reshape(-1, mix_bss.shape[1]) \
        #            - np.broadcast_to(raw_bss, (mix_bss.shape[0], raw_bss.shape[0],
        #                                        raw_bss.shape[1])).reshape(-1, mix_bss.shape[1])
        # raw_bs_tile, mix_bs_repeat = np.tile(raw_bss, (len(mix_bss), 1)), np.repeat(mix_bss, len(raw_bss), axis=0)
        # # ---------------------------------------------------------------------------------------
        # bs_diffs = (raw_bs_tile - mix_bs_repeat)  # / (2 * (raw_bs_tile + mix_bs_repeat))
        # ---------------------------------------------------------------------------------------
        diff_sd = np.std(bs_diffs, axis=0)
        diff_def = np.where((diff_sd == 0) | (np.isnan(diff_sd)), float('nan'), diff_def)
        diffs_bs_slim = np.where((diff_sd == 0) | (np.isnan(diff_sd)), float('nan'), diffs_bs_slim)  # Check this
        total_proton_data.update({total_proton: (diff_def, diff_sd, diffs_bs_slim)})
    # prof.dump_stats('C:/Users/Dylan/Desktop/get_diff_all_profile_comp2.txt')
    # prof.print_stats()

    return total_proton_data


def comp_test():
    base_path = 'F:/Research/'
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
    def_dist = np.array(abin_data.get_dist()[total_protons], dtype=np.longlong)
    def_dist_norm = np.pad(def_dist / np.sum(def_dist), (0, total_protons + 1 - def_dist.size))

    bs_dists_norm = []
    for bs_dist in abin_data.get_dist_bs():
        bs_dist = np.array(bs_dist[total_protons], dtype=np.longlong)
        bs_dist_norm = np.pad(bs_dist / np.sum(bs_dist), (0, total_protons + 1 - bs_dist.size))
        bs_dists_norm.append(bs_dist_norm)

    return def_dist_norm, bs_dists_norm


def get_norm_dists_all(path):
    abin_data = Babd(path=path)
    def_dist = abin_data.get_dist()
    def_dists_norm = {}
    for total_protons, dist in def_dist.items():
        dist = np.array(dist, dtype=np.longlong)
        def_dists_norm.update({total_protons: np.pad(dist / np.sum(dist), (0, total_protons + 1 - dist.size))})

    bs_dists_norm = {}
    for bs_dist in abin_data.get_dist_bs():
        for total_protons, dist in bs_dist.items():
            dist = np.array(dist, dtype=np.longlong)
            bs_dist_norm = np.pad(dist / np.sum(dist), (0, total_protons + 1 - dist.size))
            if total_protons in bs_dists_norm:
                bs_dists_norm[total_protons].append(bs_dist_norm)
            else:
                bs_dists_norm.update({total_protons: [bs_dist_norm]})

    return def_dists_norm, bs_dists_norm


if __name__ == '__main__':
    main()
