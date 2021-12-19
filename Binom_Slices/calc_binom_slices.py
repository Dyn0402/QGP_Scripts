#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 01 2:10 PM 2021
Created in PyCharm
Created as QGP_Scripts/calc_binom_slices.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.optimize import curve_fit as cf
from scipy.interpolate import interp1d
from multiprocessing import Pool
import itertools
import tqdm

import istarmap
from AzimuthBinData import AzimuthBinData
from Bootstrap_Az_Bin import BootstrapAzBin
from DistStats import DistStats
from Measure import Measure
from pickle_methods import *


def main():
    pars = init_pars()
    read_data(pars)

    print('donzo')


def init_pars():
    pars = {'base_path': '/home/dylan/Research/',
            'csv_path': '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df_new.csv',
            'threads': 16,
            'stats': define_stats(['standard deviation', 'skewness', 'non-excess kurtosis']),
            'datasets': define_datasets(),
            'sys_sets': define_sys_sets(),
            }

    return pars


def define_datasets():
    """
    Define sets to read for each dataset. Each dataset can contain multiple sets which will be combined to give
    systematic uncertainties. Here set parameters to define which sets will be read and combined.
    :return: list of dictionaries containing parameters for defining datasets
    """

    all_divs = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    all_energies = [7, 11, 19, 27, 39, 62]

    entry_names = ['name', 'base_ext', 'exact_keys', 'contain_keys', 'exclude_keys',
                   'set_nums', 'energies', 'cents', 'divs']
    entry_vals = [
        ['ampt_def', '_Ampt', ['default'], [], ['resample'], range(10), all_energies, [8], all_divs],
        ['ampt_resample_def', '_Ampt', ['default', 'resample'], [], [], [0], all_energies, [8], all_divs],
        # ['sim_amp01_spread01', '_Sim', ['anticlmulti', 'amp05', 'spread01'], ['single'], [], [0], [62], [8], all_divs]
    ]

    datasets = [dict(zip(entry_names, dset)) for dset in entry_vals]

    return datasets


def define_sys_sets():
    sys_sets = [{'default': {'name': 'ampt_def', 'set_nums': range(10)},
                 'sys': [{'name': 'ampt_def', 'set_nums': range(1)}]},
                {'default': {'name': 'ampt_resample_def', 'set_nums': range(0)},
                 'sys': [{'name': 'ampt_def', 'set_nums': range(1)}]}
                ]

    return sys_sets


def define_stats(stats):
    stat_methods = {'mean': get_mean_meas,
                    'standard deviation': get_sd_meas,
                    'skewness': get_skewness_meas,
                    'kurtosis': get_kurtosis_meas,
                    'non-excess kurtosis': get_nekurtosis_meas,
                    'c1': get_c1_meas,
                    'c2': get_c2_meas,
                    'c3': get_c3_meas,
                    'c4': get_c4_meas,
                    'c5': get_c5_meas,
                    'c6': get_c6_meas,
                    'k1': get_k1_meas,
                    'k2': get_k2_meas,
                    'k3': get_k3_meas,
                    'k4': get_k4_meas,
                    'k5': get_k5_meas,
                    'k6': get_k6_meas,
                    'c4/c2': get_c4_div_c2_meas,
                    'k4/k2': get_k4_div_k2_meas,
                    'c6/c2': get_c6_div_c2_meas,
                    'k6/k2': get_k6_div_k2_meas,
                    'c4/c2 - k4/k2': get_c4_div_c2_sub_k4_div_k2_meas,
                    'c6/c2 - k6/k2': get_c6_div_c2_sub_k6_div_k2_meas,
                    }

    return {stat: stat_methods[stat] for stat in stats}


def read_data(pars):
    jobs = []
    for dataset in pars['datasets']:
        print(dataset)
        jobs.extend(get_dataset_jobs(dataset, pars))

    df_subsets = []
    # pool = Pool(processes=pars['threads'])
    # for df_subset in tqdm.tqdm(pool.starmap(read_subset, jobs), total=len(jobs)):
    #     df_subsets.extend(df_subset)
    with Pool(pars['threads']) as pool:
        for df_subset in tqdm.tqdm(pool.istarmap(read_subset, jobs), total=len(jobs)):
            # print(df_subset)
            df_subsets.extend(df_subset)
            # df_subsets = pool.starmap(read_subset, jobs)

    # df = list(itertools.chain.from_iterable(df_subsets))
    df = pd.DataFrame(df_subsets)
    df_sys = get_systematics(pars, df)
    df_sys.to_csv(pars['csv_path'])


def get_dataset_jobs(dataset, pars):
    path = pars['base_path'] + 'Data' + dataset['base_ext']
    dataset_jobs = []

    for set_dir in os.listdir(path):
        good = True
        if not all(x in set_dir.split('_') for x in dataset['exact_keys']):
            good = False
        if not all(x in set_dir for x in dataset['contain_keys']):
            good = False
        if any(x in set_dir.split('_') for x in dataset['exclude_keys']):
            good = False
        if good:
            # print(set_dir)
            dataset_jobs.extend(get_set_jobs(dataset, set_dir, {'raw': path + '/', 'mix': path + '_Mix/'}, pars['stats']))

    return dataset_jobs


def get_set_jobs(dataset, set_dir, base_paths, stats):
    subset_jobs = []
    for energy in dataset['energies']:
        for div in dataset['divs']:
            for cent in dataset['cents']:
                for sub_set_dir in os.listdir(f'{base_paths["raw"]}{set_dir}'):
                    set_num = int(sub_set_dir.split('_')[-1])
                    if set_num in dataset['set_nums']:
                        path = f'{set_dir}/{sub_set_dir}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
                        if not os.path.exists(f'{base_paths["mix"]}{path}'):
                            print(f'Mixed file missing! {base_paths["mix"]}{path}')
                            continue
                        # print(path)
                        # continue
                        other_columns = {'name': dataset['name'], 'set': set_dir, 'set_num': set_num,
                                         'energy': energy, 'div': div, 'cent': cent}
                        subset_jobs.append((f'{base_paths["raw"]}{path}', f'{base_paths["mix"]}{path}', div,
                                            stats, other_columns))

    return subset_jobs


def read_subset(raw_path, mix_path, div, stats, other_columns):
    # print(raw_path)
    df_subset = []
    raw_az_data = BootstrapAzBin(div, raw_path)
    mix_az_data = BootstrapAzBin(div, mix_path)
    for total_protons in raw_az_data.get_dist():
        if total_protons in mix_az_data.get_dist():
            for stat, stat_method in stats.items():
                other_columns.update({'total_protons': total_protons, 'stat': stat})
                measures = get_div(raw_az_data, mix_az_data, total_protons, stat_method)
                for data_type, meas in zip(['raw', 'mix', 'divide'], measures):
                    df_new = {'data_type': data_type, 'val': meas.val, 'err': meas.err}
                    df_new.update(other_columns)
                    df_subset.append(df_new)
        # else:
        #     print(f'Total_protons {total_protons} exists in raw but not mixed for '
        #           f'{mix_az_data.path}')

    return df_subset


def get_div(raw_az_data, mix_az_data, total_protons, stat_method):
    raw_stat_meas, raw_bs_stats = get_stat(raw_az_data, total_protons, stat_method)
    mix_stat_meas, mix_bs_stats = get_stat(mix_az_data, total_protons, stat_method)

    div_stat_meas = raw_stat_meas / mix_stat_meas
    if raw_bs_stats and mix_bs_stats:
        div_stat_meas.err = np.std([raw / mix if abs(mix) > 0 else float('nan') for raw in raw_bs_stats
                                    for mix in mix_bs_stats])

    return raw_stat_meas, mix_stat_meas, div_stat_meas


def get_stat(az_data, total_protons, stat_method):
    """
    Get statistic from az_data distribution. For error use bootstrap standard deviation if bootstraps exist, else use
    delta theorem. Return list of bootstrap set if they exist
    :param az_data: Bootstrap_Az_Bin data
    :param total_protons: Slice of az_data to take
    :param stat_method: Stat method to get from DistData of az_data[total_protons] distribution
    :return:
    """
    stat_meas = stat_method(DistStats(az_data.get_dist()[total_protons]))
    stat_bs = None
    if len(az_data.data_bs) > 0:
        # if not all(total_protons in bs for bs in az_data.get_dist_bs()):
        #     print(f'Bootstrap missing total protons {total_protons}! {az_data.path}')
        stat_bs = [stat_method(DistStats(bs[total_protons])).val for bs in az_data.get_dist_bs() if total_protons in bs]
        stat_meas.err = np.std(stat_bs)

    return stat_meas, stat_bs


def get_systematics(pars, df):
    """
    Only implemented default and spread due to set_nums, ignore pars['sys_sets']['sys']
    :param pars:
    :param df:
    :return:
    """
    sys_df = []
    # attributes = ['energy', 'div', 'cent', 'total_protons', 'stat', 'data_type']
    # lambda df, att : np.unique(df[att])
    print(df)
    for sys_set in pars['sys_sets']:
        # default
        df_def = df[df['name'] == sys_set['default']['name']]
        for energy in np.unique(df_def['energy']):
            df_energy = df_def[df_def['energy'] == energy]
            for div in np.unique(df_energy['div']):
                df_div = df_energy[df_energy['div'] == div]
                for cent in np.unique(df_div['cent']):
                    df_cent = df_div[df_div['cent'] == cent]
                    for total_protons in np.unique(df_cent['total_protons']):
                        df_tps = df_cent[df_cent['total_protons'] == total_protons]
                        for stat in np.unique(df_tps['stat']):
                            df_stat = df_tps[df_tps['stat'] == stat]
                            for data_type in np.unique(df_stat['data_type']):
                                df_dtype = df_stat[df_stat['data_type'] == data_type]
                                # Should weigh sd with errors
                                if len(df_dtype) > 1:
                                    meases = [Measure(val, err) for val, err in zip(df_dtype['val'], df_dtype['err'])]
                                    med = np.median(meases)
                                    new_row = {'val': med.val, 'err': med.err, 'sys': np.std(df_dtype['val'])}
                                else:
                                    new_row = {'val': df_dtype['val'].iloc[0], 'err': df_dtype['err'].iloc[0]}
                                new_row.update({'name': sys_set['default']['name'], 'energy': energy, 'div': div,
                                                'cent': cent, 'total_protons': total_protons, 'stat': stat,
                                                'data_type': data_type})
                                sys_df.append(new_row)

    return pd.DataFrame(sys_df)


if __name__ == '__main__':
    main()
