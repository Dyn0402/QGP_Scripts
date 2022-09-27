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
from pandas.errors import EmptyDataError
import os

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from Bootstrap_Az_Bin import BootstrapAzBin
from Measure import Measure
from pickle_methods import *


def main():
    pars = init_pars()
    read_data(pars)

    print('donzo')


def init_pars():
    pars = {
        'base_path': 'F:/Research/',
        # 'base_path': '/media/ucla/Research/',
        # 'csv_path': 'F:/Research/Results/Azimuth_Analysis/binom_slice_stats_cent8_cfev.csv',
        'csv_path': 'F:/Research/Results/Azimuth_Analysis/binom_slice_stats_cent8_no_sim.csv',
        # 'csv_path': '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_stats_cent8_no_sim.csv',
        'csv_append': True,  # If True read dataframe from csv_path and append new datasets to it, else overwrite
        'only_new': True,  # If True check csv_path and only run missing datasets, else run all datasets
        'threads': 1,
        'stats': define_stats(['standard deviation', 'skewness', 'non-excess kurtosis']),
        'check_only': False,  # Don't do any real work, just try to read each file to check for failed reads
        'min_events': 100,  # Min number of total events per total_proton. Skip total_proton if fewer
        'min_bs': 100,  # Min number of bootstrap sets of total_proton. Skip if fewer
        'div_bs': 0,  # Number of bootstrap divide values to get
        'save_cent': False,  # Include centrality column in output dataframe
        'save_data_type': True,  # Include data_type column in output dataframe
        'save_stat': True,  # Include statistic column in output dataframe
        'save_set_num': False,  # Don't save set num if False. Currently not compatible with set nums!
        'systematics': False,  # If True run/save systematics. Currently not implemented!
    }

    pars.update({'datasets': define_datasets(pars['base_path'])})
    if pars['systematics']:
        pars.update({'sys_sets': define_sys_sets(pars['datasets'])})

    return pars


def define_datasets(base_path):
    """
    Define sets to read for each dataset. Each dataset can contain multiple sets which will be combined to give
    systematic uncertainties. Here set parameters to define which sets will be read and combined.
    :param base_path: Base path to data sets, used to find simulation sets
    :return: list of dictionaries containing parameters for defining datasets
    """

    all_divs = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    all_energies = [7, 11, 19, 27, 39, 62]
    # all_energies = [7]
    all_cents = [8]  # [0, 1, 2, 3, 4, 5, 6, 7, 8]

    entry_names = ['name', 'base_ext', 'exact_keys', 'contain_keys', 'exclude_keys',
                   'set_nums', 'energies', 'cents', 'divs']
    entry_vals = [
        # ['ampt_def', '_Ampt', ['default'], [], ['resample'], range(60), all_energies, all_cents, all_divs],
        # ['ampt_new_coal_resample_def', '_Ampt', ['default', 'resample'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_baryon_first_fix_resample_def', '_Ampt', ['default', 'resample'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_old_resample_def', '_Ampt_Old', ['default', 'resample'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_eff1_resample_def', '_Ampt', ['resample', 'Eff1'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_eff2_resample_def', '_Ampt', ['resample', 'Eff2'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_eff3_resample_def', '_Ampt', ['resample', 'Eff3'], [], [], [0], all_energies, all_cents, all_divs],
        # ['ampt_old_eff1_resample_def', '_Ampt_Old', ['resample', 'Eff1'], [], [], [0], all_energies, all_cents,
        #  all_divs],
        # ['ampt_old_eff2_resample_def', '_Ampt_Old', ['resample', 'Eff2'], [], [], [0], all_energies, all_cents,
        #  all_divs],
        # ['ampt_old_eff3_resample_def', '_Ampt_Old', ['resample', 'Eff3'], [], [], [0], all_energies, all_cents,
        #  all_divs],
        # ['bes_def', '', ['default'], [], ['resample'], range(60), all_energies, [8], all_divs],
        # ['bes_resample_def', '', ['default', 'resample'], [], ['alg3'], [0], all_energies, [8], all_divs],
        ['bes_single', '', ['default', 'single'], [], ['alg3'], [0], all_energies, [8], all_divs],
        # ['cf_resample_def', '_CF', ['default', 'resample'], [], ['alg3'], [0], all_energies, all_cents, all_divs],
        # ['cfev_resample_def', '_CFEV', ['default', 'resample'], [], ['alg3'], [0], all_energies, all_cents, all_divs],
        # ['cfevb342_resample_def', '_CFEVb342', ['default', 'resample'], [], ['alg3'], [0], all_energies, all_cents,
        #  all_divs],
    ]

    # # Anti-clustering
    # df = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'], True)
    # for amp in np.unique(df['amp']):
    #     amp_float = float(f'0.{amp}')  # For filtering if needed
    #     if amp_float not in [0.2, 0.5]:
    #         continue
    #     df_amp = df[df['amp'] == amp]
    #     for spread in np.unique(df_amp['spread']):
    #         spread_float = float(f'0.{spread}') * 10  # For filtering if needed
    #         if spread_float not in [0.5, 1]:
    #             continue
    #         entry_vals.append([f'sim_aclmul_amp{amp}_spread{spread}', '_Sim',
    #                            ['anticlmulti', f'amp{amp}', f'spread{spread}', 'resample'],
    #                            ['flat'], [], [0], [62], [8], all_divs])
    #
    # # Clustering
    # df = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'clmulti', 'resample'], ['test'], True)
    # for amp in np.unique(df['amp']):
    #     amp_float = float(f'0.{amp}')  # For filtering if needed
    #     df_amp = df[df['amp'] == amp]
    #     for spread in np.unique(df_amp['spread']):
    #         spread_float = float(f'0.{spread}') * 10  # For filtering if needed
    #         entry_vals.append([f'sim_clmul_amp{amp}_spread{spread}', '_Sim',
    #                            ['clmulti', f'amp{amp}', f'spread{spread}', 'resample'],
    #                            ['flat'], [], [0], [62], [8], all_divs])

    datasets = [dict(zip(entry_names, dset)) for dset in entry_vals]

    return datasets


def define_sys_sets(datasets):
    # sys_sets = [{'default': {'name': 'ampt_def', 'set_nums': range(10)},
    #              'sys': [{'name': 'ampt_def', 'set_nums': range(1)}]},
    #             {'default': {'name': 'ampt_resample_def', 'set_nums': range(0)},
    #              'sys': [{'name': 'ampt_resample_def', 'set_nums': range(0)}]},
    #             {'default': {'name': 'sim_aclmul_amp01_spread2', 'set_nums': range(0)},
    #              'sys': [{'name': 'sim_aclmul_amp01_spread2', 'set_nums': range(1)}]},
    #             {'default': {'name': 'sim_aclmul_amp02_spread2', 'set_nums': range(0)},
    #              'sys': [{'name': 'sim_aclmul_amp02_spread2', 'set_nums': range(1)}]},
    #             ]

    entry_names = ['sys_name', 'default', 'sys_sets']
    entry_vals = [
        ['ampt_def', {'name': 'ampt_def', 'set_nums': range(60)}, []],
        ['bes_def', {'name': 'bes_def', 'set_nums': range(60)}, []],
    ]
    systematic_sets = [dict(zip(entry_names, dset)) for dset in entry_vals]

    sys_names = [sset['sys_name'] for sset in systematic_sets]

    # For each dataset, if not already in systematic_sets and overridden, make blank systematic
    for dataset in datasets:
        if dataset['name'] not in sys_names:
            entry_val = [dataset['name'], {'name': dataset['name'], 'set_nums': range(0)}, []]
            systematic_sets.append(dict(zip(entry_names, entry_val)))

    return systematic_sets


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


def find_sim_sets(path, include_keys, exclude_keys=[], print_sets=False):
    df = []
    for file_path in os.listdir(path):
        file_keys = file_path.strip().split('_')
        if any([key in file_keys for key in exclude_keys]):
            continue
        if all([key in file_keys for key in include_keys]):
            df_i = {}
            for key in file_keys:
                for par_i in ['amp', 'spread']:
                    if key[:len(par_i)] == par_i:
                        df_i.update({par_i: key.strip(par_i)})
            df.append(df_i)

    df = pd.DataFrame(df)

    if print_sets:
        for spread in pd.unique(df['spread']):
            df_spread = df[df['spread'] == spread]
            amps = pd.unique(df_spread["amp"])
            # amps_sq = ['0', '01', '015', '02', '025', '03', '04', '045', '05', '06', '07', '09', '15', '2',
            #            '225', '25', '4', '45', '5']
            # amps = [amp for amp in amps if amp in amps_sq]
            print(f'spread {spread}: {len(amps)} {amps}')

        all_spreads, all_amps = pd.unique(df['spread']), pd.unique(df['amp'])
        print(f'\n All spreads: \n{all_spreads}\n\nAll amps: \n{all_amps}')
        all_amps = set(all_amps)

        print('\nMissing amps per spread')
        missing_sets = []
        for spread in all_spreads:
            df_spread = df[df['spread'] == spread]
            amps_spread = set(pd.unique(df_spread['amp']))
            miss_amps_spread = list(all_amps - amps_spread)
            print(f'spread {spread}: {len(miss_amps_spread)} {miss_amps_spread}')
            missing_sets.extend([(f"'{spread}'", f"'{amp}'") for amp in miss_amps_spread])

        print(f'\nMissing sets: \n{"".join(f"({spread}, {amp}), " for spread, amp in missing_sets)}')

    return df


def read_data(pars):
    if pars['csv_append'] or pars['only_new']:
        try:
            df_old = pd.read_csv(pars['csv_path'])
            old_names = pd.unique(df_old['name'])
        except (FileNotFoundError, EmptyDataError) as e:
            df_old, old_names = None, None
            print(f'{pars["csv_path"]} not found! {e}')

    jobs = []
    for dataset in pars['datasets']:
        if pars['only_new'] and old_names is not None:
            if dataset['name'] in old_names:
                continue
        # print(dataset)
        jobs.extend(get_dataset_jobs(dataset, pars))

    if pars['check_only']:  # Just check the files and return
        with Pool(pars['threads']) as pool:
            for df_subset in tqdm.tqdm(pool.istarmap(check_subset, jobs), total=len(jobs)):
                pass
        return

    df_subsets = []
    with Pool(pars['threads']) as pool:
        for df_subset in tqdm.tqdm(pool.istarmap(read_subset, jobs), total=len(jobs)):
            df_subsets.extend(df_subset)

    df = pd.DataFrame(df_subsets)
    if pars['systematics']:
        df = get_systematics(pars, df)
    if pars['csv_append']:
        if df_old is not None:
            df = df.append(df_old, ignore_index=True)
        else:
            print(f'{pars["csv_path"]} not found! Skipping read and writing new.')

    print(df)
    df.to_csv(pars['csv_path'], index=False)


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
            new_job = get_set_jobs(dataset, set_dir, {'raw': path + '/', 'mix': path + '_Mix/'}, pars)
            # pars['stats'], pars['min_events'], pars['min_bs'], pars['div_bs'])
            dataset_jobs.extend(new_job)

    return dataset_jobs


def get_set_jobs(dataset, set_dir, base_paths, pars):  # stats, min_events, min_bs, div_bs):
    subset_jobs = []
    for energy in dataset['energies']:
        for div in dataset['divs']:
            for cent in dataset['cents']:
                for sub_set_dir in os.listdir(f'{base_paths["raw"]}{set_dir}'):
                    set_num = int(sub_set_dir.split('_')[-1])
                    if set_num in dataset['set_nums']:
                        path = f'{set_dir}/{sub_set_dir}/{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
                        info_path = f'{set_dir}/{sub_set_dir}/{energy}GeV/info.txt'
                        if not os.path.exists(f'{base_paths["mix"]}{path}'):
                            print(f'Mixed file missing! {base_paths["mix"]}{path}')
                            continue
                        other_columns = {'name': dataset['name'], 'energy': energy, 'divs': div}
                        if pars['save_cent']:
                            other_columns.update({'cent': cent})
                        # 'set': set_dir, 'set_num': set_num,  # Without set_num pre-resample systematics won't work!
                        if 'sim_' in dataset['name']:
                            other_columns['energy'] = 0
                            amp, spread = get_name_amp_spread(dataset['name'])
                            other_columns.update({'amp': amp, 'spread': spread})
                        else:
                            other_columns.update({'amp': 0, 'spread': 0})
                        subset_jobs.append((f'{base_paths["raw"]}{path}', f'{base_paths["mix"]}{path}',
                                            f'{base_paths["raw"]}{info_path}', div, pars["stats"], other_columns,
                                            pars["min_events"], pars["min_bs"], pars["div_bs"],
                                            pars["save_data_type"], pars["save_stat"]))

    return subset_jobs


def read_subset(raw_path, mix_path, info_path, div, stats, other_columns, min_events, min_bs, div_bs,
                save_data_type, save_stat):
    """
    Read a single raw/mix file pair. Calculate stats with bootstrap or delta errors depending on bootstrap's existence.
    Calculate raw divided by mixed. Return as entries to convert to dataframe.
    :param raw_path: Path to raw histogram file
    :param mix_path: Path to mix histogram file
    :param info_path: Path to info text file in raw directory, for determining type of binning and min_counts
    :param div: Azimuthal division width
    :param stats: Stats to calculate for each distribution
    :param other_columns: Other columns characterizing this job, to be added to dataframe
    :param min_events: Minimum required events for total_proton dataset to be kept, otherwise don't process
    :param min_bs: Minimum required bootstrap sets for total_proton dataset to be kept
    :param div_bs: If true get divide bootstraps
    :param save_data_type: If True save datatype as column in output dataframe. If false only save divide data
    :param save_stat: If True save statistic as column in output dataframe.
    :return:
    """
    df_subset = []
    raw_az_data = BootstrapAzBin(div, raw_path)
    mix_az_data = BootstrapAzBin(div, mix_path)

    min_counts = get_min_counts(info_path, min_events, div)

    for total_protons in raw_az_data.get_dist():
        if total_protons in mix_az_data.get_dist():
            for stat, stat_method in stats.items():
                other_columns.update({'total_protons': total_protons})
                if save_stat:
                    other_columns.update({'stat': stat})
                measures = get_div(raw_az_data, mix_az_data, total_protons, stat_method, min_counts, min_bs, div_bs)
                measures, div_bs_sets = measures[:3], measures[-1]
                if any(x is None for x in measures):
                    continue
                if save_data_type:
                    for data_type, meas in zip(['raw', 'mix', 'divide'], measures):
                        df_new = {'data_type': data_type, 'val': meas.val, 'err': meas.err}
                        if div_bs:
                            df_new.update({'bs_index': -1})
                        df_new.update(other_columns)
                        df_subset.append(df_new)
                else:  # Only save divide val/err
                    df_new = {'val': measures[-1].val, 'err': measures[-1].err}
                    if div_bs:
                        df_new.update({'bs_index': -1})
                    df_new.update(other_columns)
                    df_subset.append(df_new)
                if div_bs:
                    for bs_index, div_bs_meas in enumerate(div_bs_sets):
                        df_new = {'val': div_bs_meas.val, 'err': div_bs_meas.err}  # 'data_type': f'divide_{bs_index}',
                        if save_data_type:
                            df_new.update({'data_type': 'divide'})
                        df_new.update(other_columns)
                        df_new.update({'bs_index': bs_index})
                        df_subset.append(df_new)
        # else:
        #     print(f'Total_protons {total_protons} exists in raw but not mixed for '
        #           f'{mix_az_data.path}')

    return df_subset


def check_subset(raw_path, mix_path, info_path, div, stats, other_columns, min_events, min_bs, div_bs,
                 save_data_type, save_stat):
    """
    Just read files to see if all are able to be read. If not hopefully get print to screen with bad path
    :param raw_path:
    :param mix_path:
    :param div:
    :return:
    """
    BootstrapAzBin(div, raw_path)
    BootstrapAzBin(div, mix_path)


def get_min_counts(info_path, min_events, div):
    """
    Use info text file to figure out what type of binning done. Since multiple counts taken from single event depending
    upon the binning, multiply counts per event by min_events to get min_counts.
    !!! Don't check the mixed values, just assume they are the same as raw !!!
    :param info_path: Path to the info file in the raw directory
    :param min_events: Minimum number of events required to keep dataset. To be converted to min_counts
    :param div: Azimuthal division width, needed to calculate number of counts per event for n1_ratios and default
    :return: min_counts -> Minimum number of counts required to keep dataset
    """
    min_counts = min_events
    with open(info_path, 'r') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]
        if 'resample: true' in lines:
            for line in lines:
                if line[:len('n_resamples')] == 'n_resamples':
                    min_counts *= int(line.strip().split(': ')[-1])
        elif 'n1_ratios: true' in lines:
            min_counts *= int(360 / div) - 1
        elif 'single_ratio: true':
            pass  # 1 count per event
        else:  # Default take all bins
            min_counts *= int(360 / div)

    return min_counts


def get_div(raw_az_data, mix_az_data, total_protons, stat_method, min_counts, min_bs, div_bs):
    """
    Calculate raw divided by mixed data for given statistic.
    :param raw_az_data: Raw binning histogram data, default and possibly bootstrap sets
    :param mix_az_data: Mixed binning histogram data, default and possibly bootstrap sets
    :param total_protons: Total proton number (bin) to calculate statistics for. Only one per function call
    :param stat_method: Method for calculating desired statistic on the raw and mixed distributions
    :param min_counts: Minimum number of counts in a total_protons bin necessary. Otherwise discard bin
    :param min_bs: Minimum nuber of bootstrap sets necessary. Otherwise discard
    :param div_bs: If greater than 0 get divide bootstraps. Value corresponds to number of bootstraps to get
    :return: raw, mixed, and divide measures for given statistic on distributions.
             Also set of partially bootstrapped divided values.
             Use raw_bs(x)/mix_bs(x) for default, std(raw_bs(x)/mix_bs(y), for all y) as err
    """
    raw_stat_meas, raw_bs_stats = get_stat(raw_az_data, total_protons, stat_method, min_counts, min_bs)
    mix_stat_meas, mix_bs_stats = get_stat(mix_az_data, total_protons, stat_method, min_counts, min_bs)

    if raw_stat_meas and mix_stat_meas:
        div_stat_meas = raw_stat_meas / mix_stat_meas
        if raw_bs_stats and mix_bs_stats:
            div_list = [raw / mix if abs(mix) > 0 else float('nan') for raw in raw_bs_stats for mix in mix_bs_stats]
            # print(f'total_protons: {total_protons}')
            ds = DistStats(div_list, unbinned=True)
            div_stat_meas.err = ds.get_sd().val
            # mean = np.mean(div_list)
            # div_stat_meas.err = np.std([raw / mix if abs(mix) > 0 else float('nan') for raw in raw_bs_stats
            #                             for mix in mix_bs_stats])
            if div_bs:
                div_stat_bs_list = []
                n_vals = ds.get_total_counts()
                raw2_sum = ds.get_raw_moment(2) * n_vals
                raw1_sum = ds.get_raw_moment(1) * n_vals  # Total Sum (mean * n_vals)
                n_bs = len(raw_bs_stats)
                for i in range(n_bs)[:div_bs]:  # Assume raw_bs_stats and mix_bs_stats same size. Get first div_bs
                    raw_i, mix_i = raw_bs_stats[i], mix_bs_stats[i]
                    if not abs(mix_i) > 0:
                        continue
                    div_stat_bs_meas = Measure(raw_i / mix_i)
                    raw_remove_list = [raw_i / mix_bs_stats[j] if abs(mix_bs_stats[j]) > 0
                                       else float('nan') for j in range(n_bs) if j != i]
                    mix_remove_list = [raw_bs_stats[j] / mix_i for j in range(n_bs) if j != i]
                    remove_array = np.array(raw_remove_list + mix_remove_list)
                    n_vals_new = n_vals - remove_array.size
                    raw_1 = (raw1_sum - np.sum(remove_array)) / n_vals_new
                    raw_2 = (raw2_sum - np.sum(remove_array ** 2)) / n_vals_new
                    div_stat_bs_meas.err = np.sqrt(raw_2 - raw_1 ** 2)
                    div_stat_bs_list.append(div_stat_bs_meas)

                # for index, raw in enumerate(raw_bs_stats):  # Assume raw_bs_stats and mix_bs_stats same size
                #     div_stat_bs_meas = Measure(raw / mix_bs_stats[index])
                #     remove_list = [raw / mix if abs(mix) > 0 else float('nan') for i, j in
                #                    zip(range(len(raw_bs_stats)), len(mix_bs_stats))]
                #     div_stat_bs_meas.err = np.std([raw / mix if abs(mix) > 0 else float('nan') for mix in mix_bs_stats])
                #     div_stat_bs_list.append(div_stat_bs_meas)
            else:
                div_stat_bs_list = None
        else:
            div_stat_bs_list = None
    else:
        div_stat_meas = None
        div_stat_bs_list = None

    return raw_stat_meas, mix_stat_meas, div_stat_meas, div_stat_bs_list


def get_stat(az_data, total_protons, stat_method, min_counts, min_bs):
    """
    Get statistic from az_data distribution. For error use bootstrap standard deviation if bootstraps exist, else use
    delta theorem. Return list of bootstrap set if they exist
    :param az_data: Bootstrap_Az_Bin data
    :param total_protons: Slice of az_data to take
    :param stat_method: Stat method to get from DistData of az_data[total_protons] distribution
    :param min_counts: Minimum number of counts required to keep data. Otherwise return None
    :param min_bs: Minimum number of bootstrap sets required to keep data. Otherwise return None
    :return: Measure object for data with boostrap as err if enough, otherwise delta. Also boostrap float itself
    """
    # print(f'{total_protons}: {sum(az_data.get_dist()[total_protons]) / 1440}, {az_data.get_dist()[total_protons]}')
    if sum(az_data.get_dist()[total_protons]) > min_counts:
        stat_meas = stat_method(DistStats(az_data.get_dist()[total_protons]))
    else:
        stat_meas = None
    stat_bs = None
    num_bootstraps = sum(total_protons in bs for bs in az_data.get_dist_bs())
    if stat_meas and num_bootstraps > min_bs:
        stat_bs = [stat_method(DistStats(bs[total_protons])).val for bs in az_data.get_dist_bs() if total_protons in bs]
        stat_meas.err = np.std(stat_bs)
        # if num_bootstraps < len(az_data.get_dist_bs()):
        #     print(f'Bootstrap missing total protons {total_protons}! {num_bootstraps}/{len(az_data.get_dist_bs())} '
        #           f'{az_data.path}')
        # fig, ax = plt.subplots()
        # ax.set_title(f'Total Protons {total_protons} {num_bootstraps}/{len(az_data.get_dist_bs())}')
        # ax.hist([stat_method(DistStats(bs[total_protons])).val for bs in az_data.get_dist_bs()
        #         if total_protons in bs])
        # ax.axvline(stat_meas.val, ls='--', color='green')
        # plt.show()

    return stat_meas, stat_bs


def get_systematics(pars, df):
    """
    Only implemented default and spread due to set_nums, ignore pars['sys_sets']['sys']
    :param pars: All parameters for script
    :param df: Dataframe of all individual sets/subsets
    :return: Dataframe of all systematic sets
    """
    sys_df = []

    jobs = [(df[df['name'] == sys_set['default']['name']], sys_set) for sys_set in pars['sys_sets']]
    with Pool(pars['threads']) as pool:
        for sys_set_list in tqdm.tqdm(pool.istarmap(get_sys_set, jobs), total=len(jobs)):
            sys_df.extend(sys_set_list)

    return pd.DataFrame(sys_df)


def get_sys_set(df, sys_set):
    sys_set_list = []
    df_def = df[df['name'] == sys_set['default']['name']]
    for energy in np.unique(df_def['energy']):
        df_energy = df_def[df_def['energy'] == energy]
        for div in np.unique(df_energy['divs']):
            df_div = df_energy[df_energy['divs'] == div]
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
                                # For now set sys to 0 if only one set_num. Then don't have issues with NA values
                                new_row = {'val': df_dtype['val'].iloc[0], 'err': df_dtype['err'].iloc[0], 'sys': 0}
                            new_row.update({'name': sys_set['sys_name'], 'energy': energy, 'divs': div,
                                            'cent': cent, 'total_protons': total_protons, 'stat': stat,
                                            'data_type': data_type})
                            if 'sim_' in sys_set['sys_name']:
                                amp, spread = get_name_amp_spread(sys_set['sys_name'])
                                new_row.update({'amp': amp, 'spread': spread})
                            else:
                                new_row.update({'amp': 0, 'spread': 0})
                            sys_set_list.append(new_row)

    return sys_set_list


def get_name_amp_spread(name, as_type='float'):
    name = name.split('_')
    for x in name:
        if 'amp' == x[:len('amp')]:
            if as_type == 'float':
                amp = float(f'0.{x.strip("amp")}')
            else:
                amp = x.strip('amp')
        if 'spread' == x[:len('spread')]:
            if as_type == 'float':
                spread = float(f'0.{x.strip("spread")}') * 10
            else:
                spread = x.strip('spread')

    return amp, spread


if __name__ == '__main__':
    main()
