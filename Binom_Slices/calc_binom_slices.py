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
    pars = {'base_path': 'D:/Research/',  # '/home/dylan/Research/',
            'csv_path': 'D:/Research/Results/Azimuth_Analysis/binom_slice_cent_sds_df_new.csv',  # '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv',
            'csv_append': True,  # If true read dataframe from csv_path and append new datasets to it, else overwrite
            'threads': 16,
            'stats': define_stats(['standard deviation']),  # , 'skewness', 'non-excess kurtosis']),
            'check_only': False,  # Don't do any real work, just try to read each file to check for failed reads
            'min_events': 100,  # Min number of total events per total_proton. Skip total_proton if fewer
            'min_bs': 100,  # Min number of bootstrap sets of total_proton. Skip if fewer
            'div_bs': False,  # If true get divide bootstrap values
            }

    pars.update({'datasets': define_datasets(pars['base_path'])})
    pars.update({'sys_sets': define_sys_sets(pars['datasets'])})

    return pars


def define_datasets(base_path):
    """
    Define sets to read for each dataset. Each dataset can contain multiple sets which will be combined to give
    systematic uncertainties. Here set parameters to define which sets will be read and combined.
    :param base_path: Base path to data sets
    :return: list of dictionaries containing parameters for defining datasets
    """

    all_divs = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    all_energies = [7, 11, 19, 27, 39, 62]
    all_cents = [8]  # [0, 1, 2, 3, 4, 5, 6, 7, 8]

    entry_names = ['name', 'base_ext', 'exact_keys', 'contain_keys', 'exclude_keys',
                   'set_nums', 'energies', 'cents', 'divs']
    entry_vals = [
        # ['ampt_def', '_Ampt', ['default'], [], ['resample'], range(60), all_energies, all_cents, all_divs],
        # ['ampt_resample_def', '_Ampt', ['default', 'resample'], [], [], [0], all_energies, all_cents, all_divs],
        # ['bes_def', '', ['default'], [], ['resample'], range(60), all_energies, [8], all_divs],
        # ['bes_resample_def', '', ['default', 'resample'], [], [], [0], all_energies, [8], all_divs],
    ]

    df = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'])

    for amp in np.unique(df['amp']):
        if amp not in ['275', '3', '325', '35', '375', '4', '45', '5']:
            continue
        df_amp = df[df['amp'] == amp]
        for spread in np.unique(df_amp['spread']):
            entry_vals.append([f'sim_aclmul_amp{amp}_spread{spread}', '_Sim',
                               ['anticlmulti', f'amp{amp}', f'spread{spread}', 'resample'],
                               ['flat'], [], [0], [62], [8], all_divs])

    # df_tests = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample', 'test'])
    # for amp in np.unique(df_tests['amp']):
    #     df_amp = df_tests[df_tests['amp'] == amp]
    #     for spread in np.unique(df_amp['spread']):
    #         entry_vals.append([f'sim_aclmul_amp{amp}_spread{spread}_test', '_Sim',
    #                            ['anticlmulti', f'amp{amp}', f'spread{spread}', 'resample', 'test'],
    #                            ['flat'], [], [0], [62], [8], all_divs])
    #
    # for amp in np.unique(df['amp']):
    #     if amp in np.unique(df_tests['amp']):
    #         df_amp = df[df['amp'] == amp]
    #         for spread in np.unique(df_amp['spread']):
    #             entry_vals.append([f'sim_aclmul_amp{amp}_spread{spread}', '_Sim',
    #                                ['anticlmulti', f'amp{amp}', f'spread{spread}', 'resample'],
    #                                ['flat'], ['test'], [0], [62], [8], all_divs])

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


def find_sim_sets(path, include_keys, exclude_keys=[]):
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
    for spread in np.unique(df['spread']):
        df_spread = df[df['spread'] == spread]
        print(f'spread {spread}: {np.unique(df_spread["amp"])}')

    return df


def read_data(pars):
    jobs = []
    for dataset in pars['datasets']:
        print(dataset)
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
    df_sys = get_systematics(pars, df)
    if pars['csv_append']:
        try:
            df_old = pd.read_csv(pars['csv_path'])
            df_sys = df_sys.append(df_old, ignore_index=True)
        except FileNotFoundError:
            print(f'{pars["csv_path"]} not found! Skipping read and writing new.')

    print(df_sys)
    df_sys.to_csv(pars['csv_path'], index=False)


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
            new_job = get_set_jobs(dataset, set_dir, {'raw': path + '/', 'mix': path + '_Mix/'}, pars['stats'],
                                   pars['min_events'], pars['min_bs'], pars['div_bs'])
            dataset_jobs.extend(new_job)

    return dataset_jobs


def get_set_jobs(dataset, set_dir, base_paths, stats, min_events, min_bs, div_bs):
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
                        other_columns = {'name': dataset['name'], 'set': set_dir, 'set_num': set_num,
                                         'energy': energy, 'divs': div, 'cent': cent}
                        if 'sim_' in dataset['name']:
                            other_columns['energy'] = 0
                        subset_jobs.append((f'{base_paths["raw"]}{path}', f'{base_paths["mix"]}{path}',
                                            f'{base_paths["raw"]}{info_path}', div, stats, other_columns,
                                            min_events, min_bs, div_bs))

    return subset_jobs


def read_subset(raw_path, mix_path, info_path, div, stats, other_columns, min_events, min_bs, div_bs):
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
    :return:
    """
    df_subset = []
    raw_az_data = BootstrapAzBin(div, raw_path)
    mix_az_data = BootstrapAzBin(div, mix_path)

    min_counts = get_min_counts(info_path, min_events, div)

    for total_protons in raw_az_data.get_dist():
        if total_protons in mix_az_data.get_dist():
            for stat, stat_method in stats.items():
                other_columns.update({'total_protons': total_protons, 'stat': stat})
                measures = get_div(raw_az_data, mix_az_data, total_protons, stat_method, min_counts, min_bs, div_bs)
                measures, div_bs_sets = measures[:3], measures[-1]
                if any(x is None for x in measures):
                    continue
                for data_type, meas in zip(['raw', 'mix', 'divide'], measures):
                    df_new = {'data_type': data_type, 'val': meas.val, 'err': meas.err}
                    df_new.update(other_columns)
                    df_subset.append(df_new)
                if div_bs:
                    for bs_index, div_bs_meas in enumerate(div_bs_sets):
                        df_new = {'data_type': f'divide_{bs_index}', 'val': div_bs_meas.val, 'err': div_bs_meas.err}
                        df_new.update(other_columns)
                        df_subset.append(df_new)
        # else:
        #     print(f'Total_protons {total_protons} exists in raw but not mixed for '
        #           f'{mix_az_data.path}')

    return df_subset


def check_subset(raw_path, mix_path, info_path, div, stats, other_columns, min_events, min_bs, div_bs):
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
    :param div_bs: If true get divide bootstraps
    :return: raw, mixed, and divide measures for given statistic on distributions.
             Also set of partially bootstrapped divided values.
             Use raw_bs(x)/mix_bs(x) for default, std(raw_bs(x)/mix_bs(y), for all y) as err
    """
    raw_stat_meas, raw_bs_stats = get_stat(raw_az_data, total_protons, stat_method, min_counts, min_bs)
    mix_stat_meas, mix_bs_stats = get_stat(mix_az_data, total_protons, stat_method, min_counts, min_bs)

    if raw_stat_meas and mix_stat_meas:
        div_stat_meas = raw_stat_meas / mix_stat_meas
        if raw_bs_stats and mix_bs_stats:
            div_stat_meas.err = np.std([raw / mix if abs(mix) > 0 else float('nan') for raw in raw_bs_stats
                                        for mix in mix_bs_stats])
            if div_bs:
                div_stat_bs_list = []
                for index, raw in enumerate(raw_bs_stats):  # Assume raw_bs_stats and mix_bs_stats same size
                    div_stat_bs_meas = Measure(raw / mix_bs_stats[index])
                    div_stat_bs_meas.err = np.std([raw / mix if abs(mix) > 0 else float('nan') for mix in mix_bs_stats])
                    div_stat_bs_list.append(div_stat_bs_meas)
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
    :param pars:
    :param df:
    :return:
    """
    sys_df = []
    # attributes = ['energy', 'div', 'cent', 'total_protons', 'stat', 'data_type']
    # lambda df, att : np.unique(df[att])
    # print(df)

    jobs = [(df[df['name'] == sys_set['default']['name']], sys_set) for sys_set in pars['sys_sets']]
    with Pool(pars['threads']) as pool:
        for sys_set_list in tqdm.tqdm(pool.istarmap(get_sys_set, jobs), total=len(jobs)):
            sys_df.extend(sys_set_list)

    # for sys_set in pars['sys_sets']:
    #     # default
    #     df_def = df[df['name'] == sys_set['default']['name']]
    #     for energy in np.unique(df_def['energy']):
    #         df_energy = df_def[df_def['energy'] == energy]
    #         for div in np.unique(df_energy['divs']):
    #             df_div = df_energy[df_energy['divs'] == div]
    #             for cent in np.unique(df_div['cent']):
    #                 df_cent = df_div[df_div['cent'] == cent]
    #                 for total_protons in np.unique(df_cent['total_protons']):
    #                     df_tps = df_cent[df_cent['total_protons'] == total_protons]
    #                     for stat in np.unique(df_tps['stat']):
    #                         df_stat = df_tps[df_tps['stat'] == stat]
    #                         for data_type in np.unique(df_stat['data_type']):
    #                             df_dtype = df_stat[df_stat['data_type'] == data_type]
    #                             # Should weigh sd with errors
    #                             if len(df_dtype) > 1:
    #                                 meases = [Measure(val, err) for val, err in zip(df_dtype['val'], df_dtype['err'])]
    #                                 med = np.median(meases)
    #                                 new_row = {'val': med.val, 'err': med.err, 'sys': np.std(df_dtype['val'])}
    #                             else:
    #                                 # For now set sys to 0 if only one set_num. Then don't have issues with NA values
    #                                 new_row = {'val': df_dtype['val'].iloc[0], 'err': df_dtype['err'].iloc[0], 'sys': 0}
    #                             new_row.update({'name': sys_set['sys_name'], 'energy': energy, 'divs': div,
    #                                             'cent': cent, 'total_protons': total_protons, 'stat': stat,
    #                                             'data_type': data_type})
    #                             if 'sim_' in sys_set['sys_name']:
    #                                 amp, spread = get_name_amp_spread(sys_set['sys_name'])
    #                                 new_row.update({'amp': amp, 'spread': spread})
    #                             else:
    #                                 new_row.update({'amp': 0, 'spread': 0})
    #                             sys_df.append(new_row)

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


def get_name_amp_spread(name):
    name = name.split('_')
    for x in name:
        if 'amp' in x:
            amp = float(f'0.{x.strip("amp")}')
        if 'spread' in x:
            spread = float(f'0.{x.strip("spread")}') * 10

    return amp, spread


if __name__ == '__main__':
    main()
