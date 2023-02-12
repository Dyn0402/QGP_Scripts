#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 17 4:27 PM 2022
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slice_plotter

@author: Dylan Neff, Dyn04
"""

import pandas as pd
import pickle

import itertools

import tqdm

from analyze_binom_slices import *


def main():
    # plot_sims()
    # get_sim_mapping()
    # get_sim_mapping_pm()
    # plot_star_model()
    # plot_star_model_onediv()
    # plot_vs_cent()
    # plot_closest_sims()
    # plot_vs_cent_nofit()
    # plot_vs_cent_fittest()
    # plot_all_zero_base()
    # plot_ampt_efficiency()
    # plot_flow()
    plot_ampt_v2_closure()
    print('donzo')


def plot_star_model():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'cf_tprotons_fits.csv'
    df_partitions_fits_name = 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST EV $1fm^3$']))

    # data_sets_plt = ['cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['MUSIC+FIST', 'MUSIC+FIST EV']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_sims():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # df_name = 'binom_slice_stats_cent8_sim.csv'
    df_name = 'binom_slice_stats_cent8_sim_test.csv'
    df_path = base_path + df_name
    sim_sets = []

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    amps = ['002', '006', '01']  # ['002', '006', '01']
    # spreads = ['03', '04', '05', '06', '07', '08', '09', '1', '11', '12']
    spreads = ['1', '12']
    # spreads = ['1', '15']
    # amps = ['0', '002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    #         '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5',
    #         '6', '7', '8', '9', '99']
    # spreads = ['001', '01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25',
    #            '275', '3', '325', '35', '4']  # '375' missing '08' and '09'
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    # stat_plot = 'skewness'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 7]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    if stat_plot == 'skewness':
        exclude_divs.append(180)

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={spread}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    # plt.show()
    proton_fits_div = protons_fits[protons_fits['divs'] == div_plt]
    plot_protons_fits_vs_amp(proton_fits_div, all_sets_plt, data_sets_colors, data_sets_labels)
    # plt.show()
    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    # plt.show()
    plot_sigma_fits_interp(sigma_fits)

    plt.show()


def get_sim_mapping():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_sim.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'sim_tprotons_fits.csv'
    df_partitions_fits_name = 'sim_partitions_fits.csv'
    threads = 11
    df_path = base_path + df_name
    sim_sets = []

    print(f'{base_path}{fits_out_base}{df_tproton_fits_name}')

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    # spreads = ['02', '05', '06', '065', '07']
    # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04']
    # spreads = ['05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
            '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5']
    spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25', '275', '3',
               '325', '35', '375', '4']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            # sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={round(spread, 2)}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, False, True, False, None, None)
    jobs = [(df, stat_plot, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    with Pool(threads) as pool:
        for proton_fits_div in tqdm.tqdm(pool.istarmap(stat_vs_protons, jobs), total=len(jobs)):
            # print(proton_fits_div)
            protons_fits.append(proton_fits_div)
    print(protons_fits)
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     print(f'Div {div}')
    #     protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
    #                                        plot=False, fit=True)
    #     protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    protons_fits.to_csv(f'{base_path}{fits_out_base}/{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    interpolations = plot_sigma_fits_interp(sigma_fits)

    with open(f'{base_path}{pickle_map_name}', 'ab') as pickle_file:
        pickle.dump(interpolations, pickle_file)

    # with open(f'{base_path}{pickle_map_name}', 'rb') as pickle_file:
    #     interps_pickle = pickle.load(pickle_file)

    for fig_i in plt.get_fignums():
        plt.figure(fig_i)
        window_title = plt.gcf().canvas.manager.get_window_title()
        plt.savefig(f'{out_dir}{window_title}.png')

    plt.show()


def get_sim_mapping_pm():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_simpm_test.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    df_fits_out_name = 'Base_Zero_Fits/simpm.csv'
    threads = 11
    df_path = base_path + df_name
    sim_sets = []

    aps = ['04']  # ['02', '04']
    ams = ['02']  # ['01', '02']
    sps = ['02', '04']  # ['02', '04']
    sms = ['2']  # ['1', '2']
    for (ap, am, sp, sm) in itertools.product(aps, ams, sps, sms):
        sim_sets.append(f'sim_clmultipm_ampplus{ap}_ampminus{am}_spreadplus{sp}_spreadminus{sm}')

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim']
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        ap, am, sp, sm = get_name_amp_spread_pm(sim_set)
        label = f'A+={ap} σ+={round(sp, 2)} A-={am} σ-={round(sm, 2)}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, True, True, False, None, None)
    jobs = [(df, stat_plot, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    stat_vs_protons(*jobs[0])
    stat_vs_protons(*jobs[1])
    with Pool(threads) as pool:
        for proton_fits_div in tqdm.tqdm(pool.istarmap(stat_vs_protons, jobs), total=len(jobs)):
            protons_fits.append(proton_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    print(protons_fits)
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{df_fits_out_name}', index=False)
    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_star_model_onediv():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies = [7, 11, 19, 27, 39, 62]
    energy_plt = 39
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def',
                     'cfevb342_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'olive']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$',
                                                'MUSIC+FIST EV $3.42fm^3$']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [energy_plt], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies, data_types_plt,
                                   all_sets_plt, plot=False, fit=True)
    print(protons_fits)
    plot_protons_fits_vs_energy(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, title=f'0-5% Centrality, {div_plt}° Partitions, '
                                                                         f'{samples} Samples per Event')

    plt.show()


def plot_vs_cent():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    save_fits = False
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'bes_ampt_cents_tprotons_fits.csv'
    df_partitions_fits_name = 'bes_ampt_cents_partitions_fits.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    # stat_plot = 'non-excess kurtosis'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 39
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    if stat_plot == 'skewness':
        divs_all.pop(180)
        exclude_divs.append(180)

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    # data_sets_plt = ['bes_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black']))
    # data_sets_energies_cmaps = dict(zip(data_sets_plt, ['brg']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, 8, [energy_plt], ['raw', 'mix'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    stat_vs_protons_cents(df, stat_plot, [div_plt], cents, energy_plt, data_types_plt, all_sets_plt, plot=True,
                          fit=True, plot_fit=True, data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    df_slope_fits = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy_plt, data_types_plt, all_sets_plt,
                                          plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                          data_sets_labels=data_sets_labels)

    plot_protons_fits_divs_cents(df_slope_fits, data_sets_plt, fit=True, plot=True, data_sets_colors=data_sets_colors,
                                 data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)

    df_onediv_fits, df_tp_fits, df_divs_fits = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    for energy in energies_fit:
        df_fit = stat_vs_protons_cents(df, stat_plot, [div_plt], cents, energy, data_types_plt, all_sets_plt,
                                       plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                       data_sets_labels=data_sets_labels)
        df_onediv_fits = pd.concat([df_onediv_fits, df_fit])

        df_tp_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy, data_types_plt, all_sets_plt,
                                          plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                          data_sets_labels=data_sets_labels)
        df_tp_fits = pd.concat([df_tp_fits, df_tp_fit])

        df_fits = plot_protons_fits_divs_cents(df_tp_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                               data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)
        # print(f'{energy}GeV {df_fit}')
        df_divs_fits = pd.concat([df_divs_fits, df_fits])

    if save_fits:
        df_tp_fits.to_csv(f'{base_path}{fits_out_base}/{df_tproton_fits_name}', index=False)
        df_divs_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=None,
                          ref_type=None, data_sets_energies_cmaps=data_sets_energies_cmaps)
    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_vs_cent_nofit():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 60
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 62
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_divs_fits = pd.DataFrame()
    for energy in energies_fit:
        df_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy, data_types_plt, all_sets_plt,
                                       plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                       data_sets_labels=data_sets_labels)
        df_fits = plot_protons_fits_divs_cents(df_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                               data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)
        df_divs_fits = pd.concat([df_divs_fits, df_fits])

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_vs_cent_fittest():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'ref'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 60
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 62
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy_plt, data_types_plt, all_sets_plt,
                                   plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                   data_sets_labels=data_sets_labels)
    df_fits = plot_protons_fits_divs_cents(df_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                           data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)

    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=f'{energy_plt}GeV', fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_all_zero_base():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Base_Zero_Fits/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_names = ['cf_fits.csv', 'sim.csv', 'simpm.csv', 'bes_ampt_cents.csv']  #
    df_paths = [base_path + df_name for df_name in df_names]

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 60
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energies = [7, 11, 19, 27, 39, 62, 'sim']
    data_types_plt = ['divide']

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation']))

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_',
                     'sim_clmultipm_']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation', 'Simulation 2 Gaus']))

    df_fits = pd.DataFrame()
    for df_path in df_paths:
        df = pd.read_csv(df_path)
        df = df.dropna()
        df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['data_set'] else row['energy'], axis=1)
        print(f'{df_path}\n{df.head()}\n\n')
        df_fits = pd.concat([df_fits, df])

    plot_base_zeros(df_fits, data_sets_plt, data_sets_labels, data_sets_colors)

    plt.show()


def plot_closest_sims():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Base_Zero_Fits/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_data_name = 'bes_ampt_cents_tprotons_fits.csv'
    df_sim_name = 'sim_tprotons_fits.csv'

    stat = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    cent = 8
    energy = 62
    num_sims_plt = 4
    data_name = 'ampt_new_coal_resample_def'  # 'bes_resample_def'
    data_title = 'AMPT'  # 'BES'

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation']))

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_',
                     'sim_clmultipm_']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation', 'Simulation 2 Gaus']))

    df_data = pd.read_csv(base_path + df_data_name)
    df_data = df_data.dropna()
    df_data = df_data[(df_data['name'] == data_name) & (df_data['energy'] == energy) & (df_data['cent'] == cent)]
    print(df_data)

    df_sim = pd.read_csv(base_path + df_sim_name)
    df_sim = df_sim.dropna()
    # df_sim['energy'] = df_sim.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_diffs = []
    for sim_set in pd.unique(df_sim['name']):
        df_sim_set = df_sim[df_sim['name'] == sim_set]
        diff, diff_weight = 0, 0
        for div in divs_all:
            diff += abs(df_data[df_data['divs'] == div]['slope'].iloc[0] -
                        df_sim_set[df_sim_set['divs'] == div]['slope'].iloc[0])
            weight = 10 if div == 180 else 1
            diff_weight += abs(df_data[df_data['divs'] == div]['slope'].iloc[0] -
                               df_sim_set[df_sim_set['divs'] == div]['slope'].iloc[0]) * weight
        diff_180 = abs(df_data[df_data['divs'] == 180]['slope'].iloc[0] -
                       df_sim_set[df_sim_set['divs'] == 180]['slope'].iloc[0])
        df_diffs.append({'sim_set': sim_set, 'diff': diff, 'diff_weight': diff_weight, 'diff_180': diff_180})
    df_diffs = pd.DataFrame(df_diffs)

    close_sims_diff = list(df_diffs.sort_values(by=['diff'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_diff = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_diff)]])

    close_sims_180 = list(df_diffs.sort_values(by=['diff_180'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_180 = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_180)]])

    close_sims_weight = list(df_diffs.sort_values(by=['diff_weight'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_weight = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_weight)]])

    df_plt_amp1 = pd.concat([df_data, df_sim[(df_sim['amp'] == 0.1) &
                                             (df_sim['spread'] < 3) & (df_sim['spread'] > 0.5)]])
    df_plt_amp1 = df_plt_amp1.sort_values(by=['spread'])
    df_plt_spread1 = pd.concat([df_data, df_sim[(df_sim['spread'] == 1) & (df_sim['amp'] < 0.125)]])
    df_plt_spread1 = df_plt_spread1.sort_values(by=['amp'])

    plot_vs_div_width_comp(df_plt_diff, f'{data_title} {energy} GeV Equal Weight')
    plot_vs_div_width_comp(df_plt_180, f'{data_title} {energy} GeV 180 Only')
    plot_vs_div_width_comp(df_plt_weight, f'{data_title} {energy} GeV 180 Weighted x10')
    plot_vs_div_width_comp(df_plt_amp1, f'{data_title} {energy} GeV A=0.1')
    plot_vs_div_width_comp(df_plt_spread1, f'{data_title} {energy} GeV spread=1')

    plt.show()


def plot_ampt_efficiency():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_ampt_eff.csv'
    save_fits = False
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = 'ampt_eff_tprotons_fits.csv'
    df_partitions_fits_name = 'ampt_eff_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_resample_def', 'ampt_new_coal_resample_eff1', 'ampt_new_coal_resample_eff2',
                     'ampt_new_coal_resample_eff3', 'ampt_new_coal_resample_eff4']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['100%', '90%', '80%', '70%', '60%']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=False)

    # plt.show()
    # return

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    if save_fits:
        df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    def_set = 'ampt_new_coal_resample_def'
    df_fits_ratio = []  # Divide baseline/zeros by 100% default to see difference
    for energy in pd.unique(df_fits['energy']):
        df_fits_energy = df_fits[df_fits['energy'] == energy]
        df_fits_energy_def = df_fits_energy[df_fits_energy['data_set'] == def_set]
        base_def = Measure(df_fits_energy_def['baseline'].iloc[0], df_fits_energy_def['base_err'].iloc[0])
        zero_def = Measure(df_fits_energy_def['zero_mag'].iloc[0], df_fits_energy_def['zero_mag_err'].iloc[0])
        for ds_i, data_set in enumerate(pd.unique(df_fits_energy['data_set'])):
            if data_set != def_set:
                df_fits_e_ds = df_fits_energy[df_fits_energy['data_set'] == data_set]
                base_div = Measure(df_fits_e_ds['baseline'].iloc[0], df_fits_e_ds['base_err'].iloc[0]) / base_def
                zero_div = Measure(df_fits_e_ds['zero_mag'].iloc[0], df_fits_e_ds['zero_mag_err'].iloc[0]) / zero_def
                energy_shifted = energy - 1 + ds_i / 2.0
                df_fits_ratio.append({'data_set': data_set, 'energy': energy_shifted, 'baseline': base_div.val,
                                      'base_err': base_div.err, 'zero_mag': zero_div.val, 'zero_mag_err': zero_div.err})
    plot_slope_div_fits(pd.DataFrame(df_fits_ratio), data_sets_colors, data_sets_labels, ref_line=1)

    plt.show()


def plot_flow():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_flow.csv'
    save_fits = False
    v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = 'flow_tprotons_fits.csv'
    df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    # data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res75_v207', 'flow_resample_res5_v207',
    #                  'flow_resample_res3_v207', 'flow_resample_res15_v207']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['0.99', '0.75', '0.5', '0.3', '0.15']))

    data_sets_plt = ['flow_resample_res15_v207', 'flow_resample_res15_v206', 'flow_resample_res15_v205',
                     'flow_resample_res15_v204', 'flow_resample_res15_v203', 'flow_resample_res15_v202']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green', 'olive']))
    data_sets_labels = dict(zip(data_sets_plt, ['v2=0.07', 'v2=0.06', 'v2=0.05', 'v2=0.04', 'v2=0.03', 'v2=0.02']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # for data_set in data_sets_plt:
    #     df_set = df[df['name'] == data_set]
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], ['raw', 'mix'], all_sets_plt, plot=True, fit=False,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False)

    # flow_vs_v2(df, div_plt, '15', v2_fit_out_dir)

    # plt.show()
    # return

    # for data_set in data_sets_plt:
    #     df_set = df[df['name'] == data_set]
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
    #                     y_ranges={'standard deviation': [0.946, 1.045]})
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], ['raw'], all_sets_plt, plot=True, fit=True,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
    #                     y_ranges={'standard deviation': [0.946, 1.045]})
    #
    #
    # plt.show()
    # return
    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        flow_vs_v2(df, div, '15', v2_fit_out_dir)
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    return
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    if save_fits:
        df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    plt.show()


def plot_ampt_v2_closure():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_ampt_v2_closure.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
    v2_in_dir = f'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/' \
                f'Ampt_rapid05_resample_norotate_epbins1_0/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = None  # 'cf_tprotons_fits.csv'
    df_partitions_fits_name = None  # 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 6
    energies_fit = [7, 11, 19, 27, 39, 62]
    energies_fit = [7, 11, 19, 27]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_rp', 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                     'ampt_new_coal_epbins1_v2rpcor']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['Reaction Plane', 'Flow Included', 'Flow Corrected ep',
                                                'Flow Corrected rp']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}
    v2_vals = read_v2_values(v2_in_dir)  # Dictionary of {cent: v2}
    v2_rp_vals = read_v2_values(v2_in_dir, 'v2_rp')  # Dictionary of {cent: v2}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div, 'ampt_new_coal_epbins1',
                                               'ampt_new_coal_epbins1_v2rpcor', v2_rp_vals[cent_plt],
                                               flow_cor_coefs[div])
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div, 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                                               v2_vals[cent_plt], flow_cor_coefs[div])
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


if __name__ == '__main__':
    main()
