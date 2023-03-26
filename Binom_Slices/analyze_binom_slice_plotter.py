#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 17 4:27 PM 2022
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slice_plotter

@author: Dylan Neff, Dyn04
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

import itertools

import tqdm

from analyze_binom_slices import *


def main():
    # plot_star_model_var()
    # plot_vs_cent_var()
    # plot_sims_var()
    get_sim_mapping_var()

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
    # plot_flow_k2()
    # plot_ampt_v2_closure()
    # plot_ampt_v2_closure_var()
    # plot_flow_v2_closure()
    # plot_flow_v2_closure_raw()
    # plot_flow_eff_test()
    # plot_anticl_flow_closure_test()
    # plot_efficiency_closure_tests()
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


def plot_star_model_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
                     'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    df_name = 'binom_slice_stats_cent8_var_epbins1.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def']
    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_epbins1', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST EV $1fm^3$']))

    data_sets_plt = ['bes_resample_epbins1', 'ampt_new_coal_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}
    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # df = df[df['name'].str.contains('bes')]

    df = df[df['stat'] == stat_plot]
    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])

    stat_binom_vs_protons(df, stat_plot, div_plt, cent_plt, 39, ['raw', 'mix'], 'ampt_new_coal_epbins1')

    dvar_vs_protons(pd.concat([df_raw, df_mix, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)

    dvar_vs_protons_energies(df_diff, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], all_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.0014, 0.0011])

    dsig_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsig_avgs_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets_plt, plot=False,
                                            avg=True)
        dsig_avgs_div_raw.loc[:, 'name'] = dsig_avgs_div_raw['name'] + '_raw'
        dsig_avgs.append(dsig_avgs_div_raw)

        dsig_avgs_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets_plt, plot=False,
                                            avg=True)
        dsig_avgs_div_mix.loc[:, 'name'] = dsig_avgs_div_mix['name'] + '_mix'
        dsig_avgs.append(dsig_avgs_div_mix)

        dsig_avgs_div_diff = dvar_vs_protons(df_diff, div, cent_plt, energies_fit, ['diff'], all_sets_plt, plot=False,
                                             avg=True)
        dsig_avgs_div_diff.loc[:, 'name'] = dsig_avgs_div_diff['name'] + '_sub'

        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'ampt_new_coal_epbins1_sub',
                                                  'ampt_new_coal_epbins1', v2_ampt_vals, div, cent_plt)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'bes_resample_epbins1_sub',
                                                  'bes_resample_epbins1', v2_star_vals, div, cent_plt)
        dsig_avgs.append(dsig_avgs_div_diff)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_raw', '_mix', '_sub']]
        colors = dict(zip(data_sets, ['blue', 'green', 'red']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Raw', ' Mix', ' Sub']]))
        plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels,
                            ylab=r'$\widebar{\Delta\sigma^2}$')
    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_sub', '']]
        colors = dict(zip(data_sets, ['blue', 'red']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Original', ' v2 Corrected']]))
        plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels)
    plot_dvar_avgs_divs(dsig_avgs, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels)
    # if df_tproton_fits_name:
    #     dsig_avgs.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    # print(dsig_avgs)
    # print(pd.unique(dsig_avgs['amp']))
    # print(pd.unique(dsig_avgs['spread']))
    # df_fits = plot_dsig_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                                  data_sets_labels=data_sets_labels)
    # if df_partitions_fits_name:
    #     df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
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
    plt.show()
    plot_sigma_fits_interp(sigma_fits)

    plt.show()


def plot_sims_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_stats_sim.csv'
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
            # sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    # sim_sets = sorted(sim_sets, reverse=True)
    # sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'k2'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 7]
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
        label += f'A={amp} σ={spread}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_diff = calc_dsigma(df)

    dvar_vs_protons(df_diff, div_plt, cent_plt, ['sim'], ['diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)

    dsig_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsig_avgs_div_diff = dvar_vs_protons(df_diff, div, cent_plt, energies_fit, ['diff'], all_sets_plt, plot=False,
                                             avg=True)
        dsig_avgs.append(dsig_avgs_div_diff)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    print(dsig_avgs)
    print(pd.unique(dsig_avgs['amp']))
    print(pd.unique(dsig_avgs['spread']))
    plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels)
    plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                        data_sets_labels=data_sets_labels)
    # # plt.show()
    # proton_fits_div = protons_fits[protons_fits['divs'] == div_plt]
    # plot_protons_fits_vs_amp(proton_fits_div, all_sets_plt, data_sets_colors, data_sets_labels)
    # # plt.show()
    # # print(df_fits)
    # sigma_fits = plot_slope_div_fits_simpars(df_fits)
    # plt.show()
    # plot_sigma_fits_interp(sigma_fits)

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


def get_sim_mapping_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_stats_sim.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'sim_tprotons_fits.csv'
    df_partitions_fits_name = 'sim_partitions_fits.csv'
    threads = 10
    df_path = base_path + df_name
    sim_sets = []

    print(f'{base_path}{fits_out_base}{df_tproton_fits_name}')

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    # spreads = ['02', '05', '06', '065', '07']
    # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04']
    # spreads = ['05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    #         '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5']
    # spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25', '275', '3',
    #            '325', '35', '375', '4']
    amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
            '045', '05']
    spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            # sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'k2'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 62]
    data_types_plt = ['diff']
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

    df_diff = calc_dsigma(df[df['stat'] == stat_plot])

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    dsig_avgs = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, False, True, False, None, None, None)
    jobs = [(df_diff, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    with Pool(threads) as pool:
        for dsig_avgs_div in tqdm.tqdm(pool.istarmap(dvar_vs_protons, jobs), total=len(jobs)):
            # print(dsig_avgs_div)
            dsig_avgs.append(dsig_avgs_div)
    print(dsig_avgs)
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     print(f'Div {div}')
    #     protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
    #                                        plot=False, fit=True)
    #     protons_fits.append(protons_fits_div)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    # dsig_avgs.to_csv(f'{base_path}{fits_out_base}/{df_tproton_fits_name}', index=False)
    print(dsig_avgs)
    print(pd.unique(dsig_avgs['amp']))
    print(pd.unique(dsig_avgs['spread']))
    df_fits = plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                  data_sets_labels=data_sets_labels, plt_energies=False)
    # # df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)
    #
    print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    interpolations = plot_sigma_fits_interp(sigma_fits)
    #
    # with open(f'{base_path}{pickle_map_name}', 'ab') as pickle_file:
    #     pickle.dump(interpolations, pickle_file)
    #
    # # with open(f'{base_path}{pickle_map_name}', 'rb') as pickle_file:
    # #     interps_pickle = pickle.load(pickle_file)
    #
    # for fig_i in plt.get_fignums():
    #     plt.figure(fig_i)
    #     window_title = plt.gcf().canvas.manager.get_window_title()
    #     plt.savefig(f'{out_dir}{window_title}.png')

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


def plot_vs_cent_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
                     'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    df_name = 'binom_slice_stats_cent8_var_epbins1.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'F:/Research/Results/Azimuth_Analysis/{cent_ref_name}')
    ref_type = 'refn'  # 'refn'

    sim_sets = []

    stat_plot = 'k2'
    div_plt = 120
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 39
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['diff']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_epbins1', 'ampt_new_coal_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}
    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df = df[df['stat'] == stat_plot]
    df_diff = calc_dsigma(df)

    dvar_vs_protons(df_diff, div_plt, 8, [39], ['diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons_cents(df_diff, [div_plt], cents, energy_plt, data_types_plt, all_sets_plt, plot=True,
                          avg=True, plot_avg=True, data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels,
                          y_ranges=[-0.0051, 0.0013])

    dsig_avgs = []
    for energy in energies_fit:
        print(f'Energy {energy}GeV')
        dsig_avgs_div_diff = dvar_vs_protons_cents(df_diff, divs_all, cents, energy, ['diff'], all_sets_plt,
                                                   plot=False, avg=True)
        dsig_avgs_div_diff.loc[:, 'name'] = dsig_avgs_div_diff['name'] + '_flow_included'
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'ampt_new_coal_epbins1_flow_included',
                                                  'ampt_new_coal_epbins1', v2_ampt_vals)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'bes_resample_epbins1_flow_included',
                                                  'bes_resample_epbins1', v2_star_vals)
        dsig_avgs.append(dsig_avgs_div_diff)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)

    for energy in energies_fit:
        dsig_avgs_energy = dsig_avgs[dsig_avgs['energy'] == energy]
        plot_dvar_avgs_divs_cents(dsig_avgs_energy, data_sets_plt, plot=True, data_sets_colors=data_sets_colors,
                                  data_sets_labels=data_sets_labels, exclude_divs=exclude_divs, fit=False)
        plot_dvar_avgs_divs_cents(dsig_avgs_energy, data_sets_plt, plot=True, data_sets_colors=data_sets_colors,
                                  data_sets_labels=data_sets_labels, exclude_divs=exclude_divs, fit=True)
    plt.show()

    # df_onediv_fits, df_tp_fits, df_divs_fits = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    # for energy in energies_fit:
    #     df_fit = stat_vs_protons_cents(df, stat_plot, [div_plt], cents, energy, data_types_plt, all_sets_plt,
    #                                    plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
    #                                    data_sets_labels=data_sets_labels)
    #     df_onediv_fits = pd.concat([df_onediv_fits, df_fit])
    #
    #     df_tp_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy, data_types_plt, all_sets_plt,
    #                                       plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
    #                                       data_sets_labels=data_sets_labels)
    #     df_tp_fits = pd.concat([df_tp_fits, df_tp_fit])
    #
    #     df_fits = plot_protons_fits_divs_cents(df_tp_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
    #                                            data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)
    #     # print(f'{energy}GeV {df_fit}')
    #     df_divs_fits = pd.concat([df_divs_fits, df_fits])
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
    #                       data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=None,
    #                       ref_type=None, data_sets_energies_cmaps=data_sets_energies_cmaps)
    # plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
    #                       data_sets_labels=data_sets_labels, title=None, fit=True, cent_ref=cent_ref_df,
    #                       ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
    #                       data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
    #                       ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plt.show()


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
    v2_fit_out_dir = None
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
    df_raw_div = df[df['data_type'] == 'raw']
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
        if v2_fit_out_dir:
            flow_vs_v2(df, div, '15', v2_fit_out_dir)
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    # return
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_protons_fits_divs_flow(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                                  data_sets_labels=data_sets_labels)
    # if save_fits:
    #     df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    plt.show()


def plot_flow_k2():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_v2_ck2.csv'
    save_fits = False
    v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    v2_fit_out_dir = None
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = None  # 'flow_tprotons_fits.csv'
    df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_resample_res15_v207', 'flow_resample_res15_v205', 'flow_resample_res15_v202']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue']))
    data_sets_labels = dict(zip(data_sets_plt, ['v2=0.07', 'v2=0.05', 'v2=0.02']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    df = df[(df['data_type'] == 'raw') & ((df['stat'] == 'k2') | (df['stat'] == 'c2'))]
    df['val'] = df['val'] / (df['total_protons'] * df['divs'] / 360 * (1 - df['divs'] / 360))
    df['err'] = df['err'] / (df['total_protons'] * df['divs'] / 360 * (1 - df['divs'] / 360))

    stat_vs_protons(df, stat_plot, 120, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                    plot=True, fit=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        # if v2_fit_out_dir:
        #     flow_vs_v2(df, div, '15', v2_fit_out_dir)
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    # return
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_protons_fits_divs_flow(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                                  data_sets_labels=data_sets_labels)
    # if save_fits:
    #     df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
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

    pd.set_option('max_columns', None)

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_rp', 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                     'ampt_new_coal_epbins1_v2rpcor']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['Reaction Plane', 'Flow Included', 'Flow Corrected ep',
                                                'Flow Corrected rp']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}
    v2_vals = read_flow_values(v2_in_dir)  # Dictionary of {energy: {cent: v2}}
    v2_rp_vals = read_flow_values(v2_in_dir, 'v2_rp')  # Dictionary of {energy: {cent: v2}}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    # df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))  # Placeholder, need k2
    df_raw.loc[:, 'val'] = (df_raw['val'] ** 2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    # df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'val'] = (df_mix['val'] ** 2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))
    df_diff = df[(df['data_type'] == 'diff') & (df['stat'] == stat_plot)]
    p, tp = df_diff['divs'] / 360, df_diff['total_protons']
    df_diff.loc[:, 'val'] = df_diff['val'] / (tp * (tp - 1))
    df_diff.loc[:, 'err'] = df_diff['err'] / (tp * (tp - 1))

    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    dvar_vs_protons(df_sub, div_plt, cent_plt, [39], ['sub'], all_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)
    dvar_vs_protons_energies(df_sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['sub'], all_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.0014, 0.0011])

    protons_fits = []
    v2_apparents = {energy: [] for energy in energies_fit}
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets_plt, plot=False,
                                           avg=True)
        if div != 180:
            for energy in energies_fit:
                loc_ep = (protons_fits_div['name'] == 'ampt_new_coal_epbins1') & (protons_fits_div['energy'] == energy)
                loc_rp = (protons_fits_div['name'] == 'ampt_new_coal_rp') & (protons_fits_div['energy'] == energy)
                d_avg = protons_fits_div[loc_ep]['avg_meas'].iloc[0] - protons_fits_div[loc_rp]['avg_meas'].iloc[0]
                v2_apparent = np.sqrt(d_avg / flow_cor_coefs[div])
                v2_apparents[energy].append(v2_apparent)
                # print(f'{div}, {energy},  v2: {v2_apparent}')
        protons_fits_div = ampt_v2_closure_sub_dsigma(protons_fits_div, 'ampt_new_coal_epbins1',
                                                      'ampt_new_coal_epbins1_v2rpcor', v2_rp_vals,
                                                      flow_cor_coefs[div], cent_plt)
        protons_fits_div = ampt_v2_closure_sub_dsigma(protons_fits_div, 'ampt_new_coal_epbins1',
                                                      'ampt_new_coal_epbins1_v2cor',
                                                      v2_vals, flow_cor_coefs[div], cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    # plt.show()

    proton_fits_w_flow = protons_fits[protons_fits['name'] == 'ampt_new_coal_epbins1']
    v2_ep_e_vals, v2_ep_e_errs, v2_rp_e_vals, v2_rp_e_errs, v2_op_e_vals = [], [], [], [], []
    for energy in energies_fit:
        df_e = proton_fits_w_flow[proton_fits_w_flow['energy'] == energy]
        # print(df_e)
        divs, avgs = np.array(df_e['divs']), np.array(df_e['avg_meas'])
        # print(divs, avgs)
        coefs = np.array([flow_cor_coefs[div] for div in divs])
        v2 = v2_vals[energy][cent_plt].val
        v2_rp = v2_rp_vals[energy][cent_plt].val
        fig, ax = plt.subplots()
        v2s, chi_2s = [], []
        x_divs = np.linspace(60, 300, 1000)
        for v2_i in np.linspace(v2 * 0, v2 * 1.25, 50):
            cor_avgs = avgs - v2_divs(divs * np.pi / 180, v2_i)
            # cor_avgs = avgs - coefs * v2_i ** 2
            cor_avg_vals, cor_avg_errs = [x.val for x in cor_avgs], [x.err for x in cor_avgs]
            popt, pcov = cf(quad_180, divs, cor_avg_vals, sigma=cor_avg_errs, absolute_sigma=True)
            chi_2 = np.sum((cor_avg_vals - quad_180(divs, *popt)) ** 2 / cor_avg_errs)
            ax.scatter(divs, cor_avg_vals)
            ax.plot(x_divs, quad_180(x_divs, *popt))
            v2s.append(v2_i)
            chi_2s.append(chi_2)
        chi_2s, v2s = zip(*sorted(zip(chi_2s, v2s)))
        v2_op_e_vals.append(v2s[0])
        v2_ep_e_vals.append(v2)
        v2_ep_e_errs.append(v2_vals[energy][cent_plt].err)
        v2_rp_e_vals.append(v2_rp)
        v2_rp_e_errs.append(v2_rp_vals[energy][cent_plt].err)

        fig2, ax2 = plt.subplots()
        ax2.axhline(0, color='black')
        ax2.scatter(v2s, chi_2s)
        ax2.set_title(f'{energy}GeV')
        ax2.axvline(v2, color='green', ls='--', label='v2_ep')
        ax2.axvline(v2_rp, color='red', ls='--', label='v2_rp')
        ax2.legend()
        plt.show()
    fig3, ax3 = plt.subplots()
    ax3.grid()
    ax3.set_xlabel('Energy (GeV)')
    ax3.set_ylabel('v2')
    ax3.errorbar(energies_fit, v2_ep_e_vals, v2_ep_e_errs, ls='none', marker='o', color='green', label='v2_ep')
    ax3.errorbar(energies_fit, v2_rp_e_vals, v2_rp_e_errs, ls='none', marker='o', color='red', label='v2_rp')
    ax3.scatter(energies_fit, v2_op_e_vals, color='blue', label='v2_op')
    ax3.legend()
    # plt.show()

    # plt.show()
    # return

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    v2_divs_plt = np.setdiff1d(np.unique(df['divs']), exclude_divs + [180])
    energy_colors = {7: 'blue', 11: 'orange', 19: 'green', 27: 'red', 39: 'purple', 62: 'brown'}
    avg_v2aps = []
    for energy in energies_fit:
        fig_e, ax_e = plt.subplots(dpi=144)
        avg_v2ap = np.mean(v2_apparents[energy])
        avg_v2aps.append(avg_v2ap)
        ax_e.grid()
        ax.errorbar(v2_divs_plt, [x.val for x in v2_apparents[energy]], yerr=[x.err for x in v2_apparents[energy]],
                    marker='o', ls='none', alpha=0.6, color=energy_colors[energy], label=f'{energy}GeV')
        ax.axhspan(avg_v2ap.val + avg_v2ap.err, avg_v2ap.val - avg_v2ap.err, color=energy_colors[energy], alpha=0.3)
        ax.axhline(avg_v2ap.val, ls='--', color=energy_colors[energy])
        ax_e.errorbar(v2_divs_plt, [x.val for x in v2_apparents[energy]], yerr=[x.err for x in v2_apparents[energy]],
                      marker='o', ls='none', alpha=0.9, color=energy_colors[energy])
        ax_e.set_title(f'{energy} GeV')
        ax_e.axhspan(avg_v2ap.val + avg_v2ap.err, avg_v2ap.val - avg_v2ap.err, color=energy_colors[energy], alpha=0.5)
        ax_e.axhline(avg_v2ap.val, color=energy_colors[energy], ls='--', alpha=0.9)
        fig_e.canvas.manager.set_window_title(f'Apparent v2 {energy}GeV')
        fig_e.tight_layout()
    ax.legend()
    ax.set_ylabel('Apparent v2')
    fig.tight_layout()
    fig.canvas.manager.set_window_title('Apparent v2')

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    ax.errorbar(energies_fit, [x.val for x in avg_v2aps], [x.err for x in avg_v2aps], ls='none', marker='o',
                label='v2 Apparent')
    v2_plt_vals, v2_errs, v2_rp_plt_vals, v2_rp_errs = [], [], [], []
    for energy in energies_fit:
        v2_plt_vals.append(v2_vals[energy][cent_plt].val)
        v2_errs.append(v2_vals[energy][cent_plt].err)
        v2_rp_plt_vals.append(v2_rp_vals[energy][cent_plt].val)
        v2_rp_errs.append(v2_rp_vals[energy][cent_plt].err)
    ax.errorbar(energies_fit, v2_plt_vals, v2_errs, ls='none', marker='o', label='v2 Event Plane')
    ax.errorbar(energies_fit, v2_rp_plt_vals, v2_rp_errs, ls='none', marker='o', label='v2 Reaction Plane')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('v2')
    ax.set_title(f'V2 vs Energy Cent {cent_plt}')
    ax.legend()
    fig.tight_layout()
    fig.canvas.manager.set_window_title(f'V2 vs Energy Cent {cent_plt}')

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    v2_rp_div_app = []
    for v2_app, v2_rp_val, v2_rp_err in zip(avg_v2aps, v2_rp_plt_vals, v2_rp_errs):
        v2_rp_div_app.append(v2_app / Measure(v2_rp_val, v2_rp_err))
    ax.errorbar(energies_fit, [x.val for x in v2_rp_div_app], [x.err for x in v2_rp_div_app], ls='none', marker='o')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('v2')
    ax.set_title(f'V2 Reaction Plane / V2 Apparent vs Energy Cent {cent_plt}')
    fig.tight_layout()
    fig.canvas.manager.set_window_title(f'V2 Reaction Plane / V2 Apparent vs Energy Cent {cent_plt}')

    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_dvar_avgs_divs(protons_fits, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
    #                                  data_sets_labels=data_sets_labels, verbose=False)
    # df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_ampt_v2_closure_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_stats_ampt_flow_closure.csv'
    vs_in_dir = f'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/' \
                f'Ampt_rapid05_resample_norotate_epbins1_0/'
    df_path = base_path + df_name

    pd.set_option('max_columns', None)

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    energies_fit = [7, 11, 19, 27, 39, 62]

    data_sets_plt = ['ampt_new_coal_rp', 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                     'ampt_new_coal_epbins1_v2rpcor', 'ampt_new_coal_epbins1_flow_rpcor']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red', 'purple', 'orange']))
    data_sets_labels = dict(zip(data_sets_plt, ['Reaction Plane', 'Flow Included', 'Flow Corrected ep',
                                                'Flow Corrected rp', 'Flow Corrected rp 2-6']))

    all_sets_plt = data_sets_plt

    v_orders = [2]

    # Dictionary of {order: {energy: {cent: v2}}}
    v_ep_vals = {i: read_flow_values(vs_in_dir, f'v{i}') for i in [2]}
    v_rp_vals = {i: read_flow_values(vs_in_dir, f'v{i}_rp') for i in v_orders}

    # v2_vals = read_flow_values(vs_in_dir)  # Dictionary of {energy: {cent: v2}}
    # v2_rp_vals = read_flow_values(vs_in_dir, 'v2_rp')  # Dictionary of {energy: {cent: v2}}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))  # Placeholder, need k2
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    # df_raw.loc[:, 'err'] = abs(2 * df_raw['val'] * df_raw['err'] / (tp * (tp - 1)))
    # df_raw.loc[:, 'val'] = (df_raw['val']**2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))
    # df_mix.loc[:, 'err'] = abs(2 * df_mix['val'] * df_mix['err'] / (tp * (tp - 1)))
    # df_mix.loc[:, 'val'] = (df_mix['val']**2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_diff = df[(df['data_type'] == 'diff') & (df['stat'] == stat_plot)]
    p, tp = df_diff['divs'] / 360, df_diff['total_protons']
    df_diff.loc[:, 'val'] = df_diff['val'] / (tp * (tp - 1))
    df_diff.loc[:, 'err'] = df_diff['err'] / (tp * (tp - 1))
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'

    dvar_vs_protons(df_sub, div_plt, cent_plt, [39], ['sub'], all_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'sub', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)
    dvar_vs_protons_energies(df_sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['sub'], all_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.0014, 0.0011])

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets_plt, plot=False,
                                           avg=True)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_v2rpcor', {2: v_rp_vals[2]},
                                                div, cent_plt)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_v2cor',
                                                v_ep_vals, div, cent_plt)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_flow_rpcor', v_rp_vals,
                                                div, cent_plt)
        # protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
        #                                                 'ampt_new_coal_epbins1_flow_cor',
        #                                                 v_ep_vals, div, cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    plot_dvar_avgs_divs(protons_fits, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels, alpha=0.6)

    plt.show()


def plot_flow_v2_closure():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_flow.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
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
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res15_v207', 'flow_resample_res15_v207_v2sub']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red']))
    data_sets_labels = dict(zip(data_sets_plt, ['Flow Excluded', 'Flow Included', 'Flow Subtracted']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
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
        print(protons_fits_div)
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div, 'flow_resample_res15_v207',
                                               'flow_resample_res15_v207_v2sub', 0.07, flow_cor_coefs[div])
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
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_flow_v2_closure_raw():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_flow.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
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
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    pd.set_option('max_columns', None)

    data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res15_v207', 'flow_resample_res15_v207_v2sub']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red']))
    data_sets_labels = dict(zip(data_sets_plt, ['Flow Excluded', 'Flow Included', 'Flow Subtracted']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[df['data_type'] == 'raw']
    p = df_raw['divs'] / 360
    df_raw.loc[:, 'val'] = df_raw['val'] / np.sqrt((df_raw['total_protons'] * p * (1 - p)))
    df_mix = df[df['data_type'] == 'mix']
    p = df_mix['divs'] / 360
    df_mix.loc[:, 'val'] = df_mix['val'] / np.sqrt((df_mix['total_protons'] * p * (1 - p)))
    print(df)
    print(df_raw)

    stat_vs_protons(df_raw, stat_plot, div_plt, cent_plt, [62], ['raw'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons(df_mix, stat_plot, div_plt, cent_plt, [62], ['mix'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    # stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
    #                          plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = stat_vs_protons(df_raw, stat_plot, div, cent_plt, energies_fit, ['raw'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_mix = stat_vs_protons(df_mix, stat_plot, div, cent_plt, energies_fit, ['mix'],
                                               all_sets_plt, plot=False, fit=True)
        # print(protons_fits_div_raw)
        for index, row in protons_fits_div_mix.iterrows():
            if data_sets_labels[row['name']] == 'Flow Included':
                continue
            loc_index = (protons_fits_div_raw.name == row['name']) & (protons_fits_div_raw.energy == row['energy'])
            raw_sub_meas = (protons_fits_div_raw.loc[loc_index, 'slope_meas'] - row['slope_meas']).iloc[0]
            # print(row, raw_sub_meas)
            protons_fits_div_raw.loc[loc_index, 'slope_meas'] = raw_sub_meas
            protons_fits_div_raw.loc[loc_index, 'slope'] = raw_sub_meas.val
            protons_fits_div_raw.loc[loc_index, 'slope_err'] = raw_sub_meas.err

        # print(protons_fits_div_raw)
        # input()
        # print(protons_fits_div_mix)
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div_raw, 'flow_resample_res15_v207',
                                               'flow_resample_res15_v207_v2sub', 0.07, flow_cor_coefs[div], cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels, verbose=False)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_flow_eff_test():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_var_cent8_2source_closure_tests.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_eff_res15_v207', 'flow_res15_v207']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'blue']))
    data_sets_labels = dict(zip(data_sets_plt, ['flow_eff_test', 'flow_test']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[df['data_type'] == 'raw']
    p = df_raw['divs'] / 360
    df_raw.loc[:, 'val'] = df_raw['val'] / (df_raw['total_protons'] * p * (1 - p))
    df_raw.loc[:, 'err'] = df_raw['err'] / (df_raw['total_protons'] * p * (1 - p))
    df_mix = df[df['data_type'] == 'mix']
    p = df_mix['divs'] / 360
    df_mix.loc[:, 'val'] = df_mix['val'] / (df_mix['total_protons'] * p * (1 - p))
    df_mix.loc[:, 'err'] = df_mix['err'] / (df_mix['total_protons'] * p * (1 - p))
    print(df)
    print(df_raw)

    stat_vs_protons(df_raw, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels,
                    y_ranges={'k2': [0.88, 1.09]})

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = stat_vs_protons(df_raw, stat_plot, div, cent_plt, energies_fit, ['raw'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        protons_fits_div_mix = stat_vs_protons(df_mix, stat_plot, div, cent_plt, energies_fit, ['mix'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        protons_fits.append(protons_fits_div_mix)

        sub_meas = protons_fits_div_raw['slope_meas'] - protons_fits_div_mix['slope_meas']
        protons_fits_div_sub = protons_fits_div_raw.copy()
        protons_fits_div_sub.loc[:, 'slope'] = [x.val for x in sub_meas]
        protons_fits_div_sub.loc[:, 'slope_err'] = [x.err for x in sub_meas]
        protons_fits_div_sub.loc[:, 'slope_meas'] = sub_meas
        protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_mix['name'] + '_sub'
        protons_fits.append(protons_fits_div_sub)

        protons_fits_div_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, ['divide'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_div.loc[:, 'name'] = protons_fits_div_div['name'] + '_div'
        protons_fits.append(protons_fits_div_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_raw', '_mix', '_div', '_sub']]
        colors = dict(zip(data_sets, ['blue', 'green', 'red', 'purple']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Raw', ' Mix', ' Div', ' Sub']]))
        plot_protons_fits_divs(protons_fits, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels)
    data_sets = ['flow_eff_res15_v207_raw', 'flow_eff_res15_v207_div', 'flow_res15_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Flow+Efficiency Raw', 'Flow+Efficiency Raw/Mix', 'Flow Raw']))
    plot_protons_fits_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels)

    plt.show()


def plot_anticl_flow_closure_test():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_flow_anticl_convo_test.csv'
    df_name = 'binom_slice_var_cent8_2source_closure_tests.csv'
    df_name = 'binom_slice_var_cent8_v2_anticl_closure.csv'
    # save_fits = False
    # v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    # v2_fit_out_dir = None
    # fits_out_base = 'Base_Zero_Fits/'
    # df_tproton_fits_name = None  # 'flow_tprotons_fits.csv'
    # df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    # data_sets_plt = ['flow_resample_res15_v207', 'anticlmulti_resample_s05_a05',
    #                  'anticlflow_resample_res15_v207_s05_a05']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['flow', 'anticlustering', 'both']))

    df = pd.read_csv(df_path)
    df = df.dropna()

    # Delete this block if I append v2 value to end of name in calc script
    # df_no_flow = df[df['name'].str.contains('anticlmulti')]
    # df_flow = df[df['name'].str.contains('anticlflow')]
    # df_flow.loc[:, 'name'] = df_flow['name'] + '_v207'
    # df = pd.concat([df_flow, df_no_flow], ignore_index=True)

    all_sets = pd.unique(df['name'])
    print(all_sets)

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))

    dvar_vs_protons(df_raw, div_plt, cent_plt, [62], ['raw'], all_sets, plot=True, avg=True)
    dvar_vs_protons(df_mix, div_plt, cent_plt, [62], ['mix'], all_sets, plot=True, avg=True)
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    dvar_vs_protons(df_sub, div_plt, cent_plt, [62], ['sub'], all_sets, plot=True, avg=True)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub], ignore_index=True), div_plt, cent_plt, [62],
                    ['raw', 'mix', 'sub'], all_sets, plot=True, avg=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets, plot=False,
                                               avg=True)
        # protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        # protons_fits_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets, plot=False,
        #                                        avg=True)
        # protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        # protons_fits.append(protons_fits_div_mix)
        #
        # protons_fits_div_sub = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets, plot=False,
        #                                        avg=True)
        # protons_fits.append(protons_fits_div_sub.copy())
        # protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_sub['name'] + '_sub'
        # protons_fits.append(protons_fits_div_sub)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    # plot_protons_fits_divs_flow(protons_fits, all_sets)

    for flow_set in [set_name for set_name in all_sets if 'anticlflow' in set_name and '_eff_' not in set_name]:
        no_flow_set = flow_set[:flow_set.find('_v2')].replace('anticlflow', 'anticlmulti')

        print(pd.unique(protons_fits['name']))
        print(flow_set)
        anticl_plus_v2 = protons_fits[protons_fits['name'] == flow_set]
        fits = [anticl_plus_v2.copy(), protons_fits[protons_fits['name'] == no_flow_set]]
        divs, avgs = np.array(anticl_plus_v2['divs']), np.array(anticl_plus_v2['avg_meas'])
        print(divs, avgs)
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        v2 = 0.07
        v2s, chi_2s = [], []
        x_divs = np.linspace(60, 300, 1000)
        for v2_i in np.linspace(0, v2 * 1.25, 50):
            cor_avgs = avgs - v2_divs(divs * np.pi / 180, v2_i)
            # print(v2_divs(divs * np.pi / 180, v2_i))
            # print(cor_slopes)
            cor_avg_vals, cor_avg_errs = [x.val for x in cor_avgs], [x.err for x in cor_avgs]
            ax2.errorbar(divs, cor_avg_vals, yerr=cor_avg_errs, ls='none', alpha=0.4, marker='o')
            popt, pcov = cf(quad_180, divs, cor_avg_vals, sigma=cor_avg_errs, absolute_sigma=True)
            chi_2 = np.sum((cor_avg_vals - quad_180(divs, *popt)) ** 2 / cor_avg_errs)
            ax.scatter(divs, cor_avg_vals)
            ax.plot(x_divs, quad_180(x_divs, *popt))
            v2s.append(v2_i)
            chi_2s.append(chi_2)
        fig2, ax2 = plt.subplots()
        ax2.axhline(0, color='black')
        ax2.scatter(v2s, chi_2s)
        chi_2s, v2s = zip(*sorted(zip(chi_2s, v2s)))
        # anticl_plus_v2['name'] = f'corrected_v2{str(v2s[0])[2:]}_'
        anticl_plus_v2.loc[:, 'name'] = 'corrected_v2'
        # new_slope = anticl_plus_v2['slope'] - v2_divs(anticl_plus_v2['divs'] * np.pi / 180, v2s[0])
        # anticl_plus_v2 = anticl_plus_v2.assign(slope=new_slope)
        anticl_plus_v2.loc[:, 'avg'] = anticl_plus_v2['avg'] - v2_divs(anticl_plus_v2['divs'] * np.pi / 180, v2s[0])
        print(anticl_plus_v2)
        print(v2s[0], chi_2s[0])
        fits = pd.concat([*fits, anticl_plus_v2], ignore_index=True)
        # data_sets_colors.update({'corrected_v2': 'olive'})
        plot_dsigma_fits_divs_flow(fits, [flow_set, no_flow_set] + ['corrected_v2'])

    plt.show()


def plot_efficiency_closure_tests():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_var_cent8_efficiency_closure.csv'
    df2_name = 'binom_slice_var_cent8_v2_anticlindep_closure.csv'
    # df_name = 'binom_slice_var_cent8_simpleclust_test.csv'
    df_path = base_path + df_name
    df2_path = base_path + df2_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    samples = 72  # For title purposes only

    df = pd.read_csv(df_path)
    df = pd.concat([df, pd.read_csv(df2_path)], ignore_index=True)
    df = df.dropna()

    # df = df[df['name'].str.contains('simpleclust_eff')]

    all_sets = pd.unique(df['name'])
    print(all_sets)

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))

    dvar_vs_protons(df_raw, div_plt, cent_plt, [62], ['raw'], all_sets, plot=True, avg=True)
    dvar_vs_protons(df_mix, div_plt, cent_plt, [62], ['mix'], all_sets, plot=True, avg=True)
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    dvar_vs_protons(df_sub, div_plt, cent_plt, [62], ['sub'], all_sets, plot=True, avg=True)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub], ignore_index=True), div_plt, cent_plt, [62],
                    ['raw', 'mix', 'sub'], all_sets, plot=True, avg=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets, plot=False,
                                               avg=True)
        protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        protons_fits_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets, plot=False,
                                               avg=True)
        protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        protons_fits.append(protons_fits_div_mix)

        protons_fits_div_sub = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets, plot=False,
                                               avg=True)
        protons_fits.append(protons_fits_div_sub.copy())
        protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_sub['name'] + '_sub'
        protons_fits.append(protons_fits_div_sub)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    for data_set in all_sets:
        data_sets = [data_set + x for x in ['_raw', '_mix', '_sub']]
        colors = dict(zip(data_sets, ['blue', 'green', 'red']))
        labels = dict(zip(data_sets, [data_set + x for x in [' Raw', ' Mix', ' Sub']]))
        plot_dvar_avgs_divs(protons_fits, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels,
                            plt_energies=False, title=data_set)

    data_sets = ['flow_eff_res15_v207_raw', 'flow_eff_res15_v207_sub', 'flow_res15_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Flow + Efficiency Raw', 'Flow + Efficiency Raw-Mix', 'Flow Raw']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plt_energies=False, alpha=0.6, title='Flow v2=0.07 Efficiency Correction')

    data_sets = ['simpleclust_eff_raw', 'simpleclust_eff_sub', 'simpleclust_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Simple Clustering + Efficiency Raw', 'Simple Clustering + Efficiency Raw-Mix',
                                  'Simple Clustering Raw']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plt_energies=False, alpha=0.6, title='Simple Clustering (s05, a2) Efficiency Correction')

    for s, a in [('1', '01'), ('08', '01'), ('05', '01'), ('01', '01')]:
        data_sets = [f'anticlflow_eff_s{s}_a{a}_raw', f'anticlflow_eff_s{s}_a{a}_sub', f'anticlmulti_s{s}_a{a}_raw']
        colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
        labels = dict(zip(data_sets, ['Anti-Clustering + Efficiency Raw', 'Anti-Clustering + Efficiency Raw-Mix',
                                      'Anti-Clustering Raw']))
        plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                            plt_energies=False, alpha=0.6,
                            title=f'Anti-Clustering (s{s}, a{a}) Efficiency Correction')

    for s, a in [('1', '01'), ('08', '01'), ('05', '01'), ('01', '01')]:
        data_sets = [f'anticlflowindep_eff_s{s}_a{a}_raw', f'anticlflowindep_eff_s{s}_a{a}_sub',
                     f'anticlmulti_s{s}_a{a}_raw']
        colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
        labels = dict(zip(data_sets, ['Anti-Clustering (Independent) + Efficiency Raw',
                                      'Anti-Clustering (Independent) + Efficiency Raw-Mix', 'Anti-Clustering Raw']))
        plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                            plt_energies=False, alpha=0.6,
                            title=f'Anti-Clustering (Independent) (s{s}, a{a}) Efficiency Correction')

    data_sets = [f'anticlflowindep_eff_s05_a01_v207_raw', f'anticlflowindep_eff_s05_a01_v207_sub',
                 f'anticlflowindep_s05_a01_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Anti-Clustering (Independent) + Efficiency + Flow Raw',
                                  'Anti-Clustering (Independent) + Efficiency + Flow Raw-Mix',
                                  'Anti-Clustering + Flow Raw']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plt_energies=False, alpha=0.6,
                        title=f'Anti-Clustering (Independent) (s{s}, a{a}) + Flow v2=0.07 Efficiency Correction')

    # data_sets = ['anticlflow_eff_s1_a01_v207_raw', 'anticlflow_eff_s1_a01_v207_sub', '']
    # colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    # labels = dict(zip(data_sets, ['Anti-Clustering + Flow + Efficiency Raw',
    #                               'Anti-Clustering + Flow + Efficiency Raw-Mix', 'Anti-Clustering + Flow Raw']))
    # plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
    #                     plt_energies=False, alpha=0.6, title='Flow v2=0.07 Efficiency Correction')

    plt.show()


if __name__ == '__main__':
    main()
