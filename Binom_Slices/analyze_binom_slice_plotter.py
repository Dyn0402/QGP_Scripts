#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 17 4:27 PM 2022
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slice_plotter

@author: Dylan Neff, Dyn04
"""

from analyze_binom_slices import *


def main():
    plot_sims()
    # plot_star_model()
    # plot_star_model_onediv()
    print('donzo')


def plot_star_model():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
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
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_sims():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_sim_test.csv'
    df_path = base_path + df_name
    sim_sets = []

    amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    spreads = ['03', '04', '05', '06', '07', '08', '09', '1', '11', '12']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
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

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
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
    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    plot_sigma_fits_interp(sigma_fits)

    plt.show()


def plot_star_model_onediv():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
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
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV 1fm^3',
                                                'MUSIC+FIST EV 3.42fm^3']))

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


if __name__ == '__main__':
    main()
