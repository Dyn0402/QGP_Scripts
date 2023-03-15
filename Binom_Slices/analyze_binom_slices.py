#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 18 5:16 PM 2021
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slices

@author: Dylan Neff, dylan
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
from mpl_toolkits import mplot3d
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit as cf
from scipy import odr
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import scipy.interpolate as interp
from scipy.interpolate import RectBivariateSpline
from scipy.stats import norm

from calc_binom_slices import get_name_amp_spread, get_name_amp_spread_pm

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats
from Measure import Measure


def main():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    threads = 15
    # base_path = 'F:/Research/Results/Azimuth_Analysis/'
    base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    # df_name = 'binom_slice_sds_cent8.csv'
    # df_name = 'binom_slice_stats_cent8_no_sim.csv'
    # df_name = 'binom_slice_stats_cent8_ampt_eff.csv'
    # df_name = 'binom_slice_stats_cent8_no_sim.csv'
    df_name = 'binom_slice_stats_cent8_sim_test.csv'
    # df_name = 'binom_slice_stats_cent8_ex_cl_ant_sim.csv'
    chi_df_name = 'chi_df_ampt_neks_cent8.csv'
    df_path = base_path + df_name
    # df_path = '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'F:/Research/Results/Azimuth_Analysis/binom_slice_cent_sds_df.csv'
    # df_path = 'F:/Transfer/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'C:/Users/Dyn04/Desktop/binom_slice_df.csv'
    # df_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/binom_slice_cent_sds_df.csv'
    sim_sets = []
    # sim_sets = ['sim_aclmul_amp02_spread1', 'sim_aclmul_amp2_spread1', 'sim_clmul_amp02_spread1',
    #             'sim_clmul_amp2_spread1']
    amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    spreads = ['05', '1']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    # sim_amp_pairs = [('2', '1'), ('5', '1'), ('2', '05'), ('5', '05')]
    # # sim_amp_pairs = []  # [('5', '35'), ('02', '05'), ('015', '1')]
    # for amp, spread in sim_amp_pairs:
    #     sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
    #
    # sim_amp_pairs = [('1', '1'), ('5', '1'), ('1', '05'), ('5', '05')]
    # # sim_amp_pairs = []  # [('5', '35'), ('02', '05'), ('015', '1')]
    # for amp, spread in sim_amp_pairs:
    #     sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')

    # spreads.extend(['5', '6'])
    # for amp in amps:
    #     for spread in spreads:
    #         sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}_test')
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.25), 'skewness': (0, 1.25), 'kurtosis': (-3, 2),
                'non-excess kurtosis': (0.95, 1.025)}
    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    total_protons_plt = 20
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    # energies_fit = [7, 11, 19, 27, 39, 62]  # , 11, 19, 27, 39, 62]
    energies_fit = ['sim', 7]
    energy_plt = 62
    data_types_plt = ['divide']
    data_type_plt = 'divide'
    samples = 72  # For title purposes only
    # data_sets_plt = ['ampt_eff1_resample_def', 'ampt_eff2_resample_def', 'ampt_eff3_resample_def']  # ['ampt_resample_def', 'ampt_old_resample_def', 'bes_resample_def']  # ['ampt_resample_def']
    # data_sets_plt = ['ampt_old_eff1_resample_def', 'ampt_old_eff2_resample_def', 'ampt_old_eff3_resample_def']

    # data_sets_plt = ['bes_resample_def', 'ampt_baryon_first_resample_def', 'ampt_meson_first_resample_def',
    #                  'ampt_new_coal_resample_def', 'ampt_baryon_first_fix_resample_def',
    #                  'cf_resample_def', 'cfev_resample_def', 'cfevb342_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'green', 'salmon', 'cyan', 'blue', 'purple', 'olive']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT_Baryon_First', 'AMPT_Meson_First', 'AMPT_New_Coalescence',
    #                                             'AMPT_Baryon_Fix', 'MUSIC+FIST', 'MUSIC+FIST+EV', 'MUSIC+FIST+EVb342']))

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'ampt_baryon_first_fix_resample_def',
    #                  'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'salmon', 'green', 'olive']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT_New_Coalescence', 'AMPT_Baryon_First',
    #                                             'MUSIC+FIST', 'MUSIC+FIST+EV']))

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV']))

    # data_sets_plt = ['cf_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['MUSIC+FIST']))

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

    # data_sets_plt = ['bes_resample_def', 'bes_resample_def_alg3', 'bes_single']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR Stochastic', 'STAR Even Space', 'STAR Single']))

    # data_sets_plt = ['ampt_meson_first_resample_def', 'cf_resample_def', 'cfev_resample_def', 'cfevb342_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['red', 'blue', 'green', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['AMPT', 'MUSIC+FIST', 'MUSIC+FIST+EV', 'MUSIC+FIST+EVb342']))

    # data_sets_plt = ['bes_resample_def',
    #                  'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def', 'cfevb342_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'salmon', 'blue', 'purple', 'olive']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT_New_Coalescence',
    #                                             'MUSIC+FIST', 'MUSIC+FIST+EV', 'MUSIC+FIST+EVb342']))

    # data_sets_labels = dict(zip(sim_sets, ['Repulsive Amp=0.2 Spread=1', 'Repulsive Amp=0.5 Spread=1',
    #                                        'Repulsive Amp=0.2 Spread=0.5', 'Repulsive Amp=0.5 Spread=0.5',
    #                                        'Attractive Amp=0.1 Spread=1', 'Attractive Amp=0.5 Spread=1',
    #                                        'Attractive Amp=0.1 Spread=0.5', 'Attractive Amp=0.5 Spread=0.5']))
    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))
    # df_data = df[~df['name'].str.contains('sim_')]
    # df_sim = df[df['name'].str.contains('sim_')]
    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    # df = df[df['spread'] != 1.2]

    # fit_pars = stat_vs_divs(df, stat_plot, total_protons_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt,
    #                         exclude_divs, plot=True, fit=True)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=True,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    # plt.show()
    # return

    # protons_fits = stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt,
    #                                all_sets_plt, plot=False, fit=True)
    # protons_fits = stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt,
    #                                all_sets_plt, plot=False, fit=True)
    # print(protons_fits)
    # plot_protons_fits_vs_energy(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, title=f'0-5% Centrality, {div_plt}° Partitions, '
    #                                                                      f'{samples} Samples per Event')
    # plt.show()
    # return
    # plot_protons_fits_vs_amp(protons_fits, all_sets_plt, data_sets_colors, data_sets_labels)
    #
    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    #
    # stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
    #                          plot=True, fit=False, plot_fit=False, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels)
    #
    # stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
    #                          plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], ['raw', 'mix'], ['bes_resample_def'], plot=True, fit=False)
    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, ['bes_resample_def'], plot=True, fit=False)
    # plt.show()
    # return

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
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    plot_slope_div_fits_simpars(df_fits)

    plt.show()
    return

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt, plot=True, fit=True)
    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, ['raw', 'mix'], all_sets_plt, plot=True, fit=False)
    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    for energy in [7, 11, 19, 27, 39, 62]:
        stat_vs_protons_divs(df, stat_plot, [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356], cent_plt, [energy],
                             data_types_plt, all_sets_plt, plot=True, fit=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels)
    plt.show()
    return

    # chi_res = []
    # chi2_list_all = []
    # for energy in energies_fit:
    #     chi2_sets = []
    #     jobs = []
    #     print(f'Energy {energy}GeV')
    #     for data_type in data_types_plt:
    #         for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #             jobs.append((df, stat_plot, div, cent_plt, energy, data_type, data_sets_plt))
    #     with Pool(threads) as pool:
    #         for chi2_set in tqdm.tqdm(pool.istarmap(chi2_vs_protons, jobs), total=len(jobs)):
    #             chi2_sets.extend(chi2_set)
    #     # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     #     print(f'{energy}GeV {div} bin width')
    #     #     chi2_sets.extend(chi2_vs_protons(df, stat_plot, div, cent_plt, energy, data_type_plt, data_sets_plt))
    #     data_sim_sets_plt, chi_res_e, chi_list = plot_chi2_protons(pd.DataFrame(chi2_sets), energy, n_sims_plt=4)
    #     chi2_list_all.extend(chi_list)
    #     chi_res.extend(chi_res_e)
    #
    # chi2_sets_df = pd.DataFrame(chi2_list_all)
    # chi2_sets_df.to_csv(base_path + chi_df_name, index=False)

    # chi_res = []
    # chi2_sets_df = pd.read_csv(base_path + chi_df_name)
    # for energy in np.unique(chi2_sets_df['energy']):
    #     data_sim_sets_plt, chi_res_e = plot_chi2_protons(chi2_sets_df, energy, n_sims_plt=4)
    #     chi_res.extend(chi_res_e)

    # plot_chi_sim_pars(pd.DataFrame(chi_res))
    # for data_set, sim_sets_plt in data_sim_sets_plt.items():
    #     sets = [data_set] + list(sim_sets_plt)
    #     print(sets)
    #     for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):
    #         stat_vs_protons(df, stat_plot, div, cent_plt, energies_plt, data_types_plt, sets, plot=True, fit=False)

    # stat_vs_divs(df, stat_plot, total_protons_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt,
    #              exclude_divs, plot=True, fit=True)
    # df_pars = stat_vs_divs(df, stat_plot, total_protons_plt, cent_plt, energies_plt, data_types_plt, all_sets_plt,
    #                        exclude_divs, plot=False, fit=True)
    # df_fits = interp_vs_div(df_pars, stat_plot, quad_par_line=line_yint1, r_amp_line=line_yint0, plot=True)
    # for div_plt in [240]:
    # plot_vs_protons(df, stats_plot, div_plt, cent_plt, data_type_plt, y_ranges, data_sets_plt)
    # for total_protons_plt in [20]:
    #     df_pars = fit_vs_divs(df, stats_plot, total_protons_plt, cent_plt, energies_plt, data_type_plt, y_ranges,
    #                           all_sets_plt, exclude_divs, False)
    #     df_fits = interp_vs_div(df_pars, y_ranges, quad_par_line=line_yint1, r_amp_line=line_yint0, plot=True)
    #     # plot_fits(df_fits)
    # plt.show()

    print('donzo')


def line(x, a, b):
    return a * x + b


def line_yint1(x, a):
    return a * (x - 1) + 1


def line_yint0(x, a):
    return a * x


def quad_180(x, a, c):
    return a * (x - 180) ** 2 + c


def quad_180_rad(x, a, c):
    return a * (x - np.pi) ** 2 + c


def quad_180_zparam(x, z, c):
    return -c * ((x - 180) / z) ** 2 + c


def double_quad_180_zparam(x, z1, c1, z2, c2):
    return quad_180_zparam(x, z1, c1) * quad_180_zparam(x, z2, c2)


def inv_sqrtx(x, a, c):
    return a / np.sqrt(x) + c


def inv_sqrtx_odr(pars, x):
    return pars[0] / np.sqrt(x) + pars[1]


def inv_sqrtxlin_odr(pars, x):
    return pars[0] / np.sqrt(x) + pars[1] + pars[2] * x


def inv_invxpow_odr(pars, x):
    return pars[0] / np.power(x, pars[2]) + pars[1]


def inv_invxpow_noc_odr(pars, x):
    return pars[0] / np.power(x, pars[1])


def inv_invx_odr(pars, x):
    return pars[0] / x + pars[1]


def quad(x, a, b, c):
    return a * x ** 2 + b * x + c


def x2(x, a):
    return a * x ** 2


def v2_divs_binom_norm(w, v2):
    return v2 ** 2 * np.sin(w) ** 2 / (np.pi * w * (1 - w / (2 * np.pi)))


def v2_divs(w, v2):
    return v2 ** 2 * np.sin(w) ** 2 / (2 * np.pi**2)


def stat_vs_protons(df, stat, div, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False, fit=False,
                    hist=False, data_sets_colors=None, data_sets_labels=None, star_prelim=False):
    data = []
    for data_type in data_types:
        for data_set in data_sets_plt:
            for energy in energies:
                df_pre = df
                if 'data_type' in df_pre:
                    df_pre = df_pre[df_pre['data_type'] == data_type]
                if 'cent' in df_pre:
                    df_pre = df_pre[df_pre['cent'] == cent]
                if 'stat' in df_pre:
                    df_pre = df_pre[df_pre['stat'] == stat]
                df_set = df_pre[(df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                if len(df_set) == 0:
                    continue
                if data_sets_labels is None:
                    if energy == 'sim':
                        lab = f'{data_set}_{data_type}'
                    else:
                        lab = f'{data_set}_{data_type}_{energy}GeV'
                else:
                    lab = data_sets_labels[data_set]
                df_set = df_set.sort_values(by=['total_protons'])
                if '_clmultipm_' in data_set:
                    ap, am = df_set['ampplus'].iloc[0], df_set['ampminus'].iloc[0]
                    sp, sm = df_set['spreadplus'].iloc[0], df_set['spreadminus'].iloc[0]
                    amp, spread = 0, 0
                else:
                    amp, spread = df_set['amp'].iloc[0], df_set['spread'].iloc[0]
                    ap, am, sp, sm = 0, 0, 0, 0
                data.append((df_set, lab, data_set, amp, spread, ap, am, sp, sm, energy))

    if plot:
        fig, ax = plt.subplots(figsize=(6.66, 5), dpi=144)
        fig.canvas.manager.set_window_title(f'Raw Div Mix {stat.title()} vs Total Protons {energy}GeV {div}°')
        ax.set_title(f'{energy} GeV, 0-5% Centrality, {div}° Partitions, 72 Samples per Event')
        ax.set_xlabel('Total Protons in Event')
        ax.set_ylabel(f'Raw / Mix {stat.title()}')
        ax.axhline(1, ls='--', color='black')
        if y_ranges:
            ax.set_ylim(y_ranges[stat])
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))

    fits = []
    for i, (df, lab, data_set, amp, spread, ap, am, sp, sm, energy) in enumerate(data):
        zo = len(data) - i + 4
        if plot:
            if data_sets_colors is None:
                c = next(color)
            else:
                c = data_sets_colors[data_set]
            if 'sim_' in data_set:
                ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'], label=lab, color=c,
                                alpha=0.4)
            else:
                ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                            marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                if 'sys' in df:
                    ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.4, zorder=zo)

        if fit and len(df) > 1:
            popt, pcov = cf(line_yint1, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
            fits.append({'name': data_set, 'energy': energy, 'divs': div, 'amp': amp, 'spread': spread, 'ap': ap,
                         'am': am, 'sp': sp, 'sm': sm, 'slope': popt[0], 'slope_err': np.sqrt(np.diag(pcov))[0],
                         'slope_meas': Measure(popt[0], np.sqrt(np.diag(pcov))[0]),
                         'int': 1, 'int_err': 0})
            if plot:
                ax.plot(df['total_protons'], line_yint1(df['total_protons'], *popt), ls='--', color=c)
                if hist:
                    sigs = (df['val'] - line_yint1(df['total_protons'], *popt)) / df['err']
                    fig_hist, ax_hist = plt.subplots()
                    ax_hist.set_title(f'{lab}')
                    ax_hist.set_xlabel('Standard Deviations from Linear Fit')
                    sns.histplot(sigs, stat='density', kde=True)
                    x_norm = np.linspace(min(sigs), max(sigs), 1000)
                    ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                    ax_hist.legend()
                    fig_hist.tight_layout()

    if plot:
        if star_prelim:
            ax.text(10, 0.955, 'STAR Preliminary', fontsize=15)
            eta_line = r'|$\eta$| < 1'
            pt_line = r'0.4 < $p_T$ < 2.0 GeV'
            ax.text(10, 0.962, f'Au+Au\n{eta_line}\n{pt_line}')
            ax.text(31.6, 1.041, 'statistical uncertainty only', fontsize=8)
            ax.legend(loc='upper center', framealpha=1.0).set_zorder(10)
        else:
            ax.legend(loc='upper left')
        fig.tight_layout()

    return pd.DataFrame(fits)


def stat_binom_vs_protons(df, stat, div, cent, energy, data_types, data_set_plt, y_ranges=None):
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}
    data = []
    tp_min, tp_max = 1000, -1
    for data_type in data_types:
        df_pre = df
        if 'data_type' in df_pre:
            df_pre = df_pre[df_pre['data_type'] == data_type]
        if 'cent' in df_pre:
            df_pre = df_pre[df_pre['cent'] == cent]
        if 'stat' in df_pre:
            df_pre = df_pre[df_pre['stat'] == stat]
        df_set = df_pre[(df_pre['name'] == data_set_plt) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
        if data_type == 'raw':
            lab = 'Single Event'
        elif data_type == 'mix':
            lab = 'Mixed Event'
        else:
            lab = '?'
        df_set = df_set.sort_values(by=['total_protons'])
        data.append((df_set, lab, data_type))
        set_tp_min, set_tp_max = df_set['total_protons'].min(), df_set['total_protons'].max()
        tp_min = tp_min if set_tp_min > tp_min else set_tp_min
        tp_max = tp_max if set_tp_max < tp_max else set_tp_max

    fig, ax = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig.canvas.manager.set_window_title(f'Binomial Comparison {stat.title()} vs Total Protons {energy}GeV {div}°')
    ax.set_title(f'{data_set_plt}\n{energy} GeV, {cent_map[cent]} Centrality, {div}° Partitions, 72 Samples per Event')
    ax.set_xlabel('Total Protons in Event')
    ax.set_ylabel(f'{stat.title()}')
    if y_ranges:
        ax.set_ylim(y_ranges[stat])

    xs = np.linspace(tp_min, tp_max, 1000)
    p = div / 360
    if stat in ['variance', 'k2']:
        binom_ys = xs * p * (1 - p)
    else:
        print(f'Stat {stat} not implemented!')
        binom_ys = xs * 0
    ax.plot(xs, binom_ys, color='red', label='Binomial')

    for df, lab, data_type in data:
        c = 'black'
        if data_type == 'raw':
            c = 'blue'
        elif data_type == 'mix':
            c = 'green'

        ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab, marker='o', ls='', color=c, alpha=0.7)

    ax.legend()
    fig.tight_layout()


def dvar_vs_protons(df, div, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False, avg=False,
                    hist=False, data_sets_colors=None, data_sets_labels=None, star_prelim=False):
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}
    data = []
    for data_type in data_types:
        for data_set in data_sets_plt:
            for energy in energies:
                df_pre = df
                if 'data_type' in df_pre:
                    df_pre = df_pre[df_pre['data_type'] == data_type]
                if 'cent' in df_pre:
                    df_pre = df_pre[df_pre['cent'] == cent]
                df_set = df_pre[(df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                if len(df_set) == 0:
                    continue
                if data_sets_labels is None:
                    if energy == 'sim':
                        lab = f'{data_set}_{data_type}'
                    else:
                        lab = f'{data_set}_{data_type}_{energy}GeV'
                else:
                    lab = data_sets_labels[data_set]
                if len(data_types) > 1:
                    lab += f' {data_type.capitalize()}'
                df_set = df_set.sort_values(by=['total_protons'])
                amp, spread = df_set['amp'].iloc[0], df_set['spread'].iloc[0]
                data.append((df_set, lab, data_set, data_type, amp, spread, energy))

    if plot:
        fig, ax = plt.subplots(figsize=(6.66, 5), dpi=144)
        fig.canvas.manager.set_window_title(f'dvar vs Total Protons {energy}GeV {div}°')
        ax.set_title(f'{energy} GeV, {cent_map[cent]} Centrality, {div}° Partitions, 72 Samples per Event')
        ax.set_xlabel('Total Protons in Event')
        if len(data_types) == 1:
            if data_types[0] == 'raw':
                ax.set_ylabel(rf'Single Event $\Delta \sigma^2$')
            elif data_types[0] == 'mix':
                ax.set_ylabel(rf'Mixed Event $\Delta \sigma^2$')
            else:
                ax.set_ylabel(rf'$\Delta \sigma^2$')
        else:
            ax.set_ylabel(rf'$\Delta \sigma^2$')
        ax.axhline(0, ls='-', color='gray')
        if y_ranges:
            ax.set_ylim(y_ranges)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))

    avgs = []
    for i, (df, lab, data_set, data_type, amp, spread, energy) in enumerate(data):
        zo = len(data) - i + 4
        if plot:
            if data_sets_colors is None:
                if len(data_types) > 1:
                    if data_type == 'raw':
                        c = 'blue'
                    elif data_type == 'mix':
                        c = 'green'
                    elif data_type in ['diff', 'sub']:
                        c = 'red'
                else:
                    c = next(color)
            else:
                c = data_sets_colors[data_set]
            if 'sim_' in data_set:
                ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'], label=lab, color=c,
                                alpha=0.4)
            else:
                ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                            marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                if 'sys' in df:
                    ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.4, zorder=zo)

        if avg and len(df) > 1:
            weight_avg = np.average(df['val'], weights=1 / df['err'] ** 2)
            weight_avg_err = np.sqrt(1 / np.average(1 / df['err'] ** 2))
            avgs.append({'name': data_set, 'energy': energy, 'divs': div, 'amp': amp, 'spread': spread,
                         'avg': weight_avg, 'avg_err': weight_avg_err, 'avg_meas': Measure(weight_avg, weight_avg_err)})
            if plot:
                ax.axhline(weight_avg, ls='--', color=c)
                ax.axhspan(weight_avg + weight_avg_err, weight_avg - weight_avg_err, alpha=0.4, color=c)
                if hist:
                    sigs = (df['val'] - weight_avg) / df['err']  # Should probably add weight_avg_err in quad
                    fig_hist, ax_hist = plt.subplots()
                    ax_hist.set_title(f'{lab}')
                    ax_hist.set_xlabel('Standard Deviations from Average')
                    sns.histplot(sigs, stat='density', kde=True)
                    x_norm = np.linspace(min(sigs), max(sigs), 1000)
                    ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                    ax_hist.legend()
                    fig_hist.tight_layout()

    if plot:
        if star_prelim:
            ax.text(10, 0.955, 'STAR Preliminary', fontsize=15)
            eta_line = r'|$\eta$| < 1'
            pt_line = r'0.4 < $p_T$ < 2.0 GeV'
            ax.text(10, 0.962, f'Au+Au\n{eta_line}\n{pt_line}')
            ax.text(31.6, 1.041, 'statistical uncertainty only', fontsize=8)
            ax.legend(loc='upper center', framealpha=1.0).set_zorder(10)
        else:
            ax.legend(loc='upper left')
        fig.tight_layout()

    return pd.DataFrame(avgs)


def stat_vs_protons_divs(df, stat, divs, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False,
                         fit=False, hist=False, data_sets_colors=None, data_sets_labels=None):
    div_data = []
    for div in divs:
        data = []
        for data_type in data_types:
            for data_set in data_sets_plt:
                for energy in energies:
                    df_pre = df
                    if 'data_type' in df_pre:
                        df_pre = df_pre[df_pre['data_type'] == data_type]
                    if 'cent' in df_pre:
                        df_pre = df_pre[df_pre['cent'] == cent]
                    if 'stat' in df_pre:
                        df_pre = df_pre[df_pre['stat'] == stat]
                    df_set = df_pre[
                        (df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                    if len(df_set) == 0:
                        continue
                    if data_sets_labels is None:
                        if energy == 'sim':
                            lab = f'{data_set}_{data_type}'
                        else:
                            lab = f'{data_set}_{data_type}_{energy}GeV'
                    else:
                        lab = data_sets_labels[data_set]
                    df_set = df_set.sort_values(by=['total_protons'])
                    data.append((df_set, lab, data_set, df_set['amp'].iloc[0], df_set['spread'].iloc[0]))
        div_data.append(data)

    if plot:
        fig, ax_divs = plt.subplots(4, 3, sharex=True, sharey=True, figsize=(13.33, 6.16), dpi=144)
        ax_divs = ax_divs.flat
        fig.suptitle(f'{energy}GeV')
        # ax.set_title(f'{stat} {div}° Divisions')
        for ax in ax_divs[-3:]:
            ax.set_xlabel('Total Protons in Event')
        for i, ax in enumerate(ax_divs):
            ax.axhline(1, ls='--', color='black')
            if i == 6:
                ax.set_ylabel(f'Raw / Mix {stat}')
            # if y_ranges:
            #     ax.set_ylim(y_ranges[stat])
            if stat == 'standard deviation':
                ax.set_ylim(0.921, 1.019)
            elif stat == 'skewness':
                ax.set_ylim(0.71, 1.09)
            elif stat == 'non-excess kurtosis':
                ax.set_ylim(0.92, 1.06)
            ax.set_xlim(0, 65)

    fits = []
    for data, ax, div in zip(div_data, ax_divs[:len(divs)], divs):
        if plot:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        for i, (df, lab, data_set, amp, spread) in enumerate(data):
            zo = len(data) - i + 4
            if plot:
                if data_sets_colors is None:
                    c = next(color)
                else:
                    c = data_sets_colors[data_set]
                ax.text(5, .94, f'{div}°')
                if 'sim_' in data_set:
                    ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'],
                                    label=lab, color=c, alpha=0.4)
                else:
                    ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                                marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                    if 'sys' in df:
                        ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='',
                                    elinewidth=3, color=c, alpha=0.4, zorder=zo)

            if fit and len(df) > 1:
                # popt, pcov = cf(line, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
                popt, pcov = cf(line_yint1, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
                popt = np.append(popt, 1)
                pcov = np.array(((pcov[0][0], 0), (0, 0)))
                fits.append({'name': data_set, 'divs': div, 'amp': amp, 'spread': spread, 'slope': popt[0],
                             'slope_err': np.sqrt(np.diag(pcov))[0]})
                if plot:
                    ax.plot(df['total_protons'], line(df['total_protons'], *popt), ls='--', color=c)
                    if hist:
                        sigs = (df['val'] - line(df['total_protons'], *popt)) / df['err']
                        fig_hist, ax_hist = plt.subplots()
                        ax_hist.set_title(f'{lab}')
                        ax_hist.set_xlabel('Standard Deviations from Linear Fit')
                        sns.histplot(sigs, stat='density', kde=True)
                        x_norm = np.linspace(min(sigs), max(sigs), 1000)
                        ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                        ax_hist.legend()
                        fig_hist.tight_layout()

    if plot:
        data, ax, div = div_data[-1], ax_divs[-1], divs[-1]
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        for df, lab, data_set, amp, spread in data:
            if data_sets_colors is None:
                c = next(color)
            else:
                c = data_sets_colors[data_set]
            if 'sim_' in data_set:
                ax.fill_between([], [], [], label=lab, color=c, alpha=0.4)
            else:
                ax.errorbar([], [], [], label=lab, marker='o', ls='', color=c, alpha=0.7)
                if 'sys' in df:
                    ax.errorbar([], [], [], marker='', ls='', elinewidth=3, color=c, alpha=0.4)
        ax_divs[-1].legend()
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.canvas.manager.set_window_title(f'binom_slices_{energies[0]}GeV_{stat}')

    return pd.DataFrame(fits)


def stat_vs_protons_energies(df, stat, divs, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False,
                             fit=False, plot_fit=False, hist=False, data_sets_colors=None, data_sets_labels=None,
                             star_prelim=False):
    energy_data = []
    for energy in energies:
        data = []
        for data_type in data_types:
            for data_set in data_sets_plt:
                for div in divs:
                    df_pre = df
                    if 'data_type' in df_pre:
                        df_pre = df_pre[df_pre['data_type'] == data_type]
                    if 'cent' in df_pre:
                        df_pre = df_pre[df_pre['cent'] == cent]
                    if 'stat' in df_pre:
                        df_pre = df_pre[df_pre['stat'] == stat]
                    df_set = df_pre[
                        (df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                    if len(df_set) == 0:
                        continue
                    if data_sets_labels is None:
                        if energy == 'sim':
                            lab = f'{data_set}_{data_type}'
                        else:
                            lab = f'{data_set}_{data_type}_{energy}GeV'
                    else:
                        lab = data_sets_labels[data_set]
                    df_set = df_set.sort_values(by=['total_protons'])
                    data.append((df_set, lab, data_set, df_set['amp'].iloc[0], df_set['spread'].iloc[0]))
        energy_data.append(data)

    if plot or plot_fit:
        fig, ax_energies = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(13.33, 6.16), dpi=144)
        ax_energies = ax_energies.flat
        fig.suptitle(f'0-5% Centrality, {div}° Partitions, 72 Samples per Event')
        # ax.set_title(f'{stat} {div}° Divisions')
        for ax in ax_energies[-3:]:
            ax.set_xlabel('Total Protons in Event')
        for i, ax in enumerate(ax_energies):
            ax.axhline(1, ls='-', color='gray')
            if i in [0, 3]:
                ax.set_ylabel(f'Raw / Mix {stat.title()}')
            # if y_ranges:
            #     ax.set_ylim(y_ranges[stat])
            if stat == 'standard deviation':
                ax.set_ylim(0.921, 1.019)
            elif stat == 'skewness':
                ax.set_ylim(0.71, 1.09)
            elif stat == 'non-excess kurtosis':
                ax.set_ylim(0.92, 1.06)
            ax.set_xlim(0, 69)

    fits = []
    for data, ax, energy in zip(energy_data, ax_energies[:len(energies)], energies):
        if plot or plot_fit:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        for i, (df, lab, data_set, amp, spread) in enumerate(data):
            zo = len(data) - i + 4
            if plot or plot_fit:
                if data_sets_colors is None:
                    c = next(color)
                else:
                    c = data_sets_colors[data_set]
                ax.text(12, .925, f'{energy} GeV', size='x-large')
            if plot:
                if 'sim_' in data_set:
                    ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'],
                                    label=lab, color=c, alpha=0.4)
                else:
                    ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                                marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                    if 'sys' in df:
                        ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='',
                                    elinewidth=3, color=c, alpha=0.4, zorder=zo)

            if fit and len(df) > 1:
                # popt, pcov = cf(line, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
                popt, pcov = cf(line_yint1, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
                popt = np.append(popt, 1)
                pcov = np.array(((pcov[0][0], 0), (0, 0)))
                fits.append({'name': data_set, 'divs': div, 'amp': amp, 'spread': spread, 'slope': popt[0],
                             'slope_err': np.sqrt(np.diag(pcov))[0]})
                if plot_fit:
                    if plot:
                        ax.plot(df['total_protons'], line(df['total_protons'], *popt), ls='--', color=c)
                    else:
                        ax.plot(df['total_protons'], line(df['total_protons'], *popt), ls='--', color=c, label=lab)
                    if hist:
                        sigs = (df['val'] - line(df['total_protons'], *popt)) / df['err']
                        fig_hist, ax_hist = plt.subplots()
                        ax_hist.set_title(f'{lab}')
                        ax_hist.set_xlabel('Standard Deviations from Linear Fit')
                        sns.histplot(sigs, stat='density', kde=True)
                        x_norm = np.linspace(min(sigs), max(sigs), 1000)
                        ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                        ax_hist.legend()
                        fig_hist.tight_layout()

    if plot or plot_fit:
        # data, ax, energy = energy_data[-1], ax_energies[-1], energies[-1]
        # color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        # for df, lab, data_set, amp, spread in data:
        #     if data_sets_colors is None:
        #         c = next(color)
        #     else:
        #         c = data_sets_colors[data_set]
        #     if 'sim_' in data_set:
        #         ax.fill_between([], [], [], label=lab, color=c, alpha=0.4)
        #     else:
        #         ax.errorbar([], [], [], label=lab, marker='o', ls='', color=c, alpha=0.7)
        #         if 'sys' in df:
        #             ax.errorbar([], [], [], marker='', ls='', elinewidth=3, color=c, alpha=0.4)
        ax_energies[-3].legend(loc='lower right', framealpha=1.0).set_zorder(10)
        if star_prelim:
            ax_energies[4].text(46, 0.93, 'STAR \nPreliminary', fontsize=15)
            eta_line = r'|$\eta$| < 1'
            pt_line = r'0.4 < $p_T$ < 2.0 GeV'
            ax_energies[4].text(46, 0.95, f'Au+Au\n{eta_line}\n{pt_line}')
            ax_energies[2].text(15, 1.012, 'statistical uncertainty only')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.canvas.manager.set_window_title(f'binom_slices_{divs[0]}_{stat}')

    return pd.DataFrame(fits)


def dvar_vs_protons_energies(df, divs, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False,
                             avg=False, plot_avg=False, hist=False, data_sets_colors=None, data_sets_labels=None,
                             star_prelim=False):
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}
    energy_data = []
    for energy in energies:
        data = []
        for data_type in data_types:
            for data_set in data_sets_plt:
                for div in divs:
                    df_pre = df
                    if 'data_type' in df_pre:
                        df_pre = df_pre[df_pre['data_type'] == data_type]
                    if 'cent' in df_pre:
                        df_pre = df_pre[df_pre['cent'] == cent]
                    df_set = df_pre[
                        (df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                    if len(df_set) == 0:
                        continue
                    if data_sets_labels is None:
                        if energy == 'sim':
                            lab = f'{data_set}_{data_type}'
                        else:
                            lab = f'{data_set}_{data_type}_{energy}GeV'
                    else:
                        lab = data_sets_labels[data_set]
                    if len(data_types) > 1:
                        lab += f' {data_type.capitalize()}'
                    df_set = df_set.sort_values(by=['total_protons'])
                    data.append((df_set, lab, data_set, df_set['amp'].iloc[0], df_set['spread'].iloc[0]))
        energy_data.append(data)

    if plot or plot_avg:
        fig, ax_energies = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(13.33, 6.16), dpi=144)
        ax_energies = ax_energies.flat
        fig.suptitle(f'{cent_map[cent]} Centrality, {div}° Partitions, 72 Samples per Event')
        for ax in ax_energies[-3:]:
            ax.set_xlabel('Total Protons in Event')
        for i, ax in enumerate(ax_energies):
            ax.axhline(0, ls='-', color='gray')
            if i in [0, 3]:
                if len(data_types) == 1:
                    if data_types[0] == 'raw':
                        ax.set_ylabel(rf'Single Event $\Delta \sigma^2$')
                    elif data_types[0] == 'mix':
                        ax.set_ylabel(rf'Mixed Event $\Delta \sigma^2$')
                    else:
                        ax.set_ylabel(rf'$\Delta \sigma^2$')
                else:
                    ax.set_ylabel(rf'$\Delta \sigma^2$')
            if y_ranges:
                ax.set_ylim(y_ranges)
            ax.set_xlim(0, 69)

    avgs = []
    for data, ax, energy in zip(energy_data, ax_energies[:len(energies)], energies):
        if plot or plot_avg:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        for i, (df, lab, data_set, amp, spread) in enumerate(data):
            zo = len(data) - i + 4
            if plot or plot_avg:
                if data_sets_colors is None:
                    c = next(color)
                else:
                    c = data_sets_colors[data_set]
                ax.text(0.5, 0.95, f'{energy} GeV', size='x-large',
                        horizontalalignment='center', verticalalignment='top', transform=ax.transAxes)
            if plot:
                if 'sim_' in data_set:
                    ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'],
                                    label=lab, color=c, alpha=0.4)
                else:
                    ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                                marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                    if 'sys' in df:
                        ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='',
                                    elinewidth=3, color=c, alpha=0.4, zorder=zo)

            if avg and len(df) > 1:
                weight_avg = np.average(df['val'], weights=1 / df['err'] ** 2)
                weight_avg_err = np.sqrt(1 / np.average(1 / df['err'] ** 2))
                avgs.append({'name': data_set, 'energy': energy, 'divs': div, 'amp': amp, 'spread': spread,
                             'avg': weight_avg, 'avg_err': weight_avg_err,
                             'avg_meas': Measure(weight_avg, weight_avg_err)})
                if plot_avg:
                    if plot:
                        ax.axhline(weight_avg, ls='--', color=c)
                        ax.axhspan(weight_avg + weight_avg_err, weight_avg - weight_avg_err, alpha=0.4, color=c)
                    else:
                        ax.axhline(weight_avg, ls='--', color=c, label=lab)
                        ax.axhspan(weight_avg + weight_avg_err, weight_avg - weight_avg_err, alpha=0.4, color=c)
                    if hist:
                        sigs = (df['val'] - weight_avg) / df['err']  # Should probably add weight_avg_err in quad
                        fig_hist, ax_hist = plt.subplots()
                        ax_hist.set_title(f'{lab}')
                        ax_hist.set_xlabel('Standard Deviations from Average')
                        sns.histplot(sigs, stat='density', kde=True)
                        x_norm = np.linspace(min(sigs), max(sigs), 1000)
                        ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                        ax_hist.legend()
                        fig_hist.tight_layout()

    if plot or plot_avg:
        ax_energies[-3].legend(loc='lower right', framealpha=1.0).set_zorder(10)
        if star_prelim:
            ax_energies[4].text(46, 0.93, 'STAR \nPreliminary', fontsize=15)
            eta_line = r'|$\eta$| < 1'
            pt_line = r'0.4 < $p_T$ < 2.0 GeV'
            ax_energies[4].text(46, 0.95, f'Au+Au\n{eta_line}\n{pt_line}')
            ax_energies[2].text(15, 1.012, 'statistical uncertainty only')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.canvas.manager.set_window_title(f'binom_slices_{divs[0]}')

    return pd.DataFrame(avgs)


def stat_vs_protons_cents(df, stat, divs, cents, energy, data_types, data_sets_plt, y_ranges=None, plot=False,
                          fit=False, plot_fit=False, hist=False, data_sets_colors=None, data_sets_labels=None):
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}
    cent_data = []
    for cent in cents:
        data = []
        for data_type in data_types:
            for data_set in data_sets_plt:
                for div in divs:
                    df_pre = df
                    if 'data_type' in df_pre:
                        df_pre = df_pre[df_pre['data_type'] == data_type]
                    if 'cent' in df_pre:
                        df_pre = df_pre[df_pre['cent'] == cent]
                    if 'stat' in df_pre:
                        df_pre = df_pre[df_pre['stat'] == stat]
                    df_set = df_pre[
                        (df_pre['name'] == data_set) & (df_pre['divs'] == div) & (df_pre['energy'] == energy)]
                    if len(df_set) == 0:
                        continue
                    if data_sets_labels is None:
                        if energy == 'sim':
                            lab = f'{data_set}_{data_type}'
                        else:
                            lab = f'{data_set}_{data_type}_{energy}GeV'
                    else:
                        lab = data_sets_labels[data_set]
                    df_set = df_set.sort_values(by=['total_protons'])
                    data.append((df_set, div, lab, data_set, df_set['amp'].iloc[0], df_set['spread'].iloc[0]))
        cent_data.append(data)

    if plot or plot_fit:
        fig, ax_cents = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(13.33, 6.16), dpi=144)
        ax_cents = ax_cents.flat
        fig.suptitle(f'{energy}GeV, {div}° Partitions, 72 Samples per Event')
        # ax.set_title(f'{stat} {div}° Divisions')
        for ax in ax_cents[-3:]:
            ax.set_xlabel('Total Protons in Event')
        ax_cents[3].set_ylabel(f'Raw / Mix {stat.title()}')
        for i, ax in enumerate(ax_cents):
            ax.axhline(1, ls='-', color='gray')
            # if y_ranges:
            #     ax.set_ylim(y_ranges[stat])
            if stat == 'standard deviation':
                ax.set_ylim(0.921, 1.019)
            elif stat == 'skewness':
                ax.set_ylim(0.71, 1.09)
            elif stat == 'non-excess kurtosis':
                ax.set_ylim(0.92, 1.06)
            ax.set_xlim(0, 65)

    fits = []
    for data, ax_index, cent in zip(cent_data, range(len(cents)), cents):
        if plot or plot_fit:
            ax = ax_cents[ax_index]
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
        for i, (df, div, lab, data_set, amp, spread) in enumerate(data):
            zo = len(data) - i + 4
            if plot or plot_fit:
                if data_sets_colors is None:
                    c = next(color)
                else:
                    c = data_sets_colors[data_set]
                ax.text(62, .925, f'{cent_map[cent]}', size='x-large', horizontalalignment='right')
            if plot:
                if 'sim_' in data_set:
                    ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'],
                                    label=lab, color=c, alpha=0.4)
                else:
                    ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                                marker='o', ls='', color=c, alpha=0.7, zorder=zo)
                    if 'sys' in df:
                        ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='',
                                    elinewidth=3, color=c, alpha=0.4, zorder=zo)

            if fit and len(df) > 1:
                popt, pcov = cf(line_yint1, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
                fits.append({'name': data_set, 'energy': energy, 'divs': div, 'amp': amp, 'spread': spread,
                             'slope': popt[0], 'slope_err': np.sqrt(pcov[0])[0], 'cent': cent})
                if plot_fit:
                    line_vals = line_yint1(df['total_protons'], *popt)
                    line_high = line_yint1(df['total_protons'], popt[0] + np.sqrt(pcov[0]))
                    line_low = line_yint1(df['total_protons'], popt[0] - np.sqrt(pcov[0]))
                    if plot:
                        ax.fill_between(df['total_protons'], line_high, line_low, color=c, alpha=0.3)
                        ax.plot(df['total_protons'], line_vals, ls='--', color=c)
                    else:
                        ax.fill_between(df['total_protons'], line_high, line_low, color=c, alpha=0.3)
                        ax.plot(df['total_protons'], line_yint1(df['total_protons'], *popt), ls='--', color=c,
                                label=lab)
                    if hist:
                        sigs = (df['val'] - line_yint1(df['total_protons'], *popt)) / df['err']
                        fig_hist, ax_hist = plt.subplots()
                        ax_hist.set_title(f'{lab}')
                        ax_hist.set_xlabel('Standard Deviations from Linear Fit')
                        sns.histplot(sigs, stat='density', kde=True)
                        x_norm = np.linspace(min(sigs), max(sigs), 1000)
                        ax_hist.plot(x_norm, norm(0, 1).pdf(x_norm), color='red', label='Standard Normal')
                        ax_hist.legend()
                        fig_hist.tight_layout()

    if plot or plot_fit:
        for i, (df, div, lab, data_set, amp, spread) in enumerate(data):  # Just to get legend in last axis
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
            if data_sets_colors is None:
                c = next(color)
            else:
                c = data_sets_colors[data_set]
            if 'sim_' in data_set:
                ax_cents[-1].fill_between([], [], [], label=lab, color=c, alpha=0.4)
            else:
                ax_cents[-1].errorbar([], [], [], label=lab, marker='o', ls='', color=c, alpha=0.7)
                if 'sys' in df:
                    ax_cents[-1].errorbar([], [], [], marker='', ls='', elinewidth=3, color=c, alpha=0.4)
        ax_cents[-1].legend()
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.canvas.manager.set_window_title(f'binom_slices_cents_{divs[0]}_{stat}')

    return pd.DataFrame(fits)


def plot_protons_fits_divs(df, data_sets_plt, fit=False, data_sets_colors=None, data_sets_labels=None, exclude_divs=[],
                           verbose=True, plt_energies=True, title=None, alpha=1):
    energies = pd.unique(df['energy'])
    if plt_energies:
        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title(f'Slope vs Width All Energies')
        ax.axhline(0, ls='-', color='black')
        fig_panels, ax_panels = plt.subplots(2, 3, sharex=True, sharey=True, dpi=144, figsize=(13.33, 6.16))
        fig_panels.canvas.manager.set_window_title(f'Slope vs Width Energy Panels')
        ax_panels = dict(zip(energies, ax_panels.flat))
    energy_fig_axs = {energy: plt.subplots() for energy in energies}
    markers = ['o', 's', 'P', 'D', '*', '^', 'p']
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(energies) * len(data_sets_plt))))
    fit_pars = []
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='divs')
        for energy_marker, energy in enumerate(energies):
            df_energy = df_set[df_set['energy'] == energy]
            df_energy.sort_values(by='divs')
            if data_sets_labels is None:
                lab_energy = data_set
            else:
                lab_energy = data_sets_labels[data_set]
            if len(energies) > 1:
                lab = f'{lab_energy}_{energy}GeV'
            else:
                lab = lab_energy
            if data_sets_colors is None:
                color = next(colors)
            else:
                color = data_sets_colors[data_set]
            energy_fig, energy_ax = energy_fig_axs[energy]
            energy_ax.errorbar(df_energy['divs'], df_energy['slope'], yerr=df_energy['slope_err'], ls='none',
                               marker='o', label=lab_energy, color=color, alpha=alpha)
            if plt_energies:
                ax.errorbar(df_energy['divs'], df_energy['slope'], yerr=df_energy['slope_err'], ls='none',
                            marker=markers[energy_marker], label=lab, color=color, alpha=alpha)
                ax_panels[energy].errorbar(df_energy['divs'], df_energy['slope'], yerr=df_energy['slope_err'],
                                           ls='none',
                                           marker='o', label=lab_energy, color=color, alpha=alpha)
            if fit and df_energy.size > 1:
                try:
                    df_energy = df_energy[~df_energy.divs.isin(exclude_divs)]
                    popt, pcov = cf(quad_180, df_energy['divs'], df_energy['slope'], sigma=df_energy['slope_err'],
                                    absolute_sigma=True)
                    perr = np.sqrt(np.diag(pcov))
                    fit_pars.append({'data_set': data_set, 'energy': energy, 'curvature': popt[0], 'color': color,
                                     'curve_baseline': popt[1], 'curve_err': perr[0], 'curve_base_err': perr[1],
                                     'spread': df_energy['spread'].iloc[0], 'amp': df_energy['amp'].iloc[0],
                                     'ap': df_energy['ap'].iloc[0], 'am': df_energy['am'].iloc[0],
                                     'sp': df_energy['sp'].iloc[0], 'sm': df_energy['sm'].iloc[0]})
                    x = np.linspace(0, 360, 100)
                    energy_ax.plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                    if plt_energies:
                        ax.plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                        ax_panels[energy].plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                    # p0 = [*popt, -popt[0], popt[1]]
                    # popt3, pcov3 = cf(double_quad_180_zparam, df_energy['divs'], df_energy['slope'], p0=p0,
                    #                   sigma=df_energy['slope_err'], absolute_sigma=True)
                    # perr3 = np.sqrt(np.diag(pcov))
                    # ax.plot(x, double_quad_180_zparam(x, *popt3), ls='--', color=color, alpha=0.65)
                    # energy_ax.plot(x, double_quad_180_zparam(x, *popt3), ls='--', color=color, alpha=0.65)
                    if popt[0] * popt[1] < 0:
                        popt2, pcov2 = cf(quad_180_zparam, df_energy['divs'], df_energy['slope'],
                                          sigma=df_energy['slope_err'],
                                          absolute_sigma=True)
                        perr2 = np.sqrt(np.diag(pcov2))
                        fit_pars[-1].update({'zero_mag': popt2[0], 'zero_mag_err': perr2[0], 'baseline': popt2[1],
                                             'base_err': perr2[1]})
                        if verbose:
                            print(
                                f'{data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\nz-c fit: {popt2}\n'
                                f'z-c covariance: {pcov2}')
                            print()
                    else:
                        if verbose:
                            print(f'No Zeros! {data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\n')
                            print()
                except RuntimeError as e:
                    print(f'Fitting Error, skipping data_set {data_set}, {e}')

    if title is None:
        title = ''
        if len(data_sets_plt) == 1:
            title += f'{data_set}'
        if len(energies) == 1:
            if title != '':
                title += ' '
            title += f'{energies[0]}GeV'

    if plt_energies:
        if title != '' and title is not None:
            ax.set_title(title)
        ax.set_ylabel('Slope of Raw/Mix SD vs Total Protons per Event')
        ax.set_xlabel('Azimuthal Partition Width')
        ax.legend()

    for energy_i, (energy, (energy_fig, energy_ax)) in enumerate(energy_fig_axs.items()):
        if title is None:
            energy_ax.set_title(f'{energy} GeV')
        else:
            energy_ax.set_title(title)
        energy_ax.set_xlabel('Azimuthal Partition Width (w)')
        energy_ax.set_ylabel('Slope of k2 vs Total Protons per Event')
        energy_ax.axhline(0, color='black', zorder=0)
        energy_ax.legend()
        energy_fig.tight_layout()
        if title is None:
            energy_fig.canvas.manager.set_window_title(f'Slope vs Width {energy}GeV')
        else:
            energy_fig.canvas.manager.set_window_title(f'Slope vs Width {title}')

        if plt_energies:
            ax_panels[energy].axhline(0, color='black', zorder=0)
            ax_panels[energy].text(180, 0.00025, f'{energy} GeV', ha='center', fontsize=14)
            if energy_i >= 3:
                ax_panels[energy].set_xlabel('Azimuthal Partition Width')
            if energy_i in [0, 3]:
                ax_panels[energy].set_ylabel('Slope of Raw/Mix SD vs Protons/Event')
            if energy_i == 1:
                ax_panels[energy].legend(loc='upper center', bbox_to_anchor=(0.5, 0.85), framealpha=1.0)

    if plt_energies:
        fig.tight_layout()
        fig_panels.tight_layout()
        fig_panels.subplots_adjust(wspace=0.0, hspace=0.0)

    return pd.DataFrame(fit_pars)


def plot_dvar_avgs_divs(df, data_sets_plt, fit=False, data_sets_colors=None, data_sets_labels=None, exclude_divs=[],
                        verbose=True, plt_energies=True, title=None, alpha=1):
    energies = pd.unique(df['energy'])
    if plt_energies:
        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title(f'dsigma^2 vs Width All Energies')
        ax.axhline(0, ls='-', color='black')
        fig_panels, ax_panels = plt.subplots(2, 3, sharex=True, sharey=True, dpi=144, figsize=(13.33, 6.16))
        fig_panels.canvas.manager.set_window_title(f'dsigma^2 vs Width Energy Panels')
        ax_panels = dict(zip(energies, ax_panels.flat))
    energy_fig_axs = {energy: plt.subplots() for energy in energies}
    markers = ['o', 's', 'P', 'D', '*', '^', 'p']
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(energies) * len(data_sets_plt))))
    fit_pars = []
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='divs')
        for energy_marker, energy in enumerate(energies):
            df_energy = df_set[df_set['energy'] == energy]
            df_energy.sort_values(by='divs')
            if data_sets_labels is None:
                lab_energy = data_set
            else:
                lab_energy = data_sets_labels[data_set]
            if len(energies) > 1:
                lab = f'{lab_energy}_{energy}GeV'
            else:
                lab = lab_energy
            if data_sets_colors is None:
                color = next(colors)
            else:
                color = data_sets_colors[data_set]
            energy_fig, energy_ax = energy_fig_axs[energy]
            energy_ax.errorbar(df_energy['divs'], df_energy['avg'], yerr=df_energy['avg_err'], ls='none',
                               marker='o', label=lab_energy, color=color, alpha=alpha)
            if plt_energies:
                ax.errorbar(df_energy['divs'], df_energy['avg'], yerr=df_energy['avg_err'], ls='none',
                            marker=markers[energy_marker], label=lab, color=color, alpha=alpha)
                ax_panels[energy].errorbar(df_energy['divs'], df_energy['avg'], yerr=df_energy['avg_err'], ls='none',
                                           marker='o', label=lab_energy, color=color, alpha=alpha)
            if fit and df_energy.size > 1:
                try:
                    df_energy = df_energy[~df_energy.divs.isin(exclude_divs)]
                    popt, pcov = cf(quad_180, df_energy['divs'], df_energy['avg'], sigma=df_energy['avg_err'],
                                    absolute_sigma=True)
                    perr = np.sqrt(np.diag(pcov))
                    fit_pars.append({'data_set': data_set, 'energy': energy, 'curvature': popt[0], 'color': color,
                                     'curve_baseline': popt[1], 'curve_err': perr[0], 'curve_base_err': perr[1],
                                     'spread': df_energy['spread'].iloc[0], 'amp': df_energy['amp'].iloc[0]})
                    x = np.linspace(0, 360, 100)
                    energy_ax.plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                    if plt_energies:
                        ax.plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                        ax_panels[energy].plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                    if popt[0] * popt[1] < 0:
                        popt2, pcov2 = cf(quad_180_zparam, df_energy['divs'], df_energy['avg'],
                                          sigma=df_energy['avg_err'], absolute_sigma=True)
                        perr2 = np.sqrt(np.diag(pcov2))
                        fit_pars[-1].update({'zero_mag': popt2[0], 'zero_mag_err': perr2[0], 'baseline': popt2[1],
                                             'base_err': perr2[1]})
                        if verbose:
                            print(f'{data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\n'
                                  f'z-c fit: {popt2}\nz-c covariance: {pcov2}')
                            print()
                    else:
                        if verbose:
                            print(f'No Zeros! {data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\n')
                            print()
                except RuntimeError as e:
                    print(f'Fitting Error, skipping data_set {data_set}, {e}')

    if title is None:
        title = ''
        if len(data_sets_plt) == 1:
            title += f'{data_set}'
        if len(energies) == 1:
            if title != '':
                title += ' '
            title += f'{energies[0]}GeV'

    if plt_energies:
        if title != '' and title is not None:
            ax.set_title(title)
        ax.set_ylabel(r'$\Delta\sigma^2_{single} - \Delta\sigma^2_{mix}$ vs Total Protons per Event')
        ax.set_xlabel('Azimuthal Partition Width')
        ax.legend()

    for energy_i, (energy, (energy_fig, energy_ax)) in enumerate(energy_fig_axs.items()):
        if title is None:
            energy_ax.set_title(f'{energy} GeV')
        else:
            energy_ax.set_title(title)
        energy_ax.set_xlabel('Azimuthal Partition Width (w)')
        energy_ax.set_ylabel(r'$\Delta\sigma^2_{single} - \Delta\sigma^2_{mix}$ vs Total Protons per Event')
        energy_ax.axhline(0, color='black', zorder=0)
        energy_ax.legend()
        energy_fig.tight_layout()
        if title is None:
            energy_fig.canvas.manager.set_window_title(f'dsigma^2 vs Width {energy}GeV')
        else:
            energy_fig.canvas.manager.set_window_title(f'dsigma^2 vs Width {title}')

        if plt_energies:
            ax_panels[energy].axhline(0, color='black', zorder=0)
            ax_panels[energy].text(0.5, 0.95, f'{energy} GeV', size='x-large', ha='center', va='top',
                                   transform=ax_panels[energy].transAxes)
            if energy_i >= 3:
                ax_panels[energy].set_xlabel('Azimuthal Partition Width')
            if energy_i in [0, 3]:
                ax_panels[energy].set_ylabel(r'$\Delta\sigma^2_{single} - \Delta\sigma^2_{mix}$ vs Protons/Event')
            if energy_i == 1:
                ax_panels[energy].legend(loc='upper center', bbox_to_anchor=(0.5, 0.85), framealpha=1.0)

    if plt_energies:
        fig.tight_layout()
        fig_panels.tight_layout()
        fig_panels.subplots_adjust(wspace=0.0, hspace=0.0, left=0.075, right=0.995, top=0.985, bottom=0.075)

    return pd.DataFrame(fit_pars)


def plot_protons_fits_divs_flow(df, data_sets_plt, data_sets_colors=None):
    fig, ax = plt.subplots(dpi=144)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(data_sets_plt))))
    ax.axhline(0, color='black')
    x_divs = np.linspace(0, 360, 1000)
    fit_pars = []
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='divs')
        if data_sets_colors is None:
            color = next(colors)
        else:
            color = data_sets_colors[data_set]
        v2 = None
        for element in data_set.split('_'):
            if 'v2' in element:
                v2 = float('0.' + element.strip('v2'))
                break
        print(data_set, v2)
        if 'anticlflow_' in data_set:
            # v2 = float('0.' + data_set.split('_')[-3][2:])
            lab = f'anticl + v2={v2:.2f}'
        elif 'flow_' in data_set:
            # v2 = float('0.' + data_set.split('_')[-1][2:])
            lab = f'v2={v2:.2f}'
            ax.plot(x_divs, v2_divs(x_divs / 180 * np.pi, v2), color=color)
        else:
            lab = data_set

        ax.errorbar(df_set['divs'], df_set['slope'], yerr=df_set['slope_err'], ls='none',
                    marker='o', label=lab, color=color)

        if 'anticlflow_' in data_set:
            divs = np.array(df_set['divs'])
            slopes = np.array(df_set['slope'])
            y_sub = slopes - v2_divs(divs / 180 * np.pi, v2)
            ax.errorbar(divs, y_sub, yerr=df_set['slope_err'], ls='none',
                        marker='x', label=f'{lab} - v2(div)', color=color)

    ax.set_title('V2 vs Partition Width')
    ax.set_ylabel('Slope of Raw/Mix SD vs Total Protons per Event')
    ax.set_xlabel('Azimuthal Partition Width')
    ax.legend()

    fig.canvas.manager.set_window_title('V2 vs Partition Width')
    fig.tight_layout()

    return pd.DataFrame(fit_pars)


def plot_dsigma_fits_divs_flow(df, data_sets_plt, data_sets_colors=None):
    fig, ax = plt.subplots(dpi=144)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(data_sets_plt))))
    ax.axhline(0, color='black')
    x_divs = np.linspace(0, 360, 1000)
    fit_pars = []
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='divs')
        if data_sets_colors is None:
            color = next(colors)
        else:
            color = data_sets_colors[data_set]
        v2 = 0
        for element in data_set.split('_'):
            if 'v2' in element:
                v2 = float('0.' + element.strip('v2'))
                break
        print(data_set, v2)
        if 'anticlflow_' in data_set:
            # v2 = float('0.' + data_set.split('_')[-3][2:])
            lab = f'anticl + v2={v2:.2f}'
        elif 'flow_' in data_set and abs(v2) > 0:
            # v2 = float('0.' + data_set.split('_')[-1][2:])
            lab = f'v2={v2:.2f}'
            ax.plot(x_divs, v2_divs(x_divs / 180 * np.pi, v2), color=color)
        else:
            lab = data_set

        ax.errorbar(df_set['divs'], df_set['avg'], yerr=df_set['avg_err'], ls='none',
                    marker='o', label=lab, color=color)

        if 'anticlflow_' in data_set:
            divs = np.array(df_set['divs'])
            avgs = np.array(df_set['avg'])
            y_sub = avgs - v2_divs(divs / 180 * np.pi, v2)
            ax.errorbar(divs, y_sub, yerr=df_set['avg_err'], ls='none',
                        marker='x', label=f'{lab} - v2(div)', color=color)

    ax.set_title('V2 vs Partition Width')
    ax.set_ylabel(r'$\Delta\sigma^2$')
    ax.set_xlabel('Azimuthal Partition Width')
    ax.legend()

    fig.canvas.manager.set_window_title('V2 vs Partition Width')
    fig.tight_layout()

    return pd.DataFrame(fit_pars)


def plot_protons_fits_divs_cents(df, data_sets_plt, plot=False, fit=False, data_sets_colors=None, data_sets_labels=None,
                                 exclude_divs=[]):
    if plot:
        fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(13.33, 6.16), dpi=144)
        fig.suptitle(f'Slope vs Partition Width 72 Samples per Event')
        axs = axs.flatten()
        fig.canvas.manager.set_window_title(f'Slope vs Width All Energies')
        markers = ['o', 's', 'P', 'D', '*', '^', 'p']
    energies = pd.unique(df['energy'])
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(energies) * len(data_sets_plt))))
    fit_pars = []
    for ax_index, cent in enumerate(pd.unique(df['cent'])):
        if plot:
            ax = axs[ax_index]
            ax.axhline(0, ls='-', color='black')
            ax.text(125, 0.002, f'Centrality {cent}', size='x-large')
        df_cent = df[df['cent'] == cent]
        for data_set in data_sets_plt:
            df_set = df_cent[df_cent['name'] == data_set]
            for energy_marker, energy in enumerate(energies):
                df_energy = df_set[df_set['energy'] == energy]
                # print(f'{data_set} {energy}GeV Cent {cent}:\n {df_energy}\n')
                df_energy.sort_values(by='divs')
                if data_sets_colors is None:
                    color = next(colors)
                else:
                    color = data_sets_colors[data_set]
                if plot:
                    if data_sets_labels is None:
                        lab_energy = data_set
                    else:
                        lab_energy = data_sets_labels[data_set]
                    if len(energies) > 1:
                        lab = f'{lab_energy}_{energy}GeV'
                    else:
                        lab = lab_energy
                    ax.errorbar(df_energy['divs'], df_energy['slope'], yerr=df_energy['slope_err'], ls='none',
                                marker=markers[energy_marker], label=lab, color=color)
                if fit and df_energy.size > 1:
                    df_energy = df_energy[~df_energy.divs.isin(exclude_divs)]
                    popt, pcov = cf(quad_180, df_energy['divs'], df_energy['slope'], sigma=df_energy['slope_err'],
                                    absolute_sigma=True)
                    perr = np.sqrt(np.diag(pcov))
                    fit_pars.append({'data_set': data_set, 'energy': energy, 'curvature': popt[0], 'color': color,
                                     'curve_base': popt[1], 'curve_err': perr[0], 'curve_base_err': perr[1],
                                     'cent': cent, 'spread': df_energy['spread'].iloc[0],
                                     'amp': df_energy['amp'].iloc[0]})
                    if plot:
                        x = np.linspace(0, 360, 100)
                        ax.plot(x, quad_180(x, *popt), ls='-', color=color, alpha=0.65)
                    if popt[0] * popt[1] < 0:
                        popt2, pcov2 = cf(quad_180_zparam, df_energy['divs'], df_energy['slope'],
                                          sigma=df_energy['slope_err'],
                                          absolute_sigma=True)
                        perr2 = np.sqrt(np.diag(pcov2))
                        fit_pars[-1].update({'zero_mag': popt2[0], 'zero_mag_err': perr2[0], 'baseline': popt2[1],
                                             'base_err': perr2[1]})
                        print(f'{data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\nz-c fit: {popt2}\n'
                              f'z-c covariance: {pcov2}')
                        print()
                    else:
                        print(f'No Zeros! {data_set}, {energy}GeV\na-c fit: {popt}\na-c covariance: {pcov}\n')
                        print()

    if plot:
        for i in [0, 3, 6]:
            axs[i].set_ylabel('Slope')
        for i in [6, 7, 8]:
            axs[i].set_xlabel('Azimuthal Partition Width')
        ax.legend()
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0, hspace=0.0)
        fig.canvas.manager.set_window_title('Slopes vs Partition Width Centralities')

    return pd.DataFrame(fit_pars)


def plot_vs_div_width_comp(df, title=''):
    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title(f'Closest Sims vs Partition Width {title}')
    ax.set_title(title)
    ax.set_xlabel('Azimuthal Partition Width')
    ax.set_ylabel('Raw SD / Mix SD vs Total Protons Slope')
    ax.axhline(0, color='black')

    fig2, ax2 = plt.subplots()
    fig2.canvas.manager.set_window_title(f'Closest Sims Zero Base {title}')
    ax2.set_title(title)
    ax2.set_xlabel('Baseline')
    ax2.set_ylabel('Zeros')
    ax2.set_ylim(110, 230)
    ax2.set_xlim(-0.013, 0.0005)
    ax2.grid()

    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(pd.unique(df['name'])))))
    for data_set in pd.unique(df['name']):
        color = next(colors)
        df_set = df[df['name'] == data_set]

        ax.errorbar(df_set['divs'], df_set['slope'], df_set['slope_err'], ls='none', marker='o', alpha=0.6, color=color,
                    label=data_set)
        popt, pcov = cf(quad_180_zparam, df_set['divs'], df_set['slope'], sigma=df_set['slope_err'],
                        absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        x_plt = np.linspace(0, 360, 1000)
        ax.plot(x_plt, quad_180_zparam(x_plt, *popt), color=color)
        ax2.errorbar(popt[1], popt[0], perr[0], perr[1], color=color, marker='o', ls='none', label=data_set)

    ax.legend()
    ax2.legend()
    fig.tight_layout()
    fig2.tight_layout()


def plot_slope_div_fits(df_fits, data_sets_colors=None, data_sets_labels=None, ref_line=0):
    fig_base_energy, ax_base_energy = plt.subplots()
    ax_base_energy.set_xlabel('Energy (GeV)')
    ax_base_energy.set_ylabel('Baseline')
    ax_base_energy.axhline(ref_line, color='black')
    fig_base_energy.canvas.manager.set_window_title('Slope Baseline vs Energy')

    fig_zeros_energy, ax_zeros_energy = plt.subplots()
    ax_zeros_energy.set_xlabel('Energy (GeV)')
    ax_zeros_energy.set_ylabel('Zeros')
    ax_zeros_energy.axhline(ref_line, color='black')
    fig_zeros_energy.canvas.manager.set_window_title('Slope Zeros vs Energy')

    fig_basez_zeros, ax_basez_zeros = plt.subplots()
    ax_basez_zeros.set_xlabel('Zeros')
    ax_basez_zeros.set_ylabel('Baseline')
    ax_basez_zeros.axvline(0, color='black')
    ax_basez_zeros.axhline(0, color='black')
    fig_basez_zeros.canvas.manager.set_window_title('Slope Baseline vs Zeros')

    colors = ['black', 'red', 'blue', 'green', 'purple', 'orange']
    markers = ['o', 's', '^', 'P', 'p', '2']

    data_sets = pd.unique(df_fits['data_set'])
    for data_set_index, data_set in enumerate(data_sets):
        if data_sets_labels is None:
            lab = data_set
        else:
            lab = data_sets_labels[data_set]
        if data_sets_colors is None:
            color = colors[data_set_index % len(colors)]
        else:
            color = data_sets_colors[data_set]

        df_data_set = df_fits[df_fits['data_set'] == data_set]

        ax_base_energy.errorbar(df_data_set['energy'], df_data_set['baseline'], yerr=df_data_set['base_err'], ls='none',
                                marker='o', label=lab, color=color, alpha=0.8)

        ax_zeros_energy.errorbar(df_data_set['energy'], df_data_set['zero_mag'], yerr=df_data_set['zero_mag_err'],
                                 ls='none', marker='o', label=lab, color=color, alpha=0.8)

        energies = pd.unique(df_data_set['energy'])
        for energy_index, energy in enumerate(energies):
            df_energy = df_data_set[df_data_set['energy'] == energy]
            ax_basez_zeros.errorbar(df_energy['zero_mag'], df_energy['baseline'], xerr=df_energy['zero_mag_err'],
                                    yerr=df_energy['base_err'], ls='none', color=colors[energy_index % len(colors)],
                                    marker=markers[data_set_index % len(markers)], label=f'{lab} {energy}GeV',
                                    alpha=0.8)

    ax_base_energy.legend()
    ax_zeros_energy.legend()
    ax_basez_zeros.legend(loc='center left', bbox_to_anchor=(0.05, 0.5))

    ax_base_energy.grid()
    ax_zeros_energy.grid()
    ax_basez_zeros.grid()

    fig_base_energy.tight_layout()
    fig_zeros_energy.tight_layout()
    fig_basez_zeros.tight_layout()


def plot_slope_div_fits_simpars(df_fits):
    fig_base_amp, ax_base_amp = plt.subplots()
    ax_base_amp.set_xlabel('Amplitude')
    ax_base_amp.set_ylabel('Baseline')
    ax_base_amp.axhline(0, color='black')
    ax_base_amp.axvline(0, color='black')
    fig_base_amp.canvas.manager.set_window_title('Slope Baseline vs Amplitude')

    fig_amp_base, ax_amp_base = plt.subplots()
    ax_amp_base.set_xlabel('Baseline')
    ax_amp_base.set_ylabel('Amplitude')
    ax_amp_base.axvline(0, color='black')
    ax_amp_base.axhline(0, color='black')
    fig_amp_base.canvas.manager.set_window_title('Slope Amplitude vs Baseline')

    fig_base_spread, ax_base_spread = plt.subplots()
    ax_base_spread.set_xlabel('Spread')
    ax_base_spread.set_ylabel('Baseline')
    ax_base_spread.axhline(0, color='black')
    fig_base_spread.canvas.manager.set_window_title('Slope Baseline vs Spread')

    fig_zeros_amp, ax_zeros_amp = plt.subplots()
    ax_zeros_amp.set_xlabel('Amplitude')
    ax_zeros_amp.set_ylabel('Zeros')
    ax_zeros_amp.axhline(0, color='black')
    fig_zeros_amp.canvas.manager.set_window_title('Slope Zeros vs Amplitude')

    fig_zeros_spread, ax_zeros_spread = plt.subplots()
    ax_zeros_spread.set_xlabel('Spread')
    ax_zeros_spread.set_ylabel('Zeros')
    ax_zeros_spread.axhline(0, color='black')
    fig_zeros_spread.canvas.manager.set_window_title('Slope Zeros vs Spread')

    fig_base_zeros_gamp, ax_base_zeros_gamp = plt.subplots()
    ax_base_zeros_gamp.set_ylabel('Zeros')
    ax_base_zeros_gamp.set_xlabel('Baseline')
    ax_base_zeros_gamp.axvline(0, color='black')
    ax_base_zeros_gamp.axhline(0, color='black')
    fig_base_zeros_gamp.canvas.manager.set_window_title('Slope Baseline vs Zeros Amp Sets')

    fig_base_zeros_gspread, ax_base_zeros_gspread = plt.subplots()
    ax_base_zeros_gspread.set_ylabel('Zeros')
    ax_base_zeros_gspread.set_xlabel('Baseline')
    ax_base_zeros_gspread.axvline(0, color='black')
    ax_base_zeros_gspread.axhline(0, color='black')
    fig_base_zeros_gspread.canvas.manager.set_window_title('Slope Baseline vs Zeros Spread Sets')

    amps = pd.unique(df_fits['amp'])
    spreads = pd.unique(df_fits['spread'])

    print(f'Amps: {amps}\nSpreads: {spreads}')

    for amp in amps:
        df_amp = df_fits[df_fits['amp'] == amp]
        ax_base_spread.errorbar(df_amp['spread'], df_amp['baseline'], yerr=df_amp['base_err'], ls='none', marker='o',
                                label=f'A={amp}')
        ax_zeros_spread.errorbar(df_amp['spread'], df_amp['zero_mag'], yerr=df_amp['zero_mag_err'], ls='none',
                                 marker='o', label=f'A={amp}')
        ax_base_zeros_gamp.errorbar(df_amp['baseline'], df_amp['zero_mag'], xerr=df_amp['base_err'],
                                    yerr=df_amp['zero_mag_err'], ls='none', marker='o', label=f'A={amp}')

    sigma_fits = []
    cl_type_data = {'_clmul_': {'name': 'Attractive', 'line': '--'}, '_aclmul_': {'name': 'Repulsive', 'line': ':'}}
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads))))
    for spread in spreads:
        color = next(colors)
        df_spread = df_fits[df_fits['spread'] == spread]
        for cl_type, cl_data in cl_type_data.items():
            df_set = df_spread[df_spread['data_set'].str.contains(cl_type)]
            if df_set.size == 0 or any(np.isnan(df_set['baseline'])) or any(np.isinf(df_set['baseline'])) or \
                    any(np.isnan(df_set['base_err'])) or any(np.isinf(df_set['base_err'])):
                continue
            if cl_type == '_aclmul_':
                lab = f'σ={round(spread, 2)}'
            else:
                lab = None
            ax_base_amp.errorbar(df_set['amp'], df_set['baseline'], yerr=df_set['base_err'], ls='none',
                                 marker='o', label=lab, color=color)
            popt, pcov = cf(line_yint0, df_set['amp'], df_set['baseline'], sigma=df_set['base_err'],
                            absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            x_vals = np.array([0, max(df_set['amp'])])
            ax_base_amp.plot(x_vals, line_yint0(x_vals, *popt), ls=cl_data['line'], color=color)
            ax_base_amp.fill_between(x_vals, line_yint0(x_vals, popt[0] - perr[0]),
                                     line_yint0(x_vals, popt[0] + perr[0]), color=color, alpha=0.3)

            ax_amp_base.errorbar(df_set['baseline'], df_set['amp'], xerr=df_set['base_err'], ls='none',
                                 marker='o', label=lab, color=color)
            inv_slope = 1 / Measure(popt[0], perr[0])
            x_vals = np.array([0, max(df_set['baseline']) * 1.1]) if cl_type == '_clmul_' \
                else np.array([min(df_set['baseline']) * 1.1, 0])
            ax_amp_base.plot(x_vals, line_yint0(x_vals, inv_slope.val), ls=cl_data['line'], color=color)
            ax_amp_base.fill_between(x_vals, line_yint0(x_vals, inv_slope.val - inv_slope.err),
                                     line_yint0(x_vals, inv_slope.val + inv_slope.err), color=color, alpha=0.3)

            avg_zero = np.mean(
                [Measure(val, err) for val, err in zip(df_set['zero_mag'], df_set['zero_mag_err'])])
            ax_zeros_amp.errorbar(df_set['amp'], df_set['zero_mag'], yerr=df_set['zero_mag_err'], ls='none',
                                  marker='o', label=lab, color=color)
            ax_zeros_amp.axhspan(avg_zero.val - avg_zero.err, avg_zero.val + avg_zero.err, alpha=0.3, color=color)
            ax_zeros_amp.axhline(avg_zero.val, ls=cl_data['line'], color=color)
            ax_base_zeros_gspread.errorbar(df_set['baseline'], df_set['zero_mag'], xerr=df_set['base_err'],
                                           yerr=df_set['zero_mag_err'], ls='none', marker='o', label=lab,
                                           color=color)
            sigma_fits.append({'clust_type': cl_data['name'], 'spread': spread, 'zero_mean_val': avg_zero.val,
                               'zero_mean_err': avg_zero.err, 'base_slope_val': popt[0], 'base_slope_err': perr[0],
                               'amp_slope_val': inv_slope.val, 'amp_slope_err': inv_slope.err})

    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads))))  # Go through again to plot horizontal lines
    for spread in spreads:
        color = next(colors)
        df_spread = df_fits[df_fits['spread'] == spread]
        for cl_type, cl_data in cl_type_data.items():
            df_set = df_spread[df_spread['data_set'].str.contains(cl_type)]
            if df_set.size == 0:
                continue
            avg_zero = sum([Measure(val, err) for val, err in zip(df_set['zero_mag'], df_set['zero_mag_err'])])
            avg_zero /= df_set['zero_mag'].size
            ax_base_zeros_lims = ax_base_zeros_gspread.get_xlim()
            ax_base_zeros_mid = (0 - ax_base_zeros_lims[0]) / np.diff(ax_base_zeros_lims)[0]
            if cl_type == '_aclmul_':
                ax_base_zeros_gspread.axhspan(avg_zero.val - avg_zero.err, avg_zero.val + avg_zero.err, alpha=0.3,
                                              xmax=ax_base_zeros_mid, color=color)
                ax_base_zeros_gspread.axhline(avg_zero.val, ls=cl_data['line'], color=color, xmax=ax_base_zeros_mid)
            else:
                ax_base_zeros_gspread.axhspan(avg_zero.val - avg_zero.err, avg_zero.val + avg_zero.err, alpha=0.3,
                                              xmin=ax_base_zeros_mid, color=color)
                ax_base_zeros_gspread.axhline(avg_zero.val, ls=cl_data['line'], color=color, xmin=ax_base_zeros_mid)

    ax_base_amp.legend()
    ax_amp_base.legend()
    ax_base_spread.legend()
    ax_zeros_amp.legend()
    ax_zeros_spread.legend()
    ax_base_zeros_gamp.legend()
    ax_base_zeros_gspread.legend()

    fig_base_amp.tight_layout()
    fig_amp_base.tight_layout()
    fig_base_spread.tight_layout()
    fig_zeros_amp.tight_layout()
    fig_zeros_spread.tight_layout()
    fig_base_zeros_gamp.tight_layout()
    fig_base_zeros_gspread.tight_layout()

    return pd.DataFrame(sigma_fits)


def plot_sigma_fits_interp(sigma_fits):
    fig_spread_zero, ax_spread_zero = plt.subplots()
    ax_spread_zero.set_xlabel('Mean Zero')
    ax_spread_zero.set_ylabel('Spread')
    ax_spread_zero.axhline(0, color='black')
    fig_spread_zero.canvas.manager.set_window_title('Spread vs Mean Zero')

    fig_spread_base_slope, ax_spread_base_slope = plt.subplots()
    ax_spread_base_slope.set_xlabel('Spread')
    ax_spread_base_slope.set_ylabel('Base Slope')
    ax_spread_base_slope.axhline(0, color='black')
    fig_spread_base_slope.canvas.manager.set_window_title('Base Slope vs Spread')

    fig_spread_amp_slope, ax_spread_amp_slope = plt.subplots()
    ax_spread_amp_slope.set_xlabel('Spread')
    ax_spread_amp_slope.set_ylabel('Amp Slope')
    ax_spread_amp_slope.axhline(0, color='black')
    fig_spread_amp_slope.canvas.manager.set_window_title('Amp Slope vs Spread')

    cl_types = pd.unique(sigma_fits['clust_type'])
    interpolations = {cl_type: {} for cl_type in cl_types}
    for cl_type in cl_types:
        print(cl_type)
        df_cl_type = sigma_fits[sigma_fits['clust_type'] == cl_type]
        f_spread_zero = interp1d(df_cl_type['zero_mean_val'], df_cl_type['spread'], kind='cubic')
        x_spread_zero = np.linspace(min(df_cl_type['zero_mean_val']), max(df_cl_type['zero_mean_val']), 1000)
        ax_spread_zero.errorbar(df_cl_type['zero_mean_val'], df_cl_type['spread'], xerr=df_cl_type['zero_mean_err'],
                                ls='none', label=cl_type, marker='o')
        ax_spread_zero.plot(x_spread_zero, f_spread_zero(x_spread_zero))

        f_spread_base = interp1d(df_cl_type['spread'], df_cl_type['base_slope_val'], kind='cubic')
        x_spread_base = np.linspace(min(df_cl_type['spread']), max(df_cl_type['spread']), 1000)
        ax_spread_base_slope.errorbar(df_cl_type['spread'], df_cl_type['base_slope_val'], ls='none', label=cl_type,
                                      marker='o', yerr=df_cl_type['base_slope_err'])
        ax_spread_base_slope.plot(x_spread_base, f_spread_base(x_spread_base))

        f_spread_amp = interp1d(df_cl_type['spread'], df_cl_type['amp_slope_val'], kind='cubic')
        ax_spread_amp_slope.errorbar(df_cl_type['spread'], df_cl_type['amp_slope_val'], ls='none', label=cl_type,
                                     marker='o', yerr=df_cl_type['amp_slope_err'])
        ax_spread_amp_slope.plot(x_spread_base, f_spread_amp(x_spread_base))

        interpolations[cl_type] = {'spread_zero': f_spread_zero, 'amp_slope_spread': f_spread_amp}

    ax_spread_zero.legend()
    ax_spread_base_slope.legend()
    ax_spread_amp_slope.legend()

    fig_spread_zero.tight_layout()
    fig_spread_base_slope.tight_layout()
    fig_spread_amp_slope.tight_layout()

    return interpolations


def plot_protons_fits_sim(df):
    fig_amp, ax_amp = plt.subplots()
    ax_amp.set_ylabel('slope')
    ax_amp.set_xlabel('amp')

    fig_spread, ax_spread = plt.subplots()
    ax_spread.set_ylabel('slope')
    ax_spread.set_xlabel('spread')

    spreads = np.unique(df['spread'])
    for spread in spreads:
        df_spread = df[df['spread'] == spread]
        ax_amp.errorbar(df_spread['amp'], df_spread['slope'], yerr=df_spread['slope_err'], ls='none', marker='o',
                        label=spread)

    amps = np.unique(df['amp'])
    for amp in amps:
        df_amp = df[df['amp'] == amp]
        ax_spread.errorbar(df_amp['spread'], df_amp['slope'], yerr=df_amp['slope_err'], ls='none', marker='o',
                           label=amp)

    ax_amp.legend()
    ax_spread.legend()

    fig_sa, ax_sa = plt.subplots()
    ax_sa.set_xlabel('amp')
    ax_sa.set_ylabel('spread')
    cmap = plt.cm.rainbow
    norm = matplotlib.colors.Normalize(vmin=min(df['slope']), vmax=0)
    for data_set in np.unique(df['name']):
        df_set = df[df['name'] == data_set]
        ax_sa.scatter(df_set['amp'], df_set['spread'], marker='o', color=cmap(norm(df_set['slope'])))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig_sa.colorbar(sm)

    fig_3d = plt.figure()
    ax_3d = plt.axes(projection='3d')
    df_sim = df[df['name'].str.contains('sim_')]
    ax_3d.plot_trisurf(df_sim['amp'], df_sim['spread'], df_sim['slope'], cmap='viridis', edgecolor='none')

    f = interp2d(df_sim['amp'], df_sim['spread'], df_sim['slope'], kind='cubic')
    x = np.linspace(df_sim['amp'].min(), df_sim['amp'].max(), 100)
    y = np.linspace(df_sim['spread'].min(), df_sim['spread'].max(), 100)
    xx, yy = np.meshgrid(x, y)
    fig_3d_interp = plt.figure()
    ax_3d_interp = plt.axes(projection='3d')
    z = f(x, y)
    print(z)
    ax_3d_interp.plot_surface(xx, yy, z, cmap='viridis', edgecolor='none')
    ax_3d_interp.set_xlabel('amp')
    ax_3d_interp.set_ylabel('spread')
    ax_3d_interp.set_zlabel('slope')
    ax_3d_interp.axhline(0.002)


def plot_protons_fits_vs_energy(df, data_sets_plt, data_sets_colors=None, data_sets_labels=None, title=None):
    fig_slope, ax_slope = plt.subplots(figsize=(6.66, 5), dpi=144)
    ax_slope.axhline(0, color='gray')
    fig_slope.canvas.manager.set_window_title(f'Slopes vs Energy')
    fig_int, ax_int = plt.subplots()
    ax_int.axhline(1, color='gray')
    fig_int.canvas.manager.set_window_title(f'Intercepts vs Energy')
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='energy')
        if data_sets_labels is None:
            lab = data_set
        else:
            lab = data_sets_labels[data_set]
        if data_sets_colors is None:
            ax_slope.errorbar(df_set['energy'], df_set['slope'], yerr=df_set['slope_err'], ls='none', marker='o',
                              label=lab)
            ax_int.errorbar(df_set['energy'], df_set['int'], yerr=df_set['int_err'], ls='none', marker='o', label=lab)
        else:
            ax_slope.errorbar(df_set['energy'], df_set['slope'], yerr=df_set['slope_err'], ls='none', marker='o',
                              color=data_sets_colors[data_set], label=lab)
            ax_int.errorbar(df_set['energy'], df_set['int'], yerr=df_set['int_err'], ls='none', marker='o',
                            color=data_sets_colors[data_set], label=lab)
    ax_slope.set_ylabel('Slope of Raw/Mix SD vs Total Protons per Event')
    ax_slope.set_xlabel('Energy (GeV)')
    ax_slope.grid()
    ax_int.set_ylabel('Intercept')
    ax_int.set_xlabel('Energy (GeV)')
    ax_int.grid()
    if title:
        ax_slope.set_title(title)
        ax_int.set_title(title)
    legend_slope = ax_slope.legend()
    # legend_slope.get_frame().set_alpha(0)
    legend_int = ax_int.legend()
    # legend_int.get_frame().set_alpha(0)
    fig_slope.tight_layout()
    fig_int.tight_layout()


def plot_protons_fits_vs_cent(df, data_sets_plt, data_sets_colors=None, data_sets_labels=None, title=None,
                              fit=False, cent_ref=None, ref_type=None, data_sets_energies_cmaps=None, ls='-'):
    fig_slope, ax_slope = plt.subplots(figsize=(6.66, 5), dpi=144)
    ax_slope.axhline(0, color='black')
    fig_slope.canvas.manager.set_window_title(f'Slopes vs Centrality')
    energies = pd.unique(df['energy'])
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        if data_sets_energies_cmaps is not None:
            colors = iter([plt.cm.get_cmap(data_sets_energies_cmaps[data_set])(i)
                           for i in np.linspace(1, 0.4, len(energies))])
        elif data_sets_colors is not None:
            color, colors = data_sets_colors[data_set], None
        else:
            colors, color = None
        for energy in pd.unique(df_set['energy']):
            if colors is not None:
                color = next(colors)
            df_energy = df_set[df_set['energy'] == energy]
            df_energy.sort_values(by='cent')
            if data_sets_labels is None:
                lab = data_set
            else:
                lab = data_sets_labels[data_set]
            if len(energies) > 1:
                lab += f' {energy}GeV'

            if cent_ref is None:
                x = df_energy['cent']
                x_err = None
            else:
                cent_energy = cent_ref[(cent_ref['data_set'] == data_set) & (cent_ref['energy'] == energy)]
                x = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_val'].iloc[0] for cent in
                     df_energy['cent']]
                x_err = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_sd'].iloc[0] for cent in
                         df_energy['cent']]
            ls = 'none' if fit else ls
            if colors is None and color is None:
                ax_slope.errorbar(x, df_energy['slope'], xerr=x_err, yerr=df_energy['slope_err'], marker='o', ls=ls,
                                  label=lab, alpha=0.6)
            else:
                ax_slope.errorbar(x, df_energy['slope'], xerr=x_err, yerr=df_energy['slope_err'], marker='o', ls=ls,
                                  color=color, label=lab, alpha=0.6)
            if fit:
                p0 = [-0.02, 0.0001]
                x_fit = np.linspace(min(x), max(x), 1000)
                odr_model = odr.Model(inv_sqrtx_odr)
                fit_index_highs = range(4, len(x))
                colors_fit = iter(plt.cm.rainbow(np.linspace(0, 1, len(fit_index_highs))))
                for i in fit_index_highs:
                    color_fit = next(colors_fit)
                    sx = None if x_err is None else x_err[:i]
                    odr_data = odr.RealData(x[:i], df_energy['slope'][:i], sx=sx, sy=df_energy['slope_err'][:i])
                    inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                    odr_out = inv_sqrt_odr.run()
                    ax_slope.plot(x_fit, inv_sqrtx(x_fit, *odr_out.beta), alpha=0.6, color=color_fit)
                    ax_slope.axhline(odr_out.beta[1], color=color_fit, ls='--')
                    ax_slope.axhspan(odr_out.beta[1] - odr_out.sd_beta[1], odr_out.beta[1] + odr_out.sd_beta[1],
                                     color=color_fit, alpha=0.4)
                    print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')

                # popt, pcov = cf(inv_sqrtx, x, df_energy['slope'], sigma=df_energy['slope_err'], absolute_sigma=True)
                # perr = np.sqrt(np.diag(pcov))
                # x_fit = np.linspace(min(x), max(x), 1000)
                # ax_slope.plot(x_fit, inv_sqrtx(x_fit, *popt), alpha=0.6, color=color)
                # ax_slope.axhline(popt[1], color=color, ls='--')
                # ax_slope.axhspan(popt[1] - perr[1], popt[1] + perr[1], color=color, alpha=0.4)
                # print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(popt, perr)]}')
    ax_slope.set_ylabel('Slope of Raw/Mix SD vs Total Protons per Event')
    if ref_type is None:
        ax_slope.set_xlabel('Centrality')
    else:
        ax_slope.set_xlabel('Reference Multiplicity')
    ax_slope.grid()
    if title:
        ax_slope.set_title(title)
    legend_slope = ax_slope.legend()
    # legend_slope.get_frame().set_alpha(0)
    fig_slope.tight_layout()


def plot_div_fits_vs_cent(df, data_sets_plt, data_sets_colors=None, data_sets_labels=None, title=None,
                          fit=False, cent_ref=None, ref_type=None, data_sets_energies_cmaps=None):
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}
    fig_base, ax_base = plt.subplots(figsize=(6.66, 5), dpi=144)
    ax_base.axhline(0, color='black')
    fig_base.canvas.manager.set_window_title(f'Baselines vs Centrality')
    fig_zeros, ax_zeros = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig_zeros.canvas.manager.set_window_title(f'Zeros vs Centrality')
    energies = pd.unique(df['energy'])
    bases, refs, cs = [], [], [0]  # Just to set x and y limits for plot
    consts = {}
    for data_set in data_sets_plt:
        df_set = df[df['data_set'] == data_set]
        if data_sets_energies_cmaps is not None:
            colors = iter([plt.cm.get_cmap(data_sets_energies_cmaps[data_set])(i)
                           for i in np.linspace(1, 0.4, len(energies))])
        elif data_sets_colors is not None:
            color, colors = data_sets_colors[data_set], None
        else:
            colors, color = None
        print(pd.unique(df_set['energy']))
        consts.update({data_set: {'energy': [], 'const_val': [], 'const_err': []}})
        for energy in pd.unique(df_set['energy']):
            if colors is not None:
                color = next(colors)
            df_energy = df_set[df_set['energy'] == energy]
            df_energy.sort_values(by='cent')
            if data_sets_labels is None:
                lab = data_set
            else:
                lab = data_sets_labels[data_set]
            if len(energies) > 1:
                lab += f' {energy}GeV'

            if cent_ref is None:
                x = [cent_map[cent_i] for cent_i in df_energy['cent']]
                x_err = None
            else:
                cent_energy = cent_ref[(cent_ref['data_set'] == data_set) & (cent_ref['energy'] == energy)]
                x = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_val'].iloc[0] for cent in
                     df_energy['cent']]
                x_err = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_sd'].iloc[0] for cent in
                         df_energy['cent']]
            ls = 'none' if fit else '-'
            if colors is None and color is None:
                ax_base.errorbar(x, df_energy['baseline'], xerr=x_err, yerr=df_energy['base_err'], marker='o',
                                 ls=ls, label=lab, alpha=0.6)
                ax_zeros.errobar(x, df_energy['zero_mag'], xerr=x_err, yerr=df_energy['zero_mag_err'], marker='o',
                                 ls='none', label=lab, alpha=0.6)
            else:
                ax_base.errorbar(x, df_energy['baseline'], xerr=x_err, yerr=df_energy['base_err'], marker='o',
                                 ls=ls, color=color, label=lab, alpha=0.6)
                ax_zeros.errorbar(x, df_energy['zero_mag'], xerr=x_err, yerr=df_energy['zero_mag_err'], marker='o',
                                  ls='none', color=color, label=lab, alpha=0.6)
            bases.extend(list(df_energy['baseline']))
            refs.extend(list(x))
            if fit:
                p0 = [-0.02, 0.0001]
                x_fit = np.linspace(1, 800, 2000)
                odr_model = odr.Model(inv_sqrtx_odr)
                odr_data = odr.RealData(x, df_energy['baseline'], sx=x_err, sy=df_energy['base_err'])
                inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                odr_out = inv_sqrt_odr.run()
                ax_base.axhline(odr_out.beta[1], color=color, ls='-')
                ax_base.axhspan(odr_out.beta[1] - odr_out.sd_beta[1], odr_out.beta[1] + odr_out.sd_beta[1],
                                color=color, alpha=0.4)
                ax_base.plot(x_fit, inv_sqrtx(x_fit, *odr_out.beta), alpha=0.6, color=color)
                print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')
                cs.append(odr_out.beta[1])
                consts[data_set]['energy'].append(energy)
                consts[data_set]['const_val'].append(odr_out.beta[1])
                consts[data_set]['const_err'].append(odr_out.sd_beta[1])

                p0 = [-0.02, 0.0001, 1]
                x_fit = np.linspace(1, 800, 2000)
                odr_model = odr.Model(inv_invxpow_odr)
                odr_data = odr.RealData(x, df_energy['baseline'], sx=x_err, sy=df_energy['base_err'])
                inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                odr_out = inv_sqrt_odr.run()
                ax_base.axhline(odr_out.beta[1], color=color, ls=':')
                ax_base.axhspan(odr_out.beta[1] - odr_out.sd_beta[1], odr_out.beta[1] + odr_out.sd_beta[1],
                                color=color, alpha=0.4)
                ax_base.plot(x_fit, inv_invxpow_odr(odr_out.beta, x_fit), alpha=0.6, color=color, ls=':')
                print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')
                # consts[data_set]['energy'].append(energy)
                # consts[data_set]['const_val'].append(odr_out.beta[1])
                # consts[data_set]['const_err'].append(odr_out.sd_beta[1])

                # p0 = [-0.02, 1]
                # x_fit = np.linspace(1, 800, 2000)
                # odr_model = odr.Model(inv_invxpow_noc_odr)
                # odr_data = odr.RealData(x, df_energy['baseline'], sx=x_err, sy=df_energy['base_err'])
                # inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                # odr_out = inv_sqrt_odr.run()
                # ax_base.plot(x_fit, inv_invxpow_noc_odr(odr_out.beta, x_fit), alpha=0.6, color=color, ls=':')
                # print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')

                p0 = [-0.02, 0.0001]
                x_fit = np.linspace(1, 800, 2000)
                odr_model = odr.Model(inv_invx_odr)
                odr_data = odr.RealData(x, df_energy['baseline'], sx=x_err, sy=df_energy['base_err'])
                inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                odr_out = inv_sqrt_odr.run()
                ax_base.axhline(odr_out.beta[1], color=color, ls='--')
                ax_base.axhspan(odr_out.beta[1] - odr_out.sd_beta[1], odr_out.beta[1] + odr_out.sd_beta[1],
                                color=color, alpha=0.4)
                ax_base.plot(x_fit, inv_invx_odr(odr_out.beta, x_fit), alpha=0.6, color=color, ls='--')
                print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')
                # consts[data_set]['energy'].append(energy)
                # consts[data_set]['const_val'].append(odr_out.beta[1])
                # consts[data_set]['const_err'].append(odr_out.sd_beta[1])
                # cs.append(odr_out.beta[1])
                # ax_base.set_ylim(ax_base_ylim)
                # ax_base.set_xlim(ax_base_xlim)
                # fit_index_highs = range(4, len(x))
                # colors_fit = iter(plt.cm.rainbow(np.linspace(0, 1, len(fit_index_highs))))
                # for i in fit_index_highs:
                #     color_fit = next(colors_fit)
                #     odr_data = odr.RealData(x[:i], df_energy['baseline'][:i], sx=x_err[:i],
                #                             sy=df_energy['base_err'][:i])
                #     inv_sqrt_odr = odr.ODR(odr_data, odr_model, beta0=p0, maxit=500)
                #     odr_out = inv_sqrt_odr.run()
                #     ax_base.plot(x_fit, inv_sqrtx(x_fit, *odr_out.beta), alpha=0.6, color=color_fit)
                #     ax_base.axhline(odr_out.beta[1], color=color_fit, ls='--')
                #     ax_base.axhspan(odr_out.beta[1] - odr_out.sd_beta[1], odr_out.beta[1] + odr_out.sd_beta[1],
                #                     color=color_fit, alpha=0.4)
                #     print(f'{lab} Fit: {[Measure(var, err) for var, err in zip(odr_out.beta, odr_out.sd_beta)]}')
    ax_base.set_ylim(min(bases) * 1.1, max(cs) * 1.1)
    try:
        ax_base.set_xlim(0, max(refs) * 1.1)
    except TypeError:
        pass

    ax_base.set_ylabel('Baseline of Slope vs Partition Width Fit')
    ax_zeros.set_ylabel('Zeros of Slope vs Partition Width Fit')
    if ref_type is None:
        ax_base.set_xlabel('Centrality')
        ax_zeros.set_xlabel('Centrality')
    else:
        ax_base.set_xlabel('Reference Multiplicity')
        ax_zeros.set_xlabel('Reference Multiplicity')
    ax_base.grid()
    ax_zeros.grid()
    if title:
        ax_base.set_title(title)
        ax_zeros.set_title(title)
    legend_base = ax_base.legend()
    legend_zeros = ax_zeros.legend()
    # legend_slope.get_frame().set_alpha(0)
    fig_base.tight_layout()
    fig_zeros.tight_layout()

    fig_consts_vs_energy, ax_consts_vs_energy = plt.subplots()
    ax_consts_vs_energy.set_ylabel('Constant Fit Parameter')
    ax_consts_vs_energy.set_xlabel('Energy (GeV)')
    for data_set, data in consts.items():
        ax_consts_vs_energy.errorbar(data['energy'], data['const_val'], yerr=data['const_err'], ls='none',
                                     label=data_set, marker='o')
    ax_consts_vs_energy.legend()
    ax_consts_vs_energy.grid()
    ax_consts_vs_energy.set_title('Fit Constants vs Energy')
    fig_consts_vs_energy.canvas.manager.set_window_title('Fit Constants vs Energy')
    fig_consts_vs_energy.tight_layout()


def plot_protons_fits_vs_amp(df, data_sets_plt, data_sets_colors=None, data_sets_labels=None):
    fig_slope, ax_slope = plt.subplots()
    ax_slope.axhline(0, color='gray')
    fig_slope.canvas.manager.set_window_title(f'Slopes vs Simulation Amplitude')
    fig_int, ax_int = plt.subplots()
    ax_int.axhline(1, color='gray')
    fig_int.canvas.manager.set_window_title(f'Intercepts vs Simulation Amplitude')
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='amp')
        if data_sets_labels is None:
            lab = data_set
        else:
            lab = data_sets_labels[data_set]
        if data_sets_colors is None:
            ax_slope.errorbar(df_set['amp'], df_set['slope'], yerr=df_set['slope_err'], ls='none', marker='o',
                              label=lab)
            ax_int.errorbar(df_set['amp'], df_set['int'], yerr=df_set['int_err'], ls='none', marker='o', label=lab)
        else:
            ax_slope.errorbar(df_set['amp'], df_set['slope'], yerr=df_set['slope_err'], ls='none', marker='o',
                              color=data_sets_colors[data_set], label=lab)
            ax_int.errorbar(df_set['amp'], df_set['int'], yerr=df_set['int_err'], ls='none', marker='o',
                            color=data_sets_colors[data_set], label=lab)
    ax_slope.set_ylabel('Slope')
    ax_slope.set_xlabel('Simulation Amplitude')
    ax_slope.grid()
    ax_int.set_ylabel('Intercept')
    ax_int.set_xlabel('Simulation Amplitude')
    ax_int.grid()
    legend_slope = ax_slope.legend()
    # legend_slope.get_frame().set_alpha(0)
    legend_int = ax_int.legend()
    # legend_int.get_frame().set_alpha(0)
    fig_slope.tight_layout()
    fig_int.tight_layout()

    fig_slope_fit, ax_slope_fit = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig_slope_fit.canvas.manager.set_window_title(f'Slopes vs Simulation Amplitude Fits')
    ax_slope_fit.axhline(0, color='black')
    colors = iter(['black', 'blue', 'green', 'red', 'purple', 'salmon'])
    cl_type_name = {'_clmul_': 'Attractive', '_aclmul_': 'Repulsive'}
    spreads = pd.unique(df['spread'])

    for cl_type in cl_type_name.keys():
        spreads = sorted(spreads) if cl_type == '_aclmul_' else sorted(spreads, reverse=True)
        for spread in spreads:
            color = next(colors)
            df_spread = df[(df['spread'] == spread) & (df['name'].str.contains(cl_type))]
            print(f'{cl_type}, spread {spread}, size: {df_spread.size}')
            ax_slope_fit.errorbar(df_spread['amp'], df_spread['slope'], df_spread['slope_err'], ls='none', marker='o',
                                  color=color, label=f'{cl_type_name[cl_type]} σ={spread}')
            popt, pcov = cf(line, df_spread['amp'], df_spread['slope'], sigma=df_spread['slope_err'],
                            absolute_sigma=True)
            ax_slope_fit.plot(df_spread['amp'], line(df_spread['amp'], *popt), ls='--', color=color)
    ax_slope_fit.set_xlabel('Simulation Amplitude')
    ax_slope_fit.set_ylabel('Slope of Raw/Mix SD vs Total Protons per Event')
    ax_slope_fit.legend()
    ax_slope_fit.grid()
    fig_slope_fit.tight_layout()


def chi2_vs_protons(df, stat, div, cent, energy, data_type, data_sets_plt):
    chi2_sets = []
    df = df[(df['divs'] == div) & ((df['energy'] == energy) | (df['energy'] == 'sim'))]
    if 'data_type' in df:
        df = df[df['data_type'] == data_type]
    if 'cent' in df:
        df = df[df['cent'] == cent]
    if 'stat' in df:
        df = df[df['stat'] == stat]

    df_data = df[~df['name'].str.contains('sim_')]
    df_sim = df[df['name'].str.contains('sim_')]

    data_sets = np.unique(df_data['name'])
    sim_sets = np.unique(df_sim['name'])
    for data_set in data_sets:
        if data_set not in data_sets_plt:
            continue
        df_data_set = df_data[df_data['name'] == data_set]
        df_data_set.sort_values(by=['total_protons'])
        data_vals, data_errs = np.array(df_data_set['val']), np.array(df_data_set['err'])

        for sim_set in sim_sets:
            df_sim_set = df_sim[(df_sim['name'] == sim_set)]
            # Sim should have all possible tproton values, so just filter out those that aren't in data
            df_sim_set = df_sim_set[df_sim_set['total_protons'].isin(df_data_set['total_protons'])]
            df_sim_set.sort_values(by=['total_protons'])
            if list(df_data_set['total_protons']) != list(df_sim_set['total_protons']):
                print('Data and sim don\'t match total protons!')

            sim_vals, sim_errs = np.array(df_sim_set['val']), np.array(df_sim_set['err'])
            chi2 = np.sum((data_vals - sim_vals) ** 2 / (data_errs ** 2 + sim_errs ** 2))

            chi2_sets.append({'data_name': data_set, 'sim_name': sim_set, 'chi2': chi2, 'n': len(data_vals),
                              'divs': div, 'amp': df_sim_set['amp'].iloc[0], 'spread': df_sim_set['spread'].iloc[0],
                              'energy': energy, 'data_type': data_type})

    return chi2_sets


def plot_chi2_protons(df, energy, n_sims_plt=6):
    sim_sets_plt = {}
    chi_list = []
    for data_set in np.unique(df['data_name']):
        chi2_sums = []
        df_data = df[df['data_name'] == data_set]
        for sim_set in np.unique(df_data['sim_name']):
            df_sim = df_data[df_data['sim_name'] == sim_set]
            chi2_sum = df_sim['chi2'].sum()
            chi2_sums.append({'sim_name': sim_set, 'chi2_sum': chi2_sum, 'amp': df_sim['amp'].iloc[0],
                              'spread': df_sim['spread'].iloc[0], 'data_set': data_set, 'energy': energy})
        chi_list.extend(chi2_sums)
        chi2_sums = pd.DataFrame(chi2_sums)
        chi2_sums = chi2_sums.sort_values(by='chi2_sum').head(n_sims_plt)
        sim_sets_plt.update({data_set: chi2_sums['sim_name']})

    chi_res = []
    for data_set in np.unique(df['data_name']):
        fig, ax = plt.subplots()
        ax.set_title(data_set)
        ax.set_ylabel('Chi^2')
        ax.set_xlabel('Bin Width')
        df_data = df[df['data_name'] == data_set]
        for sim_set in sim_sets_plt[data_set]:
            df_sim = df_data[df_data['sim_name'] == sim_set].sort_values(by='divs')
            ax.axhline(0, ls='--', color='black')
            ax.plot(df_sim['divs'], df_sim['chi2'], marker='o', alpha=0.8, label=f'{sim_set}')
        ax.legend()
        fig.tight_layout()

    return sim_sets_plt, chi_res, chi_list


def plot_chi_sim_pars(chi_res):
    fig, ax = plt.subplots()
    ax.set_xlabel('Amp')
    ax.set_ylabel('Spread')
    for data_set in np.unique(chi_res['data_set']):
        df_set = chi_res[chi_res['data_set'] == data_set]
        ax.scatter(df_set['amp'], df_set['spread'], label=data_set)
    for x, y, e in zip(chi_res['amp'], chi_res['spread'], chi_res['energy']):
        ax.annotate(f'{e}GeV', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
    ax.legend()
    ax.grid()
    fig.tight_layout()


def stat_vs_divs(df, stat, total_protons, cent, energies, data_types, data_sets_plt, exclude_divs=[], y_ranges=None,
                 plot=False, fit=False):
    data = []
    for data_type in data_types:
        for data_set in data_sets_plt:
            for energy in energies:
                df_pre = df
                if 'data_type' in df_pre:
                    df_pre = df_pre[df_pre['data_type'] == data_type]
                if 'cent' in df_pre:
                    df_pre = df_pre[df_pre['cent'] == cent]
                if 'stat' in df_pre:
                    df_pre = df_pre[df_pre['stat'] == stat]
                df_set = df_pre[
                    (df_pre['name'] == data_set) & (df_pre['total_protons'] == total_protons) & (
                            df_pre['energy'] == energy)]
                df_set = df_set.sort_values(by=['divs'])
                if len(df_set) == 0:
                    continue
                if energy == 'sim':
                    lab = f'{data_set}_{data_type}'
                else:
                    lab = f'{data_set}_{data_type}_{energy}GeV'
                df_set = df_set.sort_values(by=['total_protons'])
                data.append((df_set, lab, data_set, energy))

    if plot:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {total_protons} protons')
        ax.set_xlabel('Division Width')
        ax.axhline(1, ls='--', color='black')
        if y_ranges:
            ax.set_ylim(y_ranges[stat])
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))

    fit_pars = []
    for df, lab, data_set, energy in data:
        c = next(color)
        if plot:
            ax.errorbar(df['divs'], df['val'], df['err'], label=lab, marker='o', ls='', color=c)
            if 'sys' in df:
                ax.errorbar(df['divs'], df['val'], df['sys'], marker='', ls='', elinewidth=3, color=c, alpha=0.3)

        if fit and df.size > 1:
            df = df[~df.divs.isin(exclude_divs)]
            popt, pcov = cf(quad_180, df['divs'], df['val'], sigma=df['err'], absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            fit_pars.append({'data_set': data_set, 'energy': energy, 'curvature': popt[0], 'baseline': popt[1],
                             'spread': df['spread'].iloc[0], 'amp': df['amp'].iloc[0],
                             'curve_err': perr[0], 'base_err': perr[1], 'color': c})
            x = np.linspace(0, 360, 100)
            if plot:
                ax.plot(x, quad_180(x, *popt), ls='--', color=c)

    if plot:
        ax.legend()
        fig.tight_layout()

    return pd.DataFrame(fit_pars)


def interp_vs_div(df_pars, stat, quad_par_line=line_yint1, r_amp_line=line_yint0, plot=False, fit=False):
    base_curve_lin_fits = []
    spread_rs_amp = []

    df_data = df_pars[~df_pars['data_set'].str.contains('sim_')]
    df_sim = df_pars[df_pars['data_set'].str.contains('sim_')]
    spreads = np.unique(df_sim['spread'])
    data_sets = np.unique(df_data['data_set'])

    if plot:
        fig, ax = plt.subplots()
        ax.set_title(stat)
        ax.set_xlabel('baseline')
        ax.set_ylabel('curvature')
        fix_dev, ax_dev = plt.subplots()
        ax_dev.set_title(f'{stat} Deviation from linear')
        ax_dev.set_xlabel('baseline')
        ax_dev.set_ylabel('curvature deviation from fit')
        ax_dev.axhline(0, ls='--')
        fig_rs, ax_rs = plt.subplots()
        ax_rs.set_title('Rs vs Amp')
        ax_rs.set_xlabel('amp')
        ax_rs.set_ylabel('r')
        fig_ramp_slope, ax_ramp_slope = plt.subplots()
        ax_ramp_slope.set_title('Rs vs Amp Slopes')
        ax_ramp_slope.set_xlabel('Spread')
        ax_ramp_slope.set_ylabel('r_amp slope')
        # color = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads) + len(data_sets))))
    spreads_list = []
    ramp_slopes = []
    ramp_slope_errs = []
    all_rs = []
    if plot:
        for data_set in data_sets:
            df_set = df_data[df_data['data_set'] == data_set]
            # c = next(color)
            ax.errorbar(df_set['baseline'], df_set['curvature'], xerr=df_set['base_err'], yerr=df_set['curve_err'],
                        marker='o', color=df_set['color'], label=data_set, ls='none', alpha=0.6)
            for x, y, e in zip(df_data['baseline'], df_data['curvature'], df_data['energy']):
                ax.annotate(f'{e}GeV', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')

    for spread in spreads:
        # if plot:
        #     c = next(color)
        df_spread = df_sim[df_sim['spread'] == spread]
        rs = np.sqrt(np.power(df_spread['curvature'], 2) + np.power((df_spread['baseline'] - 1), 2))
        all_rs += [r for r in rs]
        if len(df_spread) > 1:
            popt, pcov = cf(quad_par_line, df_spread['baseline'], df_spread['curvature'], sigma=df_spread['curve_err'],
                            absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            if type(spread) == str:
                spread_float = float(f'0.{spread.strip("spread")}') * 10
            else:
                spread_float = float(spread)
            # New column of df equal to sqrt(curvature**2 + (baseline - 1)**2)
            if len(popt) == 1:  # Ideally just append this to original df
                base_curve_lin_fits.append({'spread': spread_float, 'slope': popt[0], 'slope_err': perr[0]})
            elif len(popt) == 2:
                base_curve_lin_fits.append({'spread': spread_float, 'slope': popt[0], 'slope_err': perr[0],
                                            'int': popt[1],
                                            'int_err': perr[1]})  # Ideally just append this to original df
            amp_floats = [float(f'0.{amp.strip("amp")}') if type(amp) == str else float(amp)
                          for amp in df_spread['amp']]
            popt1, pcov1 = cf(r_amp_line, amp_floats, rs)
            perr1 = np.sqrt(np.diag(pcov1))
            spreads_list.append(spread_float)
            ramp_slopes.append(popt1[0])
            ramp_slope_errs.append(perr1[0])
            if plot:
                x_plt = np.linspace(min(df_sim['baseline']), 1, 3)
                ax.plot(x_plt, quad_par_line(x_plt, *popt), ls='--', alpha=0.7, color=df_spread['color'])
                ax_dev.scatter(df_spread['baseline'],
                               df_spread['curvature'] - line_yint1(df_spread['baseline'], *popt),
                               color=df_spread['color'])
                x_plt1 = np.array([0, 0.15])
                ax_rs.plot(x_plt1, r_amp_line(x_plt1, *popt1), ls='--', alpha=0.8, color=df_spread['color'])
                ax_rs.scatter(amp_floats, rs, color=df_spread['color'], marker='o', alpha=0.6, label=spread)
        if plot:
            ax.errorbar(df_spread['baseline'], df_spread['curvature'], xerr=df_spread['base_err'], marker='o',
                        yerr=df_spread['curve_err'], color=df_spread['color'], label=spread, ls='none', alpha=0.6)
            # if spread == 'ampt':
            #     # print(rs)
            #     ax_rs.axhline(rs[0], color=c, ls=':', label=spread, alpha=0.7, zorder=0)
            #     ampt_slope = Measure(df_spread['curvature'][0], df_spread['curve_err'][0]) / \
            #                  Measure(df_spread['baseline'][0], df_spread['base_err'][0])

    base_curve_lin_fits = pd.DataFrame(base_curve_lin_fits)

    if plot and len(spreads_list) > 0:
        interp = interp1d(spreads_list, ramp_slopes, kind='cubic')
        x_plt = np.linspace(min(spreads_list), max(spreads_list), 1000)
        ax_ramp_slope.plot(x_plt, interp(x_plt), ls='--', color='salmon')
        ax_ramp_slope.errorbar(spreads_list, ramp_slopes, ramp_slope_errs, ls='none', marker='o', alpha=0.9)
        ax.legend()
        ax.set_ylim(top=1.05 * max(df_sim['curvature']))
        ax_rs.legend()
        ax.grid()
        ax_ramp_slope.grid()
        ax_rs.grid()
        ax_rs.set_ylim(top=1.05 * max(all_rs))
        fig.tight_layout()
        fig_rs.tight_layout()
        fig_ramp_slope.tight_layout()

    return base_curve_lin_fits


def plot_fits(df_fits):
    for stat, df in df_fits.items():
        if 'int' in df:
            fig_int, ax_int = plt.subplots()
            ax_int.errorbar(df['spread'], df['int'], df['int_err'], ls='none')
            ax_int.set_title(f'{stat} Intercepts')
            ax_int.set_xlabel('Spread')
            ax_int.set_ylabel('Intercept')

        fig_slope, ax_slope = plt.subplots()
        print(df)
        interp = interp1d(df['spread'], df['slope'], kind='cubic')
        x_plt = np.linspace(min(df['spread']), max(df['spread']), 1000)
        ax_slope.plot(x_plt, interp(x_plt), ls='--', color='salmon')
        ax_slope.errorbar(df['spread'], df['slope'], df['slope_err'], ls='none', marker='o', alpha=0.9)
        ax_slope.set_title(f'{stat} Slopes')
        ax_slope.set_xlabel('Spread')
        ax_slope.set_ylabel('Slope')
        ax_slope.grid()
        fig_slope.tight_layout()


def plot_base_zeros(df, data_sets_plt, data_sets_labels, data_sets_colors):
    fig_base_amp, ax_base_amp = plt.subplots()
    fig_amp_base, ax_amp_base = plt.subplots()

    ax_base_amp.axvline(0, color='black')
    ax_amp_base.axhline(0, color='black')

    ax_base_amp.grid()
    ax_amp_base.grid()

    for data_set in data_sets_plt:
        df_set = df[df['data_set'].str.contains(data_set)]
        print(f'{data_set}: {df_set.size}')
        ax_base_amp.errorbar(df_set['zero_mag'], df_set['baseline'], yerr=df_set['base_err'], marker='o', ls='none',
                             xerr=df_set['zero_mag_err'], color=data_sets_colors[data_set],
                             label=data_sets_labels[data_set])
        ax_amp_base.errorbar(df_set['baseline'], df_set['zero_mag'], xerr=df_set['base_err'], marker='o', ls='none',
                             yerr=df_set['zero_mag_err'], color=data_sets_colors[data_set],
                             label=data_sets_labels[data_set])

    ax_base_amp.set_xlabel('Zeros')
    ax_base_amp.set_ylabel('Baseline')
    ax_amp_base.set_ylabel('Baseline')
    ax_amp_base.set_ylabel('Zeros')

    ax_base_amp.set_xlim(110, 230)
    ax_amp_base.set_ylim(110, 230)

    ax_base_amp.legend()
    ax_amp_base.legend()

    fig_base_amp.tight_layout()
    fig_base_amp.canvas.manager.set_window_title('Baselines vs Zeros')
    fig_amp_base.tight_layout()
    fig_amp_base.canvas.manager.set_window_title('Zeros vs Baselines')


def flow_vs_v2(df, div, res, out_dir=None):
    df = df[df['name'].str.contains(f'res{res}')]
    df = df[(df['stat'] == 'standard deviation') & (df['divs'] == div) & (df['data_type'] == 'raw')]
    v2s, slope_vals, slope_errs = [], [], []
    p = div / 360.0
    q = 1 - p
    for data_set in pd.unique(df['name']):
        v2 = float('0.' + data_set[data_set.find('_v2') + 3:])
        df_v2 = df[df['name'] == data_set]
        fig, ax = plt.subplots(figsize=(6.66, 5), dpi=144)
        ax.axhline(1, color='black', ls='--')
        raw_div_binom = [Measure(val, err) / np.sqrt(n * p * q) for n, val, err in
                         zip(df_v2['total_protons'], df_v2['val'], df_v2['err'])]
        y_vals, y_errs = [x.val for x in raw_div_binom], [x.err for x in raw_div_binom]
        ax.errorbar(df_v2['total_protons'], y_vals, y_errs, marker='o', ls='none', color='blue')
        ax.set_xlabel('Total Protons in Event')
        ax.set_ylabel('Raw / Binomial')
        ax.set_title(data_set)
        popt, pcov = cf(line_yint1, df_v2['total_protons'], y_vals, sigma=y_errs, absolute_sigma=True)
        ax.plot(df_v2['total_protons'], line_yint1(df_v2['total_protons'], *popt), ls='--', color='red')
        v2s.append(v2)
        slope_vals.append(popt[0])
        slope_errs.append(np.sqrt(np.diag(pcov))[0])

    v2s, slope_vals, slope_errs = zip(*sorted(zip(v2s, slope_vals, slope_errs)))
    v2s, slope_vals, slope_errs = np.array(v2s), np.array(slope_vals), np.array(slope_errs)
    fig, ax = plt.subplots(figsize=(6.66, 5), dpi=144)
    ax.set_xlabel('v2')
    ax.set_ylabel(r'$\sigma_{raw}/\sigma_{binomial}$ vs Total Protons Slope')
    ax.axhline(0, color='black')
    ax.axvline(0, color='black')
    ax.errorbar(v2s, slope_vals, slope_errs, marker='o', ls='none', alpha=0.8)
    popt, pcov = cf(quad, v2s, slope_vals, sigma=slope_errs, absolute_sigma=True)
    print('With constant and linear terms: ', [f'{val} +- {err}' for val, err in zip(popt, np.sqrt(np.diag(pcov)))])
    popt2, pcov2 = cf(x2, v2s, slope_vals, sigma=slope_errs, absolute_sigma=True)
    print('Just quadratic term: ', [f'{val} +- {err}' for val, err in zip(popt2, np.sqrt(np.diag(pcov2)))])
    xs = np.linspace(0, 0.07, 200)
    ax.plot(xs, quad(xs, *popt), ls='--', color='red', label='3 Parameter')
    ax.plot(xs, x2(xs, *popt2), ls='--', color='green', label='1 Parameter')
    ax.legend()
    ax.grid()
    fig.tight_layout()
    if out_dir:
        fig.savefig(f'{out_dir}slope_vs_v2_div_{div}.png')
        with open(f'{out_dir}slope_vs_v2_data_div_{div}.txt', 'w') as file:
            file.write(f'y=a*x^2 fit:\t{popt2[0]} {np.sqrt(pcov2[0][0])}\n')
            fit_string = '\t'.join([f'{val} {err}' for val, err in zip(popt, np.sqrt(np.diag(pcov)))])
            file.write(f'y=a*x^2+bx+c fit:\t{fit_string}\n\n')
            file.write('v2\tslope slope_err\n')
            for v2, slope_val, slope_err in zip(v2s, slope_vals, slope_errs):
                file.write(f'{v2}\t{slope_val} {slope_err}\n')


def subtract_avgs(df_a, df_b, val_col='avg', err_col='avg_err', meas_col='avg_meas'):
    if meas_col in df_a.columns:
        index_cols = list(set(df_a.columns) - {val_col, err_col, meas_col})
        df_sub = df_a.set_index(index_cols) - df_b.set_index(index_cols)
        df_sub.loc[:, err_col] = df_sub[meas_col].apply(lambda x: x.err)
    else:
        index_cols = list(set(df_a.columns) - {val_col, err_col})
        df_sub = df_a.set_index(index_cols) - df_b.set_index(index_cols)
        df_sub.loc[:, err_col] = np.sqrt(df_a.set_index(index_cols)[err_col] ** 2 +
                                         df_b.set_index(index_cols)[err_col] ** 2)

    return df_sub.reset_index()


def read_v2_slope_coefs(in_dir):
    coefs = {}
    files = os.listdir(in_dir)
    for file_name in files:
        if 'slope_vs_v2_data_' in file_name:
            div = int(file_name.split('.')[0].split('_')[-1])
            with open(f'{in_dir}{file_name}', 'r') as file:
                coef = file.readlines()[0].split('\t')[-1].split(' ')
                coef_val, coef_err = float(coef[0]), float(coef[1])
            coefs.update({div: Measure(coef_val, coef_err)})
    return coefs


def read_v2_values(in_dir, v2_type='v2'):
    v2s = {}
    dir_names = os.listdir(in_dir)
    for dir_name in dir_names:
        energy = int(dir_name.strip('GeV'))
        v2s.update({energy: {}})
        file_path = f'{in_dir}{dir_name}/{v2_type}.txt'
        if not os.path.isfile(file_path):
            continue
        with open(file_path, 'r') as file:
            lines = file.readlines()[1:]  # Skip column headers
            for line in lines:
                line = line.split('\t')
                if len(line) > 2:
                    cent = int(line[0])
                    line = line[1].split(' ')
                    v2_val = float(line[0])
                    v2_err = float(line[1])
                    v2s[energy].update({cent: Measure(v2_val, v2_err)})
    return v2s


def ampt_v2_closure_sub(slope_df, data_set_name, new_name, v2, coef, cent):
    df_new = slope_df[slope_df['name'] == data_set_name]
    df_new = df_new.assign(name=new_name)
    new_slope = [slope - coef * v2[energy][cent] ** 2 for energy, slope in zip(df_new['energy'], df_new['slope_meas'])]
    new_slope_vals, new_slope_errs = list(zip(*[(slope.val, slope.err) for slope in new_slope]))
    df_new = df_new.assign(slope=new_slope_vals)
    df_new = df_new.assign(slope_err=new_slope_errs)

    return pd.concat([slope_df, df_new], ignore_index=True)


def map_to_sim(baseline, zeros, interpolations):
    """
    Given baseline and zeros of quadratic fit to Raw SD / Mix SD Slope vs Partition Width, map to simulation amp, spread
    parameters with input interpolations
    :param baseline: Constant term in quadratic fit of Raw SD / Mix SD Slope vs Partition Width
    :param zeros: Reparameterized curvature term in quadratic fit of Raw SD / Mix SD Slope vs Partition Width
    :param interpolations: Interpolations from simulations for converting baseline and zeros to sim amp and spread
    :return: Interpolated simulation amplitude and spread, mapped from input baseline and zeros
    """
    if baseline > 0:
        interpolations = interpolations['_clmul_']
    else:
        interpolations = interpolations['_aclmul_']
    spread = interpolations['spread_zero'][zeros]
    amp_slope = interpolations['amp_slope_spread'][spread]
    amp = amp_slope * baseline

    return amp, spread


class MyBounds:
    def __init__(self, xmax=[0, 0.06], xmin=[0, 4.2]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


if __name__ == '__main__':
    main()
