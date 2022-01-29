#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 18 5:16 PM 2021
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slices

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
from mpl_toolkits import mplot3d
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit as cf
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
import scipy.interpolate as interp
from scipy.interpolate import RectBivariateSpline
from scipy.stats import norm

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats
from Measure import Measure


def main():
    threads = 15
    base_path = 'D:/Research/Results/Azimuth_Analysis/'
    # df_name = 'binom_slice_sds_cent8.csv'
    df_name = 'binom_slice_stats_cent8.csv'
    chi_df_name = 'chi_df_ampt_sds_cent8.csv'
    df_path = base_path + df_name
    # df_path = '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'D:/Research/Results/Azimuth_Analysis/binom_slice_cent_sds_df.csv'
    # df_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'C:/Users/Dyn04/Desktop/binom_slice_df.csv'
    # df_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/binom_slice_cent_sds_df.csv'
    sim_sets = []
    amps = []  # ['2']  # ['01', '02', '03', '04', '05', '06']
    spreads = []  # ['2', '25', '3', '35']  # ['02', '05', '1', '15', '2', '25', '3', '35', '4', '45']  # ['1']  # ['02', '05', '1', '15', '2', '25']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')

    # sim_amp_pairs = [('02', '1'), ('07', '225'), ('25', '3'), ('12', '45')]
    sim_amp_pairs = [('5', '35'), ('02', '05'), ('015', '1')]
    for amp, spread in sim_amp_pairs:
        sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')

    # spreads.extend(['5', '6'])
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}_test')
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.25), 'skewness': (0, 1.25), 'kurtosis': (-3, 2),
                'non-excess kurtosis': (0.95, 1.025)}
    stat_plot = 'non-excess kurtosis'  # 'standard deviation'  # ['standard deviation']  # , 'skewness', 'non-excess kurtosis']
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    total_protons_plt = 20
    cent_plt = 8
    energies_plt = [62, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = [7, 11, 19, 27, 39, 62]  # , 11, 19, 27, 39, 62]
    energy_plt = 62
    data_types_plt = ['divide']
    data_type_plt = 'divide'
    data_sets_plt = ['ampt_resample_def']  # ['ampt_resample_def']
    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    # df_data = df[~df['name'].str.contains('sim_')]
    # df_sim = df[df['name'].str.contains('sim_')]
    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    # df = df[df['spread'] != 1.2]

    # protons_fits = stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt, all_sets_plt,
    #                                plot=False, fit=True)
    # plot_protons_fits(protons_fits)

    # protons_fits = []
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_plt, data_types_plt, all_sets_plt,
    #                                        plot=False, fit=True)
    #     protons_fits.append(protons_fits_div)
    # protons_fits = pd.concat(protons_fits, ignore_index=True)
    # print(protons_fits)
    # plot_protons_fits_divs(protons_fits, data_sets_plt)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt, plot=True, fit=True)
    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, ['raw', 'mix'], all_sets_plt, plot=True, fit=False)
    stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt, all_sets_plt, plot=True, fit=False)
    plt.show()
    return

    chi_res = []
    chi2_list_all = []
    for energy in energies_fit:
        chi2_sets = []
        jobs = []
        print(f'Energy {energy}GeV')
        for data_type in data_types_plt:
            for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
                jobs.append((df, stat_plot, div, cent_plt, energy, data_type, data_sets_plt))
        with Pool(threads) as pool:
            for chi2_set in tqdm.tqdm(pool.istarmap(chi2_vs_protons, jobs), total=len(jobs)):
                chi2_sets.extend(chi2_set)
        # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        #     print(f'{energy}GeV {div} bin width')
        #     chi2_sets.extend(chi2_vs_protons(df, stat_plot, div, cent_plt, energy, data_type_plt, data_sets_plt))
        data_sim_sets_plt, chi_res_e, chi_list = plot_chi2_protons(pd.DataFrame(chi2_sets), energy, n_sims_plt=4)
        chi2_list_all.extend(chi_list)
        chi_res.extend(chi_res_e)

    chi2_sets_df = pd.DataFrame(chi2_list_all)
    chi2_sets_df.to_csv(base_path + chi_df_name, index=False)

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
    # df_fits = interp_vs_div(df_pars, stat_plot, quad_par_line=line_xint1, r_amp_line=line_yint0, plot=True)
    # for div_plt in [240]:
    # plot_vs_protons(df, stats_plot, div_plt, cent_plt, data_type_plt, y_ranges, data_sets_plt)
    # for total_protons_plt in [20]:
    #     df_pars = fit_vs_divs(df, stats_plot, total_protons_plt, cent_plt, energies_plt, data_type_plt, y_ranges,
    #                           all_sets_plt, exclude_divs, False)
    #     df_fits = interp_vs_div(df_pars, y_ranges, quad_par_line=line_xint1, r_amp_line=line_yint0, plot=True)
    #     # plot_fits(df_fits)
    plt.show()

    print('donzo')


def line(x, a, b):
    return a * x + b


def line_xint1(x, a):
    return a * (x - 1)


def line_yint0(x, a):
    return a * x


def quad_180(x, a, c):
    return a * (x - 180) ** 2 + c


def stat_vs_protons(df, stat, div, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False, fit=False,
                    hist=False):
    data = []
    for data_type in data_types:
        for data_set in data_sets_plt:
            for energy in energies:
                if 'data_type' in df:
                    df = df[df['data_type'] == data_type]
                if 'cent' in df:
                    df = df[df['cent'] == cent]
                if 'stat' in df:
                    df = df[df['stat'] == stat]
                df_set = df[(df['name'] == data_set) & (df['divs'] == div) & (df['energy'] == energy)]
                if len(df_set) == 0:
                    continue
                if energy == 'sim':
                    lab = f'{data_set}_{data_type}'
                else:
                    lab = f'{data_set}_{data_type}_{energy}GeV'
                df_set = df_set.sort_values(by=['total_protons'])
                data.append((df_set, lab, data_set, df_set['amp'].iloc[0], df_set['spread'].iloc[0]))

    if plot:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {div}° Divisions')
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--', color='black')
        if y_ranges:
            ax.set_ylim(y_ranges[stat])
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))

    fits = []
    for df, lab, data_set, amp, spread in data:
        if plot:
            c = next(color)
            if 'sim_' in data_set:
                ax.fill_between(df['total_protons'], df['val'] - df['err'], df['val'] + df['err'], label=lab, color=c,
                                alpha=0.4)
            else:
                ax.errorbar(df['total_protons'], df['val'], df['err'], label=lab,
                            marker='o', ls='', color=c, alpha=0.7)
                if 'sys' in df:
                    ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.4)

        if fit and len(df) > 1:
            popt, pcov = cf(line, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
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
        ax.legend()
        fig.tight_layout()

    return pd.DataFrame(fits)


def plot_protons_fits_divs(df, data_sets_plt):
    fig, ax = plt.subplots()
    for data_set in data_sets_plt:
        df_set = df[df['name'] == data_set]
        df_set.sort_values(by='divs')
        ax.errorbar(df_set['divs'], df_set['slope'], yerr=df_set['slope_err'], label=data_set)
    ax.set_ylabel('Slope')
    ax.set_xlabel('Bin Width')
    fig.tight_layout()


def plot_protons_fits(df):
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
            chi2 = np.sum((data_vals - sim_vals)**2 / (data_errs**2 + sim_errs**2))

            chi2_sets.append({'data_name': data_set, 'sim_name': sim_set, 'chi2': chi2, 'n': len(data_vals),
                              'divs': div, 'amp': df_sim_set['amp'].iloc[0], 'spread': df_sim_set['spread'].iloc[0],
                              'energy': energy, 'data_type': data_type})

    return chi2_sets


def plot_chi2_protons(df, energy, n_sims_plt=6, df_write_path=None):
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
        # chi_df.update({data_set: chi2_sums, 'energy': energy})
        chi2_sums = chi2_sums.sort_values(by='chi2_sum').head(n_sims_plt)
        sim_sets_plt.update({data_set: chi2_sums['sim_name']})
    chi_df = pd.DataFrame(chi_list)
    # chi_df.to_csv(df_write_path)

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
    #
    #     df_chi = chi_df[chi_df['data_set'] == data_set]
    #     df_chi = df_chi.sort_values(by=['spread', 'amp'])
    #     fig_chisum, ax_chisum = plt.subplots()
    #
    #     color = iter(plt.cm.rainbow(np.linspace(0, 1, len(np.unique(df_chi['spread'])))))
    #     ax_chisum.grid()
    #     ax_chisum.axhline(0, color='black', ls='--')
    #     for spread in np.unique(df_chi['spread']):
    #         c = next(color)
    #         df_s = df_chi[df_chi['spread'] == spread].sort_values(by='amp')
    #         ax_chisum.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), marker='o', color=c, label=spread)
    #         # ax_chisum.plot(df_s['amp'], df_s['chi2_sum'], marker='o', color=c, alpha=0.8)
    #         f_1d = interp1d(df_s['amp'], np.log10(df_s['chi2_sum']), kind='cubic')
    #         x_interp = np.linspace(min(df_s['amp']), max(df_s['amp']), 1000)
    #         ax_chisum.plot(x_interp, f_1d(x_interp), color=c, alpha=0.8)
    #     ax_chisum.set_title(f'{data_set} {energy}GeV')
    #     ax_chisum.set_xlabel('amp')
    #     ax_chisum.set_ylabel('log10 chi2 sum')
    #     ax_chisum.legend()
    #     fig_chisum.tight_layout()
    #
    #     df_rect = df_chi[(df_chi['amp'] != 1.5) & (df_chi['spread'] != 0) & (df_chi['spread'] != 12) &
    #                      (df_chi['spread'] != 4)]
    #     # x_rect = np.unique(df_rect['amp'])
    #     # y_rect = np.unique(df_rect['spread'])
    #     # z_rect = []
    #     # for xi_rect in x_rect:
    #     #     z_rect.append([])
    #     #     for yi_rect in y_rect:
    #     #         z_rect[-1].append(df_rect[(df_rect['amp'] == xi_rect) &
    #     #                                   (df_rect['spread'] == yi_rect)].iloc[0]['chi2_sum'])
    #     # z_rect = np.asarray(z_rect)
    #     # print(x_rect.shape, y_rect.shape, z_rect.shape)
    #     # f = RectBivariateSpline(x_rect, y_rect, z_rect)
    #
    #     x_in, y_in, z_in = np.array(df_chi['amp']), np.array(df_chi['spread']), np.array(df_chi['chi2_sum'])
    #     z_in_log = np.log10(z_in)
    #     # print(x_in.shape, y_in.shape, z_in.shape)
    #     # print(x_in, y_in, z_in, z_in_log)
    #
    #     tcks = interp.bisplrep(x_in, y_in, z_in_log, kx=1, ky=1)
    #
    #     # f = interp2d(df_chi['amp'], df_chi['spread'], df_chi['chi2_sum'], kind='cubic')
    #     x = np.linspace(df_chi['amp'].min(), df_chi['amp'].max(), 120)
    #     y = np.linspace(df_chi['spread'].min(), df_chi['spread'].max(), 100)
    #     xx, yy = np.meshgrid(x, y)
    #     # z = f(x, y)
    #     z = interp.bisplev(x, y, tcks)
    #     print(z)
    #     print(z.shape)
    #     print(interp.bisplev(0.025, 3.0, tcks))
    #     print(interp.bisplev(0.155, 1.8, tcks))
    #     # z_min = np.min(z)
    #     # z = z - z_min + 1  # np.where(z > 0, z, 0.001)
    #     print(data_set)
    #     x0 = [0.04, 2.1]
    #     # min_res = minimize(lambda x_opt: f(*x_opt) - z_min + 1, x0, bounds=[(0, 0.1), (0, None)])
    #     # print(min_res)
    #     # # print(*min_res.x, min_res.fun)
    #     # mybounds = MyBounds()
    #     # bas_res = basinhopping(lambda x_opt: f(*x_opt) - z_min + 1, x0, minimizer_kwargs={'method': 'L-BFGS-B'},
    #     #                        niter=100, accept_test=mybounds)
    #     # print(bas_res)
    #     # print('local_min: ', *min_res.x, min_res.fun)
    #     # print('basin min: ', *bas_res.x, bas_res.fun)
    #     # chi_res.append({'data_set': data_set, 'energy': energy, 'amp': bas_res.x[0], 'spread': bas_res.x[1]})
    #
    #     fig_chisum2, ax_chisum2 = plt.subplots()
    #
    #     color = iter(plt.cm.rainbow(np.linspace(0, 1, len(np.unique(df_chi['spread'])))))
    #     ax_chisum2.grid()
    #     ax_chisum2.axhline(0, color='black', ls='--')
    #     for spread in np.unique(df_chi['spread']):
    #         c = next(color)
    #         df_s = df_chi[df_chi['spread'] == spread].sort_values(by='amp')
    #         ax_chisum2.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), marker='o', color=c, label=spread)
    #         i = 0
    #         while y[i] < spread:
    #             i += 1
    #         ax_chisum2.plot(x, z.T[i], color=c, alpha=0.8)
    #     ax_chisum2.set_title(f'{data_set} {energy}GeV')
    #     ax_chisum2.set_xlabel('amp')
    #     ax_chisum2.set_ylabel('log10 chi2 sum')
    #     ax_chisum2.legend()
    #     fig_chisum2.tight_layout()
    #
    #     # fig_interp1d, ax_interp1d = plt.subplots()
    #     # norm = matplotlib.colors.Normalize(vmin=y.min(), vmax=y.max())
    #     # cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.plasma)
    #     # # cm = plt.cm.get_cmap('plasma')
    #     # # color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data))))
    #     # for spread in y[::5]:
    #     #     ax_interp1d.plot(x, f(x, spread), color=cmap.to_rgba(spread))
    #     # ax_interp1d.set_title(f'{data_set} {energy}GeV')
    #     # ax_interp1d.set_xlabel('amp')
    #     # ax_interp1d.set_ylabel('spread')
    #     # fig_interp1d.colorbar(cmap, ticks=np.arange(y.min(), y.max(), 1))
    #
    #     fig_2d, ax_2d = plt.subplots()
    #     im_obj = ax_2d.imshow(z.T, extent=[min(x), max(x), min(y), max(y)], aspect='auto', origin='lower',
    #                           cmap='jet')
    #     ax_2d.contour(z.T, origin='lower', cmap='plasma', extent=[min(x), max(x), min(y), max(y)])
    #     ax_2d.scatter(x_in, y_in, c=z_in_log, marker='o', cmap='jet')
    #     ax_2d.scatter(*x0, color='black', label='Local Start')
    #     # ax_2d.scatter(*min_res.x, color='orange', label='Local Min')
    #     # ax_2d.scatter(*bas_res.x, color='red', label='Global Min')
    #     ax_2d.set_title(f'{data_set} {energy}GeV')
    #     ax_2d.set_xlabel('amp')
    #     ax_2d.set_ylabel('spread')
    #     ax_2d.legend()
    #     fig_2d.colorbar(im_obj, ax=ax_2d)
    #     fig_2d.tight_layout()
    #
    #     fig_3d_scatter = plt.figure()
    #     ax_3d_scatter = plt.axes(projection='3d')
    #     ax_3d_scatter.scatter3D(x_in, y_in, z_in_log, c=z_in_log, cmap='viridis')
    #     ax_3d_scatter.set_xlabel('amp')
    #     ax_3d_scatter.set_ylabel('spread')
    #     ax_3d_scatter.set_zlabel('log chi2_sum')
    #
    #     fig_3d_interp = plt.figure()
    #     ax_3d_interp = plt.axes(projection='3d')
    #     ax_3d_interp.plot_surface(xx, yy, z.T, cmap='viridis', edgecolor='none')
    #     ax_3d_interp.scatter3D(x_in, y_in, z_in_log, color='black', linewidth=0.5, alpha=0.5)
    #     # ax_3d_interp.scatter(*x0, f(*x0), color='black', label='Start')
    #     # ax_3d_interp.scatter(*min_res.x, min_res.fun, color='orange', label='Local Min')
    #     # ax_3d_interp.scatter(*bas_res.x, bas_res.fun, color='red', label='Global Min')
    #     ax_3d_interp.set_xlabel('amp')
    #     ax_3d_interp.set_ylabel('spread')
    #     ax_3d_interp.set_zlabel('log chi2')
    #     # ax_3d_interp.set_zlim(0, 10000)
    #     ax_3d_interp.set_title(f'{data_set} {energy}GeV')
    #     # ax_3d_interp.legend()
    #     fig_3d_interp.tight_layout()
    #
    #     # fig_3d_interp_log = plt.figure()
    #     # ax_3d_interp_log = plt.axes(projection='3d')
    #     # ax_3d_interp_log.plot_surface(xx, yy, np.log10(z), cmap='viridis', edgecolor='none')
    #     # # ax_3d_interp_log.scatter(*x0, np.log10(f(*x0)), color='black', label='Start')
    #     # # ax_3d_interp_log.scatter(*min_res.x, np.log10(min_res.fun), color='orange', label='Local Min')
    #     # # ax_3d_interp_log.scatter(*bas_res.x, np.log10(bas_res.fun), color='red', label='Global Min')
    #     # ax_3d_interp_log.set_xlabel('amp')
    #     # ax_3d_interp_log.set_ylabel('spread')
    #     # ax_3d_interp_log.set_zlabel('log10 chi2')
    #     # ax_3d_interp_log.set_title(data_set)
    #     # # ax_3d_interp_log.legend()
    #     # fig_3d_interp_log.tight_layout()

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
                df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['energy'] == energy) &
                            (df['total_protons'] == total_protons) & (df['cent'] == cent) & (df['stat'] == stat)]
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
            ax.errorbar(df['divs'], df['val'], df['sys'], marker='', ls='', elinewidth=3, color=c, alpha=0.3)

        if fit and len(df) > 1:
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


def interp_vs_div(df_pars, stat, quad_par_line=line_xint1, r_amp_line=line_yint0, plot=False, fit=False):
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
                               df_spread['curvature'] - line_xint1(df_spread['baseline'], *popt),
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
