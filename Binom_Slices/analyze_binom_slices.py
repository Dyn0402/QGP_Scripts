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
import pandas as pd
from scipy.optimize import curve_fit as cf
from scipy.interpolate import interp1d

from DistStats import DistStats
from Measure import Measure


def main():
    base_path = '/home/dylan/Research/'
    df_path = '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    sim_sets = []
    amps = ['01', '02', '03', '04', '05', '06']
    spreads = ['02', '05', '1', '15', '2', '25']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.25), 'skewness': (0, 1.25), 'kurtosis': (-3, 2),
                'non-excess kurtosis': (0.95, 1.025)}
    stats_plot = ['standard deviation']  # ['standard deviation']  # , 'skewness', 'non-excess kurtosis']
    div_plt = 120
    total_protons_plt = 20
    cent_plt = 8
    energies_plt = [7]  #[7, 11, 19, 27, 39, 62]
    energy_plt = 62
    data_types_plt = ['divide']
    data_type_plt = 'divide'
    data_sets_plt = ['bes_resample_def', 'ampt_resample_def']  # , 'ampt_def_resample']  # 'sim_aclmul_amp01_spread2', 'sim_aclmul_amp02_spread2'
    all_sets_plt = data_sets_plt + sim_sets

    df = pd.read_csv(df_path)
    # print(df['name'])
    df = df.dropna()

    plot_vs_protons(df, stats_plot, div_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt)
    # plot_vs_protons(df, stats_plot, div_plt, cent_plt, energy_plt, ['raw', 'mix'], data_sets_plt)
    plot_vs_divs(df, stats_plot, total_protons_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt)
    # for div_plt in [240]:
    # plot_vs_protons(df, stats_plot, div_plt, cent_plt, data_type_plt, y_ranges, data_sets_plt)
    fit_vs_protons(df, stats_plot, div_plt, cent_plt, energies_plt, data_type_plt, y_ranges, data_sets_plt, True)
    for total_protons_plt in [20]:
        df_pars = fit_vs_divs(df, stats_plot, total_protons_plt, cent_plt, energies_plt, data_type_plt, y_ranges, all_sets_plt, False)
        df_fits = interp_vs_div(df_pars, y_ranges, quad_par_line=line_xint1, r_amp_line=line_yint0, plot=True)
        plot_fits(df_fits)
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


def plot_vs_protons(df, stat_names, div, cent, energies, data_types, data_set_plt, y_ranges=None):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {div}° Divisions')
        ax.set_xlabel('Total Protons in Event')
        ax.axhline(1, ls='--')
        if y_ranges:
            ax.set_ylim(y_ranges[stat])
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(data_set_plt) * len(data_types) * len(energies))))

        for data_type in data_types:
            for data_set in data_set_plt:
                for energy in energies:
                    c = next(color)
                    lab = f'{data_set}_{data_type}_{energy}GeV'
                    df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['div'] == div) &
                                (df['cent'] == cent) & (df['stat'] == stat) & (df['energy'] == energy)]
                    df_set = df_set.sort_values(by=['total_protons'])
                    if 'bes_' in data_set:
                        scale = 1.0 / 0.8
                    else:
                        scale = 1
                    ax.errorbar(df_set['total_protons'] * scale, df_set['val'], df_set['err'], label=lab,
                                marker='o', ls='', color=c, alpha=0.4)
                    ax.errorbar(df_set['total_protons'] * scale, df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.3)
                    # ax.fill_between(df_set['total_protons'], df_set['val'] - df_set['err'],
                    #                 df_set['val'] + df_set['err'], label=f'{data_set}_{data_type}', color=c)
        ax.legend()
        fig.tight_layout()


def plot_vs_divs(df, stat_names, total_protons, cent, energies, data_types, data_sets_plt, y_ranges=None):
    for stat in stat_names:
        fig, ax = plt.subplots()
        ax.set_title(f'{stat} {total_protons} protons')
        ax.set_xlabel('Division Width')
        ax.axhline(1, ls='--')
        if y_ranges:
            ax.set_ylim(y_ranges[stat])
        n_sets = len(data_sets_plt) * len(data_types) * len(energies)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, n_sets)))

        shift = 0
        for data_type in data_types:
            for data_set in data_sets_plt:
                for energy in energies:
                    c = next(color)
                    lab = f'{data_set}_{data_type}_{energy}GeV'
                    df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['energy'] == energy) &
                                (df['total_protons'] == total_protons) & (df['cent'] == cent) & (df['stat'] == stat)]
                    df_set = df_set.sort_values(by=['div'])
                    # if 'ampt' in data_set:
                    ax.errorbar(df_set['div'] + shift, df_set['val'], df_set['err'], label=lab,
                                marker='o', ls='', color=c)
                    ax.errorbar(df_set['div'] + shift, df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.3)
                    # elif data_set in data_sets_plt:
                    #     ax.fill_between(df_set['div'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                    #                     label=data_set, color=c)
                    shift += 0.2
        ax.legend()
        fig.tight_layout()

    # plt.show()


def fit_vs_protons(df, stat_names, div, cent, energies, data_type, y_ranges, data_set_plt, plot=False):
    for stat in stat_names:
        if plot:
            fig, ax = plt.subplots()
            ax.set_title(f'{stat} {div}° Divisions')
            ax.set_xlabel('Total Protons in Event')
            ax.axhline(1, ls='--')
            ax.set_ylim(y_ranges[stat])
            n_sets = len(data_set_plt) * len(energies)
            color = iter(plt.cm.rainbow(np.linspace(0, 1, n_sets)))

        for data_set in data_set_plt:
            for energy in energies:
                df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['div'] == div) &
                            (df['cent'] == cent) & (df['stat'] == stat) & (df['energy'] == energy)]
                df_set = df_set.sort_values(by=['total_protons'])
                lab = f'{data_set}_{energy}GeV'
                if plot:
                    c = next(color)
                    # if 'ampt' in data_set:
                    ax.errorbar(df_set['total_protons'], df_set['val'], df_set['err'], label=lab, marker='o', ls='',
                                color=c, alpha=0.4)
                    ax.errorbar(df_set['total_protons'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.2)
                    # elif data_set in data_set_plt:
                    #     ax.fill_between(df_set['total_protons'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                    #                     label=lab, color=c, alpha=0.8)
                if len(df_set['val']) > 1:
                    # print(f'df_set[total_protons]: {df_set["total_protons"]}, df_set[val]: {df_set["val"]}')
                    popt, pcov = cf(line, df_set['total_protons'], df_set['val'], sigma=df_set['err'], absolute_sigma=True)
                    if 'ampt' in data_set or data_set in data_set_plt and plot:
                        # print(df_set['total_protons'])
                        # print(df_set['val'])
                        # print(df_set['err'])
                        ax.plot(df_set['total_protons'], line(df_set['total_protons'], *popt), ls='--', color=c)
        if plot:
            ax.legend()
            fig.tight_layout()


def fit_vs_divs(df, stat_names, total_protons, cent, energies, data_type, y_ranges, data_set_plt, plot=False):
    df_pars = {}
    for stat in stat_names:
        if plot:
            fig, ax = plt.subplots()
            ax.set_title(f'{stat} {total_protons} protons')
            ax.set_xlabel('Division Width')
            ax.axhline(1, ls='--')
            ax.set_ylim(y_ranges[stat])

            n_sets = len(data_set_plt) + (sum(1 for x in data_set_plt if 'sim_' not in x) * len(energies))
            color = iter(plt.cm.rainbow(np.linspace(0, 1, n_sets)))
        fit_pars = []
        for data_set in data_set_plt:
            for energy in energies:
                df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['energy'] == energy) &
                            (df['total_protons'] == total_protons) & (df['cent'] == cent) & (df['stat'] == stat) &
                            (df['div'] != 356)]
                if len(df_set) == 0:
                    continue
                df_set = df_set.sort_values(by=['div'])
                lab = f'{data_set}_{energy}GeV'

                if plot:
                    c = next(color)
                    # if 'ampt' in data_set:
                    ax.errorbar(df_set['div'], df_set['val'], df_set['err'], label=lab, marker='o', ls='',
                                color=c, alpha=0.8)
                    ax.errorbar(df_set['div'], df_set['val'], df_set['sys'], marker='', ls='', elinewidth=3,
                                color=c, alpha=0.3)
                    # elif data_set in data_set_plt:
                    #     ax.fill_between(df_set['div'], df_set['val'] - df_set['err'], df_set['val'] + df_set['err'],
                    #                     label=lab, color=c, alpha=0.8)

                if len(df_set['val'] > 1):
                    popt, pcov = cf(quad_180, df_set['div'], df_set['val'], sigma=df_set['err'], absolute_sigma=True)
                    perr = np.sqrt(np.diag(pcov))
                    fit_pars.append({'data_set': data_set, 'energy': energy, 'curvature': popt[0], 'baseline': popt[1],
                                     'spread': df_set['spread'].iloc[0], 'amp': df_set['amp'].iloc[0],
                                     'curve_err': perr[0], 'base_err': perr[1]})
                    x = np.linspace(0, 360, 100)
                    if plot and ('ampt' in data_set or data_set in data_set_plt):
                        ax.plot(x, quad_180(x, *popt), ls='--', color=c)
        if plot:
            ax.legend()
            fig.tight_layout()

        df_pars.update({stat: pd.DataFrame(fit_pars)})

    return df_pars


def interp_vs_div(df_pars_stats, y_ranges, quad_par_line=line_xint1, r_amp_line=line_yint0, plot=False):
    base_curve_lin_fits = {}
    for stat, df_pars in df_pars_stats.items():
        fit_df = []
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
            color = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads) + len(data_sets))))
        spreads_list = []
        ramp_slopes = []
        ramp_slope_errs = []
        all_rs = []
        if plot:
            for data_set in data_sets:
                df_set = df_data[df_data['data_set'] == data_set]
                c = next(color)
                ax.errorbar(df_set['baseline'], df_set['curvature'], xerr=df_set['base_err'], yerr=df_set['curve_err'],
                            marker='o', color=c, label=data_set, ls='none', alpha=0.6)
                for x, y, e in zip(df_data['baseline'], df_data['curvature'], df_data['energy']):
                    ax.annotate(f'{e}GeV', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')

        for spread in spreads:
            if plot:
                c = next(color)
            df_spread = df_sim[df_sim['spread'] == spread]
            rs = np.sqrt(np.power(df_spread['curvature'], 2) + np.power((df_spread['baseline'] - 1), 2))
            all_rs += [r for r in rs]
            if len(df_spread['curvature']) > 1:
                popt, pcov = cf(quad_par_line, df_spread['baseline'], df_spread['curvature'], sigma=df_spread['curve_err'],
                                absolute_sigma=True)
                perr = np.sqrt(np.diag(pcov))
                if type(spread) == str:
                    spread_float = float(f'0.{spread.strip("spread")}') * 10
                else:
                    spread_float = float(spread)
                # New column of df equal to sqrt(curvature**2 + (baseline - 1)**2)
                if len(popt) == 1:  # Ideally just append this to original df
                    fit_df.append({'spread': spread_float, 'slope': popt[0], 'slope_err': perr[0]})
                elif len(popt) == 2:
                    fit_df.append({'spread': spread_float, 'slope': popt[0], 'slope_err': perr[0],
                                   'int': popt[1], 'int_err': perr[1]})  # Ideally just append this to original df
                amp_floats = [float(f'0.{amp.strip("amp")}') if type(amp) == str else float(amp)
                              for amp in df_spread['amp']]
                popt1, pcov1 = cf(r_amp_line, amp_floats, rs)
                perr1 = np.sqrt(np.diag(pcov1))
                spreads_list.append(spread_float)
                ramp_slopes.append(popt1[0])
                ramp_slope_errs.append(perr1[0])
                if plot:
                    x_plt = np.linspace(min(df_sim['baseline']), 1, 3)
                    ax.plot(x_plt, quad_par_line(x_plt, *popt), ls='--', alpha=0.7, color=c)
                    ax_dev.scatter(df_spread['baseline'],
                                   df_spread['curvature'] - line_xint1(df_spread['baseline'], *popt), color=c)
                    x_plt1 = np.array([0, 0.15])
                    ax_rs.plot(x_plt1, r_amp_line(x_plt1, *popt1), ls='--', alpha=0.8, color=c)
                    ax_rs.scatter(amp_floats, rs, color=c, marker='o', alpha=0.6, label=spread)
            if plot:
                ax.errorbar(df_spread['baseline'], df_spread['curvature'], xerr=df_spread['base_err'], marker='o',
                            yerr=df_spread['curve_err'], color=c, label=spread, ls='none', alpha=0.6)
                # if spread == 'ampt':
                #     # print(rs)
                #     ax_rs.axhline(rs[0], color=c, ls=':', label=spread, alpha=0.7, zorder=0)
                #     ampt_slope = Measure(df_spread['curvature'][0], df_spread['curve_err'][0]) / \
                #                  Measure(df_spread['baseline'][0], df_spread['base_err'][0])

        fit_df = pd.DataFrame(fit_df)

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
        base_curve_lin_fits.update({stat: fit_df})

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


if __name__ == '__main__':
    main()
