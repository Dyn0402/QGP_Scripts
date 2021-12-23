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
from mpl_toolkits import mplot3d
import pandas as pd
from scipy.optimize import curve_fit as cf
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

from DistStats import DistStats
from Measure import Measure


def main():
    # df_path = '/home/dylan/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'D:/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    # df_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/binom_slice_df.csv'
    df_path = 'C:/Users/Dyn04/Desktop/binom_slice_df.csv'
    sim_sets = []
    amps = ['01', '02', '03', '04', '05', '06']
    spreads = ['02', '05', '1', '15', '2', '25']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
    y_ranges = {'mean': (0.8, 1.2), 'standard deviation': (0.8, 1.25), 'skewness': (0, 1.25), 'kurtosis': (-3, 2),
                'non-excess kurtosis': (0.95, 1.025)}
    stat_plot = 'standard deviation'  # ['standard deviation']  # , 'skewness', 'non-excess kurtosis']
    div_plt = 120
    exclude_divs = [356]
    total_protons_plt = 20
    cent_plt = 8
    energies_plt = [62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energy_plt = 62
    data_types_plt = ['divide']
    data_type_plt = 'divide'
    data_sets_plt = ['bes_resample_def', 'ampt_resample_def']
    all_sets_plt = data_sets_plt + sim_sets

    df = pd.read_csv(df_path)
    df = df.dropna()
    # df_data = df[~df['name'].str.contains('sim_')]
    # df_sim = df[df['name'].str.contains('sim_')]
    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

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

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies_plt, data_types_plt, data_sets_plt, plot=True, fit=True)
    chi2_sets = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        chi2_sets.extend(chi2_vs_protons(df, stat_plot, div, cent_plt, energy_plt, data_type_plt, data_sets_plt))
    data_sim_sets_plt = plot_chi2_protons(pd.DataFrame(chi2_sets), n_sims_plt=4)
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


def stat_vs_protons(df, stat, div, cent, energies, data_types, data_sets_plt, y_ranges=None, plot=False, fit=False):
    data = []
    for data_type in data_types:
        for data_set in data_sets_plt:
            for energy in energies:
                df_set = df[(df['name'] == data_set) & (df['data_type'] == data_type) & (df['divs'] == div) &
                            (df['cent'] == cent) & (df['stat'] == stat) & (df['energy'] == energy)]
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
                ax.errorbar(df['total_protons'], df['val'], df['sys'], marker='', ls='', elinewidth=3,
                            color=c, alpha=0.4)

        if fit and len(df) > 1:
            popt, pcov = cf(line, df['total_protons'], df['val'], sigma=df['err'], absolute_sigma=True)
            fits.append({'name': data_set, 'divs': div, 'amp': amp, 'spread': spread, 'slope': popt[0],
                         'slope_err': np.sqrt(np.diag(pcov))[0]})
            if plot:
                ax.plot(df['total_protons'], line(df['total_protons'], *popt), ls='--', color=c)

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
    df_data = df[~df['name'].str.contains('sim_')]
    df_sim = df[df['name'].str.contains('sim_')]
    data_sets = np.unique(df_data['name'])
    sim_sets = np.unique(df_sim['name'])
    for data_set in data_sets:
        if data_set not in data_sets_plt:
            continue
        df_data_set = df_data[(df_data['name'] == data_set) & (df_data['data_type'] == data_type) &
                              (df_data['divs'] == div) & (df_data['cent'] == cent) & (df_data['stat'] == stat) &
                              (df_data['energy'] == energy)]
        df_data_set.sort_values(by=['total_protons'])
        data_meases = df_data_set.apply(lambda row: Measure(row['val'], row['err']), axis=1)
        for sim_set in sim_sets:
            df_sim_set = df_sim[(df_sim['name'] == sim_set) & (df_sim['data_type'] == data_type) &
                                (df_sim['divs'] == div) & (df_sim['cent'] == cent) & (df_sim['stat'] == stat) &
                                (df_sim['energy'] == 'sim')]
            # Sim should have all possible tproton values, so just filter out those that aren't in data
            df_sim_set = df_sim_set[df_sim_set['total_protons'].isin(df_data_set['total_protons'])]
            df_sim_set.sort_values(by=['total_protons'])
            # print('data_protons', df_data_set['total_protons'])
            # print('sim_protons', df_sim_set['total_protons'])
            if list(df_data_set['total_protons']) != list(df_sim_set['total_protons']):
                print('Data and sim don\'t match total protons!')
            sim_meases = df_sim_set.apply(lambda row: Measure(row['val'], row['err']), axis=1)
            # print(f'{data_set}, {sim_set}: #data: {len(data_meases)}, #sim: {len(sim_meases)}')

            # print('data_meas: ', type(data_meases), len(data_meases))
            # print('sim_meas: ', type(sim_meases), len(sim_meases))
            diff = np.array(data_meases) - np.array(sim_meases)
            # print(diff)
            chi2 = sum([(x.val / x.err) ** 2 for x in diff])
            chi2_sets.append({'data_name': data_set, 'sim_name': sim_set, 'chi2': chi2, 'n': len(diff), 'divs': div,
                              'amp': df_sim['amp'], 'spread': df_sim['spread']})

    return chi2_sets


def plot_chi2_protons(df, n_sims_plt=6):
    sim_sets_plt = {}
    chi_df = {}
    for data_set in np.unique(df['data_name']):
        chi2_sums = []
        df_data = df[df['data_name'] == data_set]
        for sim_set in np.unique(df_data['sim_name']):
            df_sim = df_data[df_data['sim_name'] == sim_set]
            chi2_sum = df_sim['chi2'].sum()
            chi2_sums.append({'sim_name': sim_set, 'chi2_sum': chi2_sum, 'amp': df_sim['amp'],
                              'spread': df_sim['spread']})
        chi2_sums = pd.DataFrame(chi2_sums)
        chi_df.update({data_set: chi2_sums})
        chi2_sums = chi2_sums.sort_values(by='chi2_sum').head(n_sims_plt)
        sim_sets_plt.update({data_set: chi2_sums['sim_name']})
    print(sim_sets_plt)

    for data_set in np.unique(df['data_name']):
        fig, ax = plt.subplots()
        ax.set_title(data_set)
        ax.set_ylabel('Chi^2')
        ax.set_xlabel('Bin Width')
        df_data = df[df['data_name'] == data_set]
        for sim_set in sim_sets_plt[data_set]:
            df_sim = df_data[df_data['sim_name'] == sim_set]
            ax.axhline(0, ls='--', color='black')
            ax.plot(df_sim['divs'], df_sim['chi2'], marker='o', alpha=0.8, label=f'{sim_set}')
        ax.legend()
        fig.tight_layout()

        df_chi = chi_df[data_set]
        f = interp2d(df_chi['amp'], df_chi['spread'], df_chi['chi2_sum'], kind='cubic')
        x = np.linspace(df_chi['amp'].min(), df_chi['amp'].max(), 100)
        y = np.linspace(df_chi['spread'].min(), df_chi['spread'].max(), 100)
        xx, yy = np.meshgrid(x, y)
        fig_3d_interp = plt.figure()
        ax_3d_interp = plt.axes(projection='3d')
        z = f(x, y)
        ax_3d_interp.plot_surface(xx, yy, z, cmap='viridis', edgecolor='none')
        ax_3d_interp.set_xlabel('amp')
        ax_3d_interp.set_ylabel('spread')
        ax_3d_interp.set_zlabel('chi2')

    return sim_sets_plt


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


if __name__ == '__main__':
    main()
