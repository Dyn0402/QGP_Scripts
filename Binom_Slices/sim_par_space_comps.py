#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 09 4:18 PM 2022
Created in PyCharm
Created as QGP_Scripts/sim_par_space_comps.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.interpolate as interp
from scipy.optimize import basinhopping
from scipy.optimize import curve_fit
import pandas as pd

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm


def main():
    chi2_sd_slices_bs()
    print('donzo')


class Pars:
    def __init__(self):
        self.base_path = 'F:/Research/Results/Azimuth_Analysis/'
        self.data_sets = ['bes', 'ampt_new_coal', 'cfev']
        self.data_sets_colors = dict(zip(self.data_sets, ['black', 'red', 'blue']))
        self.data_sets_labels = dict(zip(self.data_sets, ['STAR', 'AMPT', 'MUSIC+FIST EV']))
        self.data_sets_colors = dict(zip(self.data_sets, ['MUSIC+FIST EV']))
        self.energies = [7, 11, 19, 27, 39, 62]
        self.x0 = [0.1, 1.0]
        self.bounds = ((0.0001, 0.99), (0.0, 5.0))
        self.max_physical_spread = 2
        self.average_spread_range = (0.5, 1)
        self.plot = ['group_2d']  # 'indiv_2d', 'group_2d', 'vs_spread_indiv', 'vs_amp_indiv'
        self.get_2d_mins = False  # Get 2 dimensional minima in parameter space if True
        self.exclude_spreads = []  # Spread values to exclude
        self.threads = 13


def chi2_sd_slices_bs():
    pars = Pars()  # Define parameters for analysis
    datasets_data = DatasetsData(pars)  # Initialize data containers for combined data sets

    for data_set in pars.data_sets:
        chi_df_name = f'{data_set}_chi2_sum_dist_bs.csv'
        df_all = pd.read_csv(pars.base_path + chi_df_name)
        # for spread in pd.unique(df_all['spread']):
        #     print(f'spread {spread}, amps: {sorted(list(pd.unique(df_all[df_all["spread"] == spread]["amp"])))}')
        df_all = df_all[~df_all['spread'].isin(pars.exclude_spreads)]

        energies_data = EnergiesData()  # Initiate data containers for combined energies
        for energy in pars.energies:
            print(f'\nStarting {data_set} {energy}GeV')
            df = df_all[df_all['energy'] == energy]
            if df.size > 0:
                energies_data.energies.append(energy)
                datasets_data.energies[data_set].append(energy)
            else:
                print(f'No data for {data_set} {energy}GeV. Skipping')
                continue
            df = df.sort_values(by=['spread', 'amp'])

            x, y, X, Y, Z, f, f_jac, f_jac_bss, bs_cols = get_2d_interp(df)  # Get 2d interpolations
            bas_res = get_2d_mins(x, y, f_jac, pars)  # Get 2d minima from interpolation
            energies_data.min_amps.append(bas_res.x[0])
            energies_data.min_spreads.append(bas_res.x[1])

            if 'group_2d' in pars.plot or 'indiv_2d' in pars.plot:
                plot_2ds(energy, data_set, x, y, X, Y, Z, f, bas_res, pars.plot, energies_data)

            energies_data, datasets_data = spread_1d_interps(df, data_set, energy, pars, f, f_jac, f_jac_bss, bs_cols,
                                                             bas_res, energies_data, datasets_data)

        plot_energy_mins(data_set, pars, energies_data)
        datasets_data = plot_energy_mins_fit(data_set, energies_data, datasets_data)
        plot_mins_over_avg(pars, energies_data)
        datasets_data = plot_avg_mins_over_avg(data_set, pars, energies_data, datasets_data)
        plot_mins_vs_energy(data_set, pars, energies_data)

        if 'group_2d' in pars.plot:
            plot_2d_group(energies_data)

    plot_avg_norm_mins(pars, datasets_data)
    plot_avg_mins(pars, datasets_data)
    plot_fit_base(pars, datasets_data)

    plot_comp_mins(pars, datasets_data)

    plt.show()


class DatasetsData:
    def __init__(self, pars):

        self.avg_norm_min, self.avg_norm_min_sd, self.avg_min, self.avg_min_sd, self.baselines, self.baselines_sd = \
            [{data_set: [] for data_set in pars.data_sets} for i in range(6)]
        self.energies = {data_set: [] for data_set in pars.data_sets}
        self.amp_min_data = {data_set: {energy: {} for energy in pars.energies} for data_set in pars.data_sets}


class EnergiesData:
    def __init__(self):
        self.fig1_data, self.fig2_data = {}, {}
        self.energies = []
        self.min_amps, self.min_spreads = [], []
        self.energy_spreads, self.amp_mins_per_spread, self.amp_min_sds_per_spread, self.amp_min_bss_per_spread = \
            [], [], [], []
        self.avgs_per_spread = []
        self.z_min, self.z_max = 100, 0


class MyInterp2D:
    def __init__(self, x, y, z, kx=3, ky=3, jac=False):
        if len(x) == len(y) == len(z):
            df = pd.DataFrame({'x': x, 'y': y, 'z': z})
            df = df.sort_values(by=['x', 'y'])
            x, y = np.sort(pd.unique(df['x'])), np.sort(pd.unique(df['y']))
            z = np.array(df['z']).reshape(x.size, y.size)
        self.rbs_interp = interp.RectBivariateSpline(x, y, z, kx=kx, ky=ky)
        self.min = np.array([min(x), min(y), min(z.ravel())])
        self.max = np.array([max(x), max(y), max(z.ravel())])
        self.jac = jac

    def __call__(self, *coords):
        z = None
        if len(coords) == 1:
            coords = coords[0]
        if len(coords) == 2:
            # if any(coords[i] > self.max[i] for i in range(2)) or any(coords[i] < self.min[i] for i in range(2)):
            #     cent = (self.max + self.min)[:2] / 2  # Get center
            #     r_norm = np.sqrt(sum((np.array(coords) - cent)**2)) / np.sqrt(min((self.max[:2] - cent)**2))
            #     # print(f'coords: {coords}\nmin: {self.min}\nmax: {self.max}\ncent: {cent} r_norm: {r_norm}')
            #     return self.max[2] * 2 * r_norm
            z = self.rbs_interp.ev(*coords)
            if self.jac:
                dx = self.rbs_interp.ev(*coords, dx=1, dy=0)
                dy = self.rbs_interp.ev(*coords, dx=0, dy=1)
                z = (z, np.array([dx, dy]))
        else:
            print("Need 2D coords")
            return None

        return z


class BasinBounds:
    def __init__(self, xmax=[1.1, 1.1], xmin=[-1.1, -1.1]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


def basin_min(f_jac_bss, bs_cols, df_s, min_kwargs, bas_bnds, spread, niter=30, stepsize=0.4, T=0.2):
    bs_mins = []
    for f_jac_bs, bs_col in zip(f_jac_bss, bs_cols):
        def f_jac_slice_bs(amp_i):
            r = f_jac_bs(amp_i, spread)
            return r[0], r[1][0]

        bas_min_amp_bs = basinhopping(f_jac_slice_bs, df_s['amp'][df_s[bs_col].idxmin()], stepsize=stepsize, T=T,
                                      minimizer_kwargs=min_kwargs, accept_test=bas_bnds, niter=niter)
        bs_mins.append(bas_min_amp_bs.x[0])
    return spread, bs_mins


def get_edges(centers):
    centers = np.array(centers)
    edges = (centers[1:] + centers[:-1]) / 2
    edges = np.insert(edges, edges.size, centers[-1] + (centers[-1] - edges[-1]))
    edges = np.insert(edges, 0, centers[0] - (edges[0] - centers[0]))

    return edges


def trim(xs, x_low, x_high, *ys):
    return (np.array(array) for array in zip(*[(xi, *yi) for xi, *yi in zip(xs, *ys) if x_low < xi < x_high]))


def amp_spread_param(x, a, b, c):
    return c * (b * x ** 2 + a) / (1 - np.exp(-x))


def amp_spread_param2(x, a, b, c, d, e):
    return a / (1 - np.exp(-b * x)) + c * np.exp(d * x) + e


def amp_spread_param3(x, a, b, c, d, e, f):
    return a * np.exp(-b * x) / (1 - np.exp(-c * x)) + d + e * np.exp(-f * x)


def get_2d_interp(df):
    x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], df['chi2_avg']]]

    x_unq, y_unq = np.unique(x), np.unique(y)
    X, Y = np.meshgrid(get_edges(x_unq), get_edges(y_unq))
    Z = []
    for spread in np.unique(df['spread']):
        Z.append(list(np.log10(df[df['spread'] == spread]['chi2_avg'])))

    print(' Starting 2D interpolations...')
    f = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_avg']))
    f_jac = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_avg']), jac=True)
    bs_cols = [col for col in df if 'chi2_avg_bs' in col]
    f_jac_bss = []
    for bs_col in bs_cols:
        f_jac_bss.append(MyInterp2D(df['amp'], df['spread'], np.log10(df[bs_col]), jac=True))
    print(' 2D interpolations finished')
    return x, y, X, Y, Z, f, f_jac, f_jac_bss, bs_cols


def get_2d_mins(x, y, f_jac, pars):
    print(' Starting 2D minimizations...')
    # min_res = minimize(lambda x_coords: f(*x_coords), pars.x0, bounds=pars.bounds)
    bas_bnds = BasinBounds(xmax=[max(x), max(y)], xmin=[min(x), min(y)])
    min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': pars.bounds}
    bas_res = basinhopping(lambda x_coords: f_jac(*x_coords), pars.x0, minimizer_kwargs=min_kwargs,
                           accept_test=bas_bnds, niter=2000, stepsize=1, T=0.2)
    print(' 2D minimizations finished')

    return bas_res


def plot_2ds(energy, data_set, x, y, X, Y, Z, f, bas_res, plot, energies_data):
    x_intp = np.linspace(min(x), max(x), 200)
    y_intp = np.linspace(min(y), max(y), 200)
    X_intp, Y_intp = np.meshgrid(get_edges(x_intp), get_edges(y_intp))
    Z_intp = []
    for yi in y_intp:
        Zi_intp = []
        for xi in x_intp:
            Zi_intp.append(f(xi, yi))
        Z_intp.append(Zi_intp)

    if 'group_2d' in plot:
        energies_data.fig1_data.update({energy: {'X': X, 'Y': Y, 'Z': Z, 'x': x, 'y': y}})
        energies_data.fig2_data.update({energy: {'X': X_intp, 'Y': Y_intp, 'Z': Z_intp, 'basin_min': bas_res.x}})
        for z in [Z, Z_intp]:
            energies_data.z_min = np.min(z) if np.min(z) < energies_data.z_min else energies_data.z_min
            energies_data.z_max = np.max(z) if np.max(z) > energies_data.z_max else energies_data.z_max

    if 'indiv_2d' in plot:
        fig2 = plt.figure()
        fig2.canvas.manager.set_window_title(f'{data_set} {energy} GeV Interpolate 2D')
        pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
        # plt.scatter(*x0, color='white', marker='s', s=60)
        # plt.scatter(*x0, color='black', marker='s', s=25, label='Initial')
        # plt.scatter(*min_res.x, color='white', marker='*', s=60)
        # plt.scatter(*min_res.x, color='black', marker='*', s=25, label='Local Min')
        # plt.scatter(*bas_res.x, color='white', marker='^', s=60)
        # plt.scatter(*bas_res.x, color='black', marker='^', s=25, label='Basin Min')

        # plt.axhline(bounds[1][0], ls='--', color='black')
        # plt.axhline(bounds[1][1], ls='--', color='black')
        # plt.axvline(bounds[0][0], ls='--', color='black')
        # plt.axvline(bounds[0][1], ls='--', color='black')

        plt.xlabel('amp')
        plt.ylabel('spread')
        # plt.legend()
        fig2.colorbar(pcm2)
        fig2.tight_layout()
        plt.subplots_adjust(left=0.084, right=1, bottom=0.096, top=0.986)

        levels = [-2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6]
        fig3 = plt.figure()
        fig3.canvas.manager.set_window_title(f'{data_set} {energy} GeV Interpolate 2D Contours')
        plt.contour(x_intp, y_intp, Z_intp, levels)
        plt.xlabel('amp')
        plt.ylabel('spread')
        fig3.tight_layout()


def spread_1d_interps(df, data_set, energy, pars, f, f_jac, f_jac_bss, bs_cols, bas_res, energies_data, datasets_data):
    spreads = np.unique(df['spread'])
    if 'vs_amp_indiv' in pars.plot:
        color = set_vs_ampt_indiv_plot(data_set, energy, spreads, bas_res)

    bas_bnds = BasinBounds(xmax=pars.bounds[0][-1], xmin=pars.bounds[0][0])
    min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': (pars.bounds[0],)}

    energies_data.amp_mins_per_spread.append([])
    energies_data.energy_spreads.append(spreads)
    bs_jobs = []
    for spread in spreads:
        df_s = df[df['spread'] == spread]

        def f_jac_slice(amp_i):  # Use slices of the 2D interpolation
            r = f_jac(amp_i, spread)
            return r[0], r[1][0]

        bas_min_amp = basinhopping(f_jac_slice, df_s['amp'][df_s['chi2_avg'].idxmin()], stepsize=0.4, T=0.2,
                                   minimizer_kwargs=min_kwargs, accept_test=bas_bnds, niter=30)
        if 'vs_amp_indiv' in pars.plot:
            plot_vs_ampt_indiv(spread, df_s, f, bas_min_amp, color)
        energies_data.amp_mins_per_spread[-1].append(bas_min_amp.x[0])

        bs_jobs.append((f_jac_bss, bs_cols, df_s, min_kwargs, bas_bnds, spread))

    if 'vs_amp_indiv' in pars.plot:
        plt.legend(bbox_to_anchor=(1.05, 1))
        plt.tight_layout()

    print(' Getting bootstrap minima')
    bs_spreads, bs_mins = [], []
    with Pool(pars.threads) as pool:
        for bs_spread, bs_min in tqdm.tqdm(pool.istarmap(basin_min, bs_jobs), total=len(bs_jobs)):
            bs_spreads.append(bs_spread)
            bs_mins.append(bs_min)
    bs_spreads, bs_mins = zip(*sorted(zip(bs_spreads, bs_mins)))
    energies_data.amp_min_bss_per_spread.append(bs_mins)
    energies_data.amp_min_sds_per_spread.append(np.nanstd(bs_mins, axis=1))
    datasets_data.amp_min_data[data_set][energy] = {'spreads': spreads, 'mins': energies_data.amp_mins_per_spread[-1],
                                                    'mins_sd': energies_data.amp_min_sds_per_spread[-1]}
    if 'vs_spread_indiv' in pars.plot:
        plot_vs_spread_indiv(data_set, energy, df, f, bas_res)

    return energies_data, datasets_data


def set_vs_ampt_indiv_plot(data_set, energy, spreads, bas_res):
    fig3 = plt.figure()
    fig3.canvas.manager.set_window_title(f'{data_set} {energy} GeV Vs Amp')
    color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
    plt.axhline(bas_res.fun, ls='--', color='gray')
    plt.xlabel('amp')
    plt.ylabel('chi2_avg')

    return color


def plot_vs_ampt_indiv(spread, df_s, f, bas_min_amp, color):
    c = next(color)
    plt.scatter(df_s['amp'], np.log10(df_s['chi2_avg']), color=c, label=f'spread {spread}')
    zs = []
    amps = np.linspace(df_s['amp'].min(), df_s['amp'].max(), 10000)
    for amp in amps:
        zs.append(f(amp, spread))
    plt.plot(amps, zs, color=c)

    f_1d_spread = interp.interp1d(df_s['amp'], np.log10(df_s['chi2_avg']), kind='cubic')

    def f_slice(amp_i):
        return f(amp_i, spread)

    plt.plot(amps, f_1d_spread(amps), color=c, ls='--', alpha=0.6)
    plt.scatter(bas_min_amp.x[0], f_slice(bas_min_amp.x[0]), color=c, marker='P')


def plot_vs_spread_indiv(data_set, energy, df, f, bas_res):
    fig4 = plt.figure()
    fig4.canvas.manager.set_window_title(f'{data_set} {energy} GeV Vs Spread')
    amps = np.unique(df['amp'])
    color = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
    for amp in amps:
        c = next(color)
        df_a = df[df['amp'] == amp]
        plt.scatter(df_a['spread'], np.log10(df_a['chi2_avg']), color=c, label=f'amp {amp}')
        zs = []
        spreads = np.linspace(min(df_a['spread']), max(df_a['spread']), 10000)
        for spread in spreads:
            zs.append(f(amp, spread))
        plt.plot(spreads, zs, color=c)
        f_1d_amp = interp.interp1d(df_a['spread'], np.log10(df_a['chi2_avg']), kind='cubic')
        plt.plot(spreads, f_1d_amp(spreads), color=c, ls='--', alpha=0.6)
    plt.axhline(bas_res.fun, ls='--', color='gray')
    plt.xlabel('spread')
    plt.ylabel('chi2_avg')
    plt.legend(bbox_to_anchor=(1.05, 1))
    fig4.tight_layout()


def plot_energy_mins(data_set, pars, energies_data):
    fig_mins = plt.figure()
    fig_mins.canvas.manager.set_window_title(f'Mins in Amp-Spread Space')
    plt.grid()
    plt.xlim((0, 0.05))
    plt.ylim((0, 2.5))
    plt.xlabel('amp')
    plt.ylabel('spread')
    plt.axhline(pars.max_physical_spread, ls='--', color='gray')
    plt.title(data_set)
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies_data.energies))))
    for x, y, am, ame, es, e in \
            zip(energies_data.min_amps, energies_data.min_spreads, energies_data.amp_mins_per_spread,
                energies_data.amp_min_sds_per_spread, energies_data.energy_spreads, energies_data.energies):
        c = next(color)
        am, ame = np.array(am), np.array(ame)
        plt.scatter(x, y, marker='+', color=c)
        plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center', color=c)
        plt.plot(am, es, color=c, alpha=0.7, label=f'{e}GeV')
        plt.fill_betweenx(es, am - ame, am + ame, alpha=0.2, color=c)

    for spread_index in range(len(energies_data.energy_spreads[0])):
        e_sum = 0
        for energy_mins in energies_data.amp_mins_per_spread:
            e_sum += energy_mins[spread_index]
        energies_data.avgs_per_spread.append(e_sum / len(energies_data.energies))

    plt.plot(energies_data.avgs_per_spread, energies_data.energy_spreads[0], ls='--', color='black', label='Average')

    fig_mins.legend()
    fig_mins.tight_layout()

    return energies_data


def plot_energy_mins_fit(data_set, energies_data, datasets_data):
    fig_mins_fit = plt.figure()
    fig_mins_fit.canvas.manager.set_window_title(f'Mins in Amp-Spread Space with Parameterization')
    plt.grid()
    # plt.xlim((0, 0.05))
    plt.ylim((0, 1))
    plt.ylabel('amp')
    plt.xlabel('spread')
    # p0 = (0.0001, 0.005, 10)
    p0 = (0.05, 6, 50, 0.01, 0.0003, -2)
    max_spread_fit = 3.1
    plt.axvline(max_spread_fit, ls='--', color='gray')
    plt.title(data_set)
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies_data.energies))))
    # Below junky, assuming all energys have same spreads
    xs = np.linspace(min(energies_data.energy_spreads[0]), min(max(energies_data.energy_spreads[0]), max_spread_fit),
                     1000)
    plt.plot(xs, amp_spread_param3(xs, *p0), ls=':', color='gray', label='Initial Fit Guess')
    for y, x, am, ame, es, e in \
            zip(energies_data.min_amps, energies_data.min_spreads, energies_data.amp_mins_per_spread,
                energies_data.amp_min_sds_per_spread, energies_data.energy_spreads, energies_data.energies):
        c = next(color)
        # am, ame = np.array(am), np.array(ame)
        plt.scatter(x, y, marker='+', color=c)
        plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center', color=c)
        plt.plot(es, am, color=c, alpha=0.7, label=f'{e}GeV')
        # print(', '.join([str(x) for x in es]))
        # print(', '.join([str(x) for x in am]))
        # print(', '.join([str(x) for x in ame]))
        es, am, ame = trim(es, -np.inf, max_spread_fit, am, ame)
        # es, am, ame = list(zip(*[(esi, ami, amei) for esi, ami, amei in zip(es, am, ame)
        #                          if esi < max_spread_fit]))
        # es, am, ame = np.array(es), np.array(am), np.array(ame)
        try:
            popt, pcov = curve_fit(amp_spread_param3, es, am, p0=p0, sigma=ame, absolute_sigma=True)
            # print(popt)
            datasets_data.baselines[data_set].append(popt[3])
            datasets_data.baselines_sd[data_set].append(np.sqrt(np.diag(pcov))[3])
            plt.plot(xs, amp_spread_param3(xs, *popt), color=c, ls='--')
            plt.fill_between(es, am - ame, am + ame, alpha=0.2, color=c)
        except RuntimeError:
            datasets_data.baselines[data_set].append(0)
            datasets_data.baselines_sd[data_set].append(0)
            print('Parameterized fit failed')

    fig_mins_fit.legend()
    fig_mins_fit.tight_layout()

    return datasets_data


def plot_mins_over_avg(pars, energies_data):
    fig_mins_over_avg = plt.figure()
    fig_mins_over_avg.canvas.manager.set_window_title('Amp Mins divided by Average')
    plt.grid()
    plt.xlabel('spread')
    plt.ylabel('minimum amplitude / energy averaged minimum amplitude')
    plt.axvline(pars.max_physical_spread, ls='--', color='gray')
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies_data.energies))))
    for energy_mins, energy_mins_sd, energy, spreads in \
            zip(energies_data.amp_mins_per_spread, energies_data.amp_min_sds_per_spread, energies_data.energies,
                energies_data.energy_spreads):
        c = next(color)
        plt.plot(spreads, np.array(energy_mins) / np.array(energies_data.avgs_per_spread), color=c,
                 label=f'{energy}GeV')
        plt.fill_between(spreads, (np.array(energy_mins) - energy_mins_sd) / np.array(energies_data.avgs_per_spread),
                         (np.array(energy_mins) + energy_mins_sd) / np.array(energies_data.avgs_per_spread), color=c,
                         alpha=0.2)
    plt.axhline(1, color='black', ls='--')
    fig_mins_over_avg.legend()
    fig_mins_over_avg.tight_layout()


def plot_avg_mins_over_avg(data_set, pars, energies_data, datasets_data):
    fig_avg_mins_over_avg = plt.figure()
    fig_avg_mins_over_avg.canvas.manager.set_window_title('Average Amp Mins divided by Average')
    plt.grid()
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Average minimum amplitude / energy averaged minimum amplitude')
    plt.title(f'{data_set.upper()} Average Normalized Minimum Amplitude')
    for energy_mins, energy_mins_bss, spreads in \
            zip(energies_data.amp_mins_per_spread, energies_data.amp_min_bss_per_spread, energies_data.energy_spreads):
        energy_mins_avg = np.where((pars.average_spread_range[0] < np.array(spreads)) &
                                   (np.array(spreads) < pars.average_spread_range[1]),
                                   energy_mins, float('nan'))
        energy_mins_avg_bss = np.where((pars.average_spread_range[0] < np.array(spreads)) &
                                       (np.array(spreads) < pars.average_spread_range[1]),
                                       np.array(energy_mins_bss).T, float('nan'))
        # energy_mins_avg_bss = np.where((np.array(spreads) < average_spread_range[1]),
        #                                np.array(energy_mins_avg_bss).T, float('nan'))
        energy_mins = np.where(np.array(spreads) < pars.max_physical_spread, energy_mins, float('nan'))
        energy_mins_bss = np.where(np.array(spreads) < pars.max_physical_spread, np.array(energy_mins_bss).T,
                                   float('nan'))
        datasets_data.avg_norm_min[data_set].append(np.nanmean(energy_mins / np.array(energies_data.avgs_per_spread)))
        # avg_norm_min_sd[data_set].append(np.sqrt(np.nanmean(energy_mins_sd ** 2)))
        datasets_data.avg_norm_min_sd[data_set].append(np.nanstd(np.nanmean(energy_mins_bss /
                                                                            np.array(energies_data.avgs_per_spread),
                                                                            axis=1)))
        # avg_min[data_set].append(np.nanmean(energy_mins_avg))
        datasets_data.avg_min_sd[data_set].append(np.nanstd(np.nanmean(energy_mins_avg_bss, axis=1)))

        energy_mins_avg = np.ma.MaskedArray(energy_mins_avg, mask=np.isnan(energy_mins_avg))
        weights = 1 / np.nanstd(energy_mins_bss, axis=0) ** 2
        weights = np.ma.MaskedArray(weights, mask=np.isnan(weights))
        datasets_data.avg_min[data_set].append(np.ma.average(energy_mins_avg, weights=weights))
    plt.axhline(1, color='black', ls='--')
    # print(energies_data.energies)
    # print(datasets_data.avg_norm_min[data_set])
    # print(datasets_data.avg_norm_min_sd[data_set])
    plt.errorbar(energies_data.energies, datasets_data.avg_norm_min[data_set], yerr=datasets_data.avg_norm_min_sd[data_set],
                 marker='o', ls='none')
    fig_avg_mins_over_avg.tight_layout()

    return datasets_data


def plot_mins_vs_energy(data_set, pars, energies_data):
    fig_mins_vs_e = plt.figure()
    fig_mins_vs_e.canvas.manager.set_window_title(f'Amp Mins vs Energy')
    plt.grid()
    plt.xlim((0, 80))
    plt.ylim((0, 0.025))
    plt.xlabel('energy')
    plt.ylabel('amp_min')
    plt.title(data_set)
    spreads = energies_data.energy_spreads[0]  # Assume all energies have all spreads
    max_spread = 2.5
    color = iter(cm.rainbow(np.linspace(0, 1, len([x for x in spreads if x <= max_spread]))))
    inv_amp_mins = np.array(energies_data.amp_mins_per_spread).T
    for spread, min_amps in zip(spreads, inv_amp_mins):
        if spread > max_spread:
            continue
        c = next(color)
        plt.plot(energies_data.energies, min_amps, marker='o', color=c, alpha=0.7, label=f'spread {spread:.2f}')
    fig_mins_vs_e.legend()
    fig_mins_vs_e.tight_layout()


def plot_2d_group(energies_data):
    fig1_all, axes1 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig1_all.canvas.manager.set_window_title(f'Raw 2D Chi2')
    axes1 = np.ravel(axes1)
    for i, energy in enumerate(energies_data.fig1_data):
        cbar = axes1[i].pcolormesh(energies_data.fig1_data[energy]['X'], energies_data.fig1_data[energy]['Y'],
                                   energies_data.fig1_data[energy]['Z'], cmap='jet', vmin=energies_data.z_min,
                                   vmax=energies_data.z_max)
        axes1[i].scatter(energies_data.fig1_data[energy]['x'], energies_data.fig1_data[energy]['y'], color='black',
                         marker='o', s=(72. / fig1_all.dpi) ** 2)
        axes1[i].text(0.3, 1, f'{energy}GeV', color='white', fontsize='x-large')
    axes1[3].set_xlabel('amp')
    axes1[4].set_xlabel('amp')
    axes1[5].set_xlabel('amp')
    axes1[0].set_ylabel('spread')
    axes1[3].set_ylabel('spread')
    fig1_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)

    fig2_all, axes2 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig2_all.canvas.manager.set_window_title(f'Interpolated 2D Chi2')
    axes2 = np.ravel(axes2)
    for i, energy in enumerate(energies_data.fig2_data):
        cbar = axes2[i].pcolormesh(energies_data.fig2_data[energy]['X'], energies_data.fig2_data[energy]['Y'],
                                   energies_data.fig2_data[energy]['Z'], cmap='jet', vmin=energies_data.z_min,
                                   vmax=energies_data.z_max)
        axes2[i].scatter(*energies_data.fig2_data[energy]['basin_min'], color='white', marker='^', s=60)
        axes2[i].scatter(*energies_data.fig2_data[energy]['basin_min'], color='black', marker='^', s=25,
                         label='Global Minimum')
        axes2[i].text(0.3, 1, f'{energy}GeV', color='white', fontsize='x-large')
        if i + 1 == len(energies_data.fig2_data):
            axes2[i].legend()
    axes2[3].set_xlabel('amp')
    axes2[4].set_xlabel('amp')
    axes2[5].set_xlabel('amp')
    axes2[0].set_ylabel('spread')
    axes2[3].set_ylabel('spread')
    fig2_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)


def plot_avg_norm_mins(pars, datasets_data):
    fig_avg_norm_mins = plt.figure()
    fig_avg_norm_mins.canvas.manager.set_window_title('Average Norm Mins per Data Set')
    plt.grid()
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Average minimum amplitude / energy averaged minimum amplitude')
    plt.title(f'Average Normalized Minimum Amplitude')
    plt.axhline(1, color='black', ls='--')
    for data_set in pars.data_sets:
        # print(f'pars.energies: {pars.energies}')
        # print(f'datasets_data.avg_norm_min[data_set]: {datasets_data.avg_norm_min[data_set]}')
        plt.errorbar(datasets_data.energies[data_set], datasets_data.avg_norm_min[data_set],
                     yerr=datasets_data.avg_norm_min_sd[data_set], marker='o', ls='none', label=data_set.upper())
    plt.legend()
    fig_avg_norm_mins.tight_layout()


def plot_avg_mins(pars, datasets_data):
    fig_avg_mins = plt.figure()
    fig_avg_mins.canvas.manager.set_window_title('Average Mins per Data Set')
    plt.grid()
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Average minimum amplitude')
    plt.title(f'Average Minimum Amplitude')
    plt.axhline(0, color='black', ls='--')
    for data_set in pars.data_sets:
        plt.errorbar(datasets_data.energies[data_set], datasets_data.avg_min[data_set],
                     yerr=datasets_data.avg_min_sd[data_set], marker='o', ls='none',
                     label=pars.data_sets_labels[data_set], color=pars.data_sets_colors[data_set])
    plt.legend()
    fig_avg_mins.tight_layout()


def plot_fit_base(pars, datasets_data):
    fig_fit_base = plt.figure()
    fig_fit_base.canvas.manager.set_window_title('Baseline per Data Set')
    plt.grid()
    plt.xlabel('Energy (GeV)')
    plt.ylabel('Fit Amp Baseline')
    plt.title(f'Amp Baselines from Fits')
    plt.axhline(0, color='black', ls='--')
    for data_set in pars.data_sets:
        plt.errorbar(datasets_data.energies[data_set], datasets_data.baselines[data_set],
                     yerr=datasets_data.baselines_sd[data_set], marker='o', ls='none', label=data_set.upper())
    plt.legend()
    fig_fit_base.tight_layout()


def plot_comp_mins(pars, datasets_data):
    # amp_min_data = pd.DataFrame(amp_min_data)
    for energy in pars.energies:
        fig, ax = plt.subplots()
        for val in pars.average_spread_range:
            ax.axhline(val, ls='--', color='gray')
        # ax.axhline(2, ls='--', color='gray')
        for data_set in pars.data_sets:
            if energy not in datasets_data.energies[data_set]:
                continue
            # df = amp_min_data[(amp_min_data['energy'] == energy) & (amp_min_data['data_set'] == data_set)]
            # df = df.sort_values(by=['spreads'])
            # mins, spreads, mins_sd = list(df['mins']), list(df['min_sd']), list(df['spreads'])
            mins = datasets_data.amp_min_data[data_set][energy]['mins']
            mins_sd = datasets_data.amp_min_data[data_set][energy]['mins_sd']
            spreads = datasets_data.amp_min_data[data_set][energy]['spreads']
            # if len(mins) == 1 and len(spreads) == 1 and len(mins_sd) == 1:
            #     mins, spreads, mins_sd = mins[0], spreads[0], mins_sd[0]
            # else:
            #     print(f'{data_set}, {energy}GeV bad')
            #     continue
            # print(df['mins'][0])
            # print(type(df['mins'][0]))
            # print(list(df['mins'][0]))
            # print(df['spreads'][0])
            # print(type(df['spreads'][0]))
            # print(list(df['spreads'][0]))
            # print(mins, spreads, mins_sd)
            # ax.plot(df['mins'][0], df['spreads'][0], label=data_set)
            # ax.fill_betweenx(df['spreads'][0], df['mins'][0] - df['min_sd'][0], df['mins'][0] + df['min_sd'][0],
            #                  alpha=0.5)
            ax.plot(mins, spreads, label=data_set)
            ax.fill_betweenx(spreads, np.array(mins) - np.array(mins_sd), np.array(mins) + np.array(mins_sd), alpha=0.5)
        ax.set_ylim((0, 2.5))
        ax.set_xlim((0, 0.05))
        ax.grid()
        ax.set_title(f'{energy}GeV')
        ax.set_xlabel('amp')
        ax.set_ylabel('spread')
        ax.legend()
        fig.tight_layout()
        fig.canvas.manager.set_window_title(f'{energy}GeV data set comparison mins')


if __name__ == '__main__':
    main()
