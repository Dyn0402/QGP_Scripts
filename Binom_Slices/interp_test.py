#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 07 9:40 PM 2022
Created in PyCharm
Created as QGP_Scripts/interp_test

@author: Dylan Neff, Dyn04
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.interpolate as interp
from scipy.optimize import minimize
from scipy.optimize import basinhopping
import pandas as pd

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm


def main():
    # scipy_doc_test_edit()
    # chi2_data_test()
    # chi2_data_test2()
    # chi2_sd_slices()
    chi2_sd_slices_bs()

    print('donzo')


def get_edges(centers):
    centers = np.array(centers)
    edges = (centers[1:] + centers[:-1]) / 2
    edges = np.insert(edges, edges.size, centers[-1] + (centers[-1] - edges[-1]))
    edges = np.insert(edges, 0, centers[0] - (edges[0] - centers[0]))

    return edges


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


def chi2_data_test():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    'F:/Research/Results/Azimuth_Analysis/chi2_sum_dist_test.csv'
    # chi_df_name = 'chi_df_ampt_skews_cent8.csv'
    energies = [7, 11, 19, 27, 39, 62]
    x0 = [0.1, 1.0]
    bounds = ((0, 1.0), (0.0, 5.0))
    stats = {'Standard Deviation': 'sd', 'Skewness': 'skew', 'Non-excess Kurtosis': 'nek'}
    # chi_df_name = 'chi2_dist_test.csv'
    fig3_data = {}
    for stat, stat_file_code in stats.items():
        chi_df_name = f'chi_df_ampt_{stat_file_code}s_cent8.csv'
        df_all = pd.read_csv(base_path + chi_df_name)
        for spread in pd.unique(df_all['spread']):
            print(f'spread {spread}, amps: {sorted(list(pd.unique(df_all[df_all["spread"] == spread]["amp"])))}')
        exclude_spreads = []

        df_all = df_all[~df_all['spread'].isin(exclude_spreads)]
        # df_all = df_all[(df_all['spread'] != 0) & (df_all['spread'] != 4)]
        # df = df.sort_values(by=['spread', 'amp'])

        # a = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
        # print(a(1, 2))
        # print(a([1, 3]))
        # print(a.max)
        # print(a(1, 3))
        # return

        fig1_data = {}
        fig2_data = {}
        z_min, z_max = 100, 0
        min_amps = []
        min_spreads = []
        for energy in energies:
            df = df_all[df_all['energy'] == energy]
            df = df.sort_values(by=['spread', 'amp'])

            x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], np.log10(df['chi2_sum'])]]

            # fig1 = plt.figure()
            # fig1.canvas.manager.set_window_title(f'{energy} GeV Raw 2D')
            x_unq, y_unq = np.unique(x), np.unique(y)
            X, Y = np.meshgrid(get_edges(x_unq), get_edges(y_unq))
            Z = []
            for spread in np.unique(df['spread']):
                Z.append(list(np.log10(df[df['spread'] == spread]['chi2_sum'])))
            # pcm1 = plt.pcolormesh(X, Y, Z, cmap='jet')
            # plt.scatter(x, y, color='black', marker='o', s=(72./fig1.dpi)**2)
            # # plt.scatter(x, y, color='white', s=40)
            # # plt.scatter(x, y, c=z, s=25, cmap='jet')
            # plt.xlabel('amp')
            # plt.ylabel('spread')
            # fig1.colorbar(pcm1)
            # fig1.tight_layout()

            fig1_data.update({energy: {'X': X, 'Y': Y, 'Z': Z, 'x': x, 'y': y}})

            f = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
            f_jac = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']), jac=True)

            # tcks = interp.bisplrep(x, y, z, kx=5, ky=5)
            # f = lambda x_ev, y_ev: interp.bisplev(x_ev, y_ev, tcks)

            # f = interp.interp2d(x, y, z)

            min_res = minimize(lambda x_coords: f(*x_coords), x0, bounds=bounds)
            print(min_res)
            bas_bnds = BasinBounds(xmax=[max(x), max(y)], xmin=[min(x), min(y)])
            min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': bounds}
            bas_res = basinhopping(lambda x_coords: f_jac(*x_coords), x0, minimizer_kwargs=min_kwargs,
                                   accept_test=bas_bnds)
            print(bas_res)
            min_amps.append(bas_res.x[0])
            min_spreads.append(bas_res.x[1])

            # fig2 = plt.figure()
            # fig2.canvas.manager.set_window_title(f'{energy} GeV Interpolate 2D')
            x_intp = np.linspace(min(x), max(x), 200)
            y_intp = np.linspace(min(y), max(y), 200)
            # x_intp = np.linspace(-0.2, 0.4, 200)
            # y_intp = np.linspace(-1.5, 6, 200)
            X_intp, Y_intp = np.meshgrid(get_edges(x_intp), get_edges(y_intp))
            Z_intp = []
            for yi in y_intp:
                Zi_intp = []
                for xi in x_intp:
                    Zi_intp.append(f(xi, yi))
                Z_intp.append(Zi_intp)
            # pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
            # plt.scatter(*x0, color='white', marker='s', s=60)
            # plt.scatter(*x0, color='black', marker='s', s=25, label='Initial')
            # plt.scatter(*min_res.x, color='white', marker='*', s=60)
            # plt.scatter(*min_res.x, color='black', marker='*', s=25, label='Local Min')
            # plt.scatter(*bas_res.x, color='white', marker='^', s=60)
            # plt.scatter(*bas_res.x, color='black', marker='^', s=25, label='Basin Min')
            #
            # plt.axhline(bounds[1][0], ls='--', color='black')
            # plt.axhline(bounds[1][1], ls='--', color='black')
            # plt.axvline(bounds[0][0], ls='--', color='black')
            # plt.axvline(bounds[0][1], ls='--', color='black')
            # plt.axhline(np.pi, ls=':', color='gray')
            #
            # # plt.axhline(0.2, ls=':', color='black')
            # # plt.axhline(4.5, ls=':', color='black')
            # # plt.axvline(0.0, ls=':', color='black')
            # # plt.axvline(0.2, ls=':', color='black')
            # # amps = ['025', '035', '045', '175', '225', '25']
            # # for amp in amps:
            # #     plt.axvline(float('0.' + amp), color='black', ls=':')
            # # spreads = ['225', '275', '325', '375', '5']
            # # for spread in spreads:
            # #     plt.axhline(float('0.' + spread) * 10, color='black', ls=':')
            # # plt.axhline(float('0.' + spread) * 10, color='black', ls=':', label='Simulation in Progress')
            #
            # plt.xlabel('amp')
            # plt.ylabel('spread')
            # plt.legend()
            # fig2.colorbar(pcm2)
            # fig2.tight_layout()
            # plt.subplots_adjust(left=0.084, right=1, bottom=0.096, top=0.986)

            fig2_data.update({energy: {'X': X_intp, 'Y': Y_intp, 'Z': Z_intp, 'basin_min': bas_res.x}})

            for z in [Z, Z_intp]:
                z_min = np.min(z) if np.min(z) < z_min else z_min
                z_max = np.max(z) if np.max(z) > z_max else z_max

            # fig3 = plt.figure()
            # fig3.canvas.manager.set_window_title(f'{energy} GeV Vs Amp')
            # spreads = np.unique(df['spread'])
            # color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
            # for spread in spreads:
            #     c = next(color)
            #     df_s = df[df['spread'] == spread]
            #     plt.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), color=c, label=f'spread {spread}')
            #     zs = []
            #     amps = np.linspace(min(df_s['amp']), max(df_s['amp']), 1000)
            #     for amp in amps:
            #         zs.append(f(amp, spread))
            #     plt.plot(amps, zs, color=c)
            #     f_1d_spread = interp.interp1d(df_s['amp'], np.log10(df_s['chi2_sum']), kind='cubic')
            #     plt.plot(amps, f_1d_spread(amps), color=c, ls='--', alpha=0.6)
            # plt.xlabel('amp')
            # plt.ylabel('chi2_sum')
            # plt.legend(bbox_to_anchor=(1.05, 1))
            # fig3.tight_layout()
            #
            # fig4 = plt.figure()
            # fig4.canvas.manager.set_window_title(f'{energy} GeV Vs Spread')
            # amps = np.unique(df['amp'])
            # color = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
            # for amp in amps:
            #     c = next(color)
            #     df_a = df[df['amp'] == amp]
            #     plt.scatter(df_a['spread'], np.log10(df_a['chi2_sum']), color=c, label=f'amp {amp}')
            #     zs = []
            #     spreads = np.linspace(min(df_a['spread']), max(df_a['spread']), 1000)
            #     for spread in spreads:
            #         zs.append(f(amp, spread))
            #     plt.plot(spreads, zs, color=c)
            #     f_1d_amp = interp.interp1d(df_a['spread'], np.log10(df_a['chi2_sum']), kind='cubic')
            #     plt.plot(spreads, f_1d_amp(spreads), color=c, ls='--', alpha=0.6)
            # plt.xlabel('spread')
            # plt.ylabel('chi2_sum')
            # plt.legend(bbox_to_anchor=(1.05, 1))
            # fig4.tight_layout()

        # fig_mins = plt.figure()
        # plt.grid()
        # plt.scatter(min_amps, min_spreads)
        # plt.xlim(bounds[0])
        # plt.ylim(bounds[1])
        # plt.xlabel('amp')
        # plt.ylabel('spread')
        # plt.title('AMPT')
        # for x, y, e in zip(min_amps, min_spreads, energies):
        #     plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
        # fig_mins.tight_layout()
        fig3_data.update({stat: [min_amps, min_spreads]})

        fig1_all, axes1 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
        fig1_all.canvas.manager.set_window_title(f'Raw 2D Chi2 {stat}')
        axes1 = np.ravel(axes1)
        for i, energy in enumerate(fig1_data):
            cbar = axes1[i].pcolormesh(fig1_data[energy]['X'], fig1_data[energy]['Y'], fig1_data[energy]['Z'],
                                       cmap='jet', vmin=z_min, vmax=z_max)
            axes1[i].scatter(fig1_data[energy]['x'], fig1_data[energy]['y'], color='black', marker='o',
                             s=(72. / fig1_all.dpi) ** 2)
            axes1[i].text(0.3, 1, f'{energy}GeV', color='white')
        axes1[3].set_xlabel('amp')
        axes1[4].set_xlabel('amp')
        axes1[5].set_xlabel('amp')
        axes1[0].set_ylabel('spread')
        axes1[3].set_ylabel('spread')
        # fig1_all.colorbar(cbar)
        fig1_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)
        # fig1_all.tight_layout()

        fig2_all, axes2 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
        fig2_all.canvas.manager.set_window_title(f'Interpolated 2D Chi2 {stat}')
        axes2 = np.ravel(axes2)
        for i, energy in enumerate(fig2_data):
            cbar = axes2[i].pcolormesh(fig2_data[energy]['X'], fig2_data[energy]['Y'], fig2_data[energy]['Z'],
                                       cmap='jet',
                                       vmin=z_min, vmax=z_max)
            axes2[i].scatter(*fig2_data[energy]['basin_min'], color='white', marker='^', s=60)
            axes2[i].scatter(*fig2_data[energy]['basin_min'], color='black', marker='^', s=25, label='Basin Min')
            axes2[i].text(0.3, 1, f'{energy}GeV', color='white')
        axes2[3].set_xlabel('amp')
        axes2[4].set_xlabel('amp')
        axes2[5].set_xlabel('amp')
        axes2[0].set_ylabel('spread')
        axes2[3].set_ylabel('spread')
        # fig2_all.colorbar(cbar)
        fig2_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)

    fig_mins = plt.figure()
    plt.grid()
    plt.xlim(bounds[0])
    plt.ylim(bounds[1])
    plt.xlabel('amp')
    plt.ylabel('spread')
    plt.title('AMPT')
    for stat, stat_data in fig3_data.items():
        min_amps, min_spreads = stat_data
        plt.scatter(min_amps, min_spreads, label=stat)
        for x, y, e in zip(min_amps, min_spreads, energies):
            plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
    plt.legend()
    fig_mins.tight_layout()

    plt.show()


def chi2_data_test2():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # chi_df_name = 'ampt_new_chi2_sum_dist.csv'
    chi_df_name = 'bes_chi2_sum_dist.csv'
    energies = [7, 11, 19, 27, 39, 62]
    x0 = [0.1, 1.0]
    bounds = ((0, 1.0), (0.0, 5.0))
    # chi_df_name = 'chi2_dist_test.csv'
    df_all = pd.read_csv(base_path + chi_df_name)
    for spread in pd.unique(df_all['spread']):
        print(f'spread {spread}, amps: {sorted(list(pd.unique(df_all[df_all["spread"] == spread]["amp"])))}')
    exclude_spreads = []

    df_all = df_all[~df_all['spread'].isin(exclude_spreads)]
    # df_all = df_all[(df_all['spread'] != 0) & (df_all['spread'] != 4)]
    # df = df.sort_values(by=['spread', 'amp'])

    # a = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
    # print(a(1, 2))
    # print(a([1, 3]))
    # print(a.max)
    # print(a(1, 3))
    # return

    fig1_data = {}
    fig2_data = {}
    z_min, z_max = 100, 0
    min_amps = []
    min_spreads = []
    for energy in energies:
        print(energy)
        df = df_all[df_all['energy'] == energy]
        df = df.sort_values(by=['spread', 'amp'])

        x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], np.log10(df['chi2_sum'])]]

        # fig1 = plt.figure()
        # fig1.canvas.manager.set_window_title(f'{energy} GeV Raw 2D')
        x_unq, y_unq = np.unique(x), np.unique(y)
        X, Y = np.meshgrid(get_edges(x_unq), get_edges(y_unq))
        Z = []
        for spread in np.unique(df['spread']):
            Z.append(list(np.log10(df[df['spread'] == spread]['chi2_sum'])))
        # pcm1 = plt.pcolormesh(X, Y, Z, cmap='jet')
        # plt.scatter(x, y, color='black', marker='o', s=(72./fig1.dpi)**2)
        # # plt.scatter(x, y, color='white', s=40)
        # # plt.scatter(x, y, c=z, s=25, cmap='jet')
        # plt.xlabel('amp')
        # plt.ylabel('spread')
        # fig1.colorbar(pcm1)
        # fig1.tight_layout()

        fig1_data.update({energy: {'X': X, 'Y': Y, 'Z': Z, 'x': x, 'y': y}})

        f = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
        f_jac = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']), jac=True)

        # tcks = interp.bisplrep(x, y, z, kx=5, ky=5)
        # f = lambda x_ev, y_ev: interp.bisplev(x_ev, y_ev, tcks)

        # f = interp.interp2d(x, y, z)

        min_res = minimize(lambda x_coords: f(*x_coords), x0, bounds=bounds)
        print(min_res)
        bas_bnds = BasinBounds(xmax=[max(x), max(y)], xmin=[min(x), min(y)])
        min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': bounds}
        bas_res = basinhopping(lambda x_coords: f_jac(*x_coords), x0, minimizer_kwargs=min_kwargs,
                               accept_test=bas_bnds, niter=2000, stepsize=1, T=0.2)
        print(bas_res)
        min_amps.append(bas_res.x[0])
        min_spreads.append(bas_res.x[1])

        fig2 = plt.figure()
        fig2.canvas.manager.set_window_title(f'{energy} GeV Interpolate 2D')
        x_intp = np.linspace(min(x), max(x), 200)
        y_intp = np.linspace(min(y), max(y), 200)
        # x_intp = np.linspace(-0.2, 0.4, 200)
        # y_intp = np.linspace(-1.5, 6, 200)
        X_intp, Y_intp = np.meshgrid(get_edges(x_intp), get_edges(y_intp))
        Z_intp = []
        for yi in y_intp:
            Zi_intp = []
            for xi in x_intp:
                Zi_intp.append(f(xi, yi))
            Z_intp.append(Zi_intp)
        pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
        plt.scatter(*x0, color='white', marker='s', s=60)
        plt.scatter(*x0, color='black', marker='s', s=25, label='Initial')
        plt.scatter(*min_res.x, color='white', marker='*', s=60)
        plt.scatter(*min_res.x, color='black', marker='*', s=25, label='Local Min')
        plt.scatter(*bas_res.x, color='white', marker='^', s=60)
        plt.scatter(*bas_res.x, color='black', marker='^', s=25, label='Basin Min')

        plt.axhline(bounds[1][0], ls='--', color='black')
        plt.axhline(bounds[1][1], ls='--', color='black')
        plt.axvline(bounds[0][0], ls='--', color='black')
        plt.axvline(bounds[0][1], ls='--', color='black')

        # plt.axhline(0.2, ls=':', color='black')
        # plt.axhline(4.5, ls=':', color='black')
        # plt.axvline(0.0, ls=':', color='black')
        # plt.axvline(0.2, ls=':', color='black')
        # amps = ['025', '035', '045', '175', '225', '25']
        # for amp in amps:
        #     plt.axvline(float('0.' + amp), color='black', ls=':')
        # spreads = ['225', '275', '325', '375', '5']
        # for spread in spreads:
        #     plt.axhline(float('0.' + spread) * 10, color='black', ls=':')
        # plt.axhline(float('0.' + spread) * 10, color='black', ls=':', label='Simulation in Progress')

        plt.xlabel('amp')
        plt.ylabel('spread')
        plt.legend()
        fig2.colorbar(pcm2)
        fig2.tight_layout()
        plt.subplots_adjust(left=0.084, right=1, bottom=0.096, top=0.986)

        fig2_data.update({energy: {'X': X_intp, 'Y': Y_intp, 'Z': Z_intp, 'basin_min': bas_res.x}})

        for z in [Z, Z_intp]:
            z_min = np.min(z) if np.min(z) < z_min else z_min
            z_max = np.max(z) if np.max(z) > z_max else z_max

        fig3 = plt.figure()
        fig3.canvas.manager.set_window_title(f'{energy} GeV Vs Amp')
        spreads = np.unique(df['spread'])
        color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
        for spread in spreads:
            c = next(color)
            df_s = df[df['spread'] == spread]
            plt.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), color=c, label=f'spread {spread}')
            zs = []
            amps = np.linspace(min(df_s['amp']), max(df_s['amp']), 10000)
            for amp in amps:
                zs.append(f(amp, spread))
            plt.plot(amps, zs, color=c)
            f_1d_spread = interp.interp1d(df_s['amp'], np.log10(df_s['chi2_sum']), kind='cubic')
            plt.plot(amps, f_1d_spread(amps), color=c, ls='--', alpha=0.6)
        plt.axhline(bas_res.fun, ls='--', color='gray')
        plt.xlabel('amp')
        plt.ylabel('chi2_sum')
        plt.legend(bbox_to_anchor=(1.05, 1))
        fig3.tight_layout()

        fig4 = plt.figure()
        fig4.canvas.manager.set_window_title(f'{energy} GeV Vs Spread')
        amps = np.unique(df['amp'])
        color = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
        for amp in amps:
            c = next(color)
            df_a = df[df['amp'] == amp]
            plt.scatter(df_a['spread'], np.log10(df_a['chi2_sum']), color=c, label=f'amp {amp}')
            zs = []
            spreads = np.linspace(min(df_a['spread']), max(df_a['spread']), 10000)
            for spread in spreads:
                zs.append(f(amp, spread))
            plt.plot(spreads, zs, color=c)
            f_1d_amp = interp.interp1d(df_a['spread'], np.log10(df_a['chi2_sum']), kind='cubic')
            plt.plot(spreads, f_1d_amp(spreads), color=c, ls='--', alpha=0.6)
        plt.axhline(bas_res.fun, ls='--', color='gray')
        plt.xlabel('spread')
        plt.ylabel('chi2_sum')
        plt.legend(bbox_to_anchor=(1.05, 1))
        fig4.tight_layout()

    fig_mins = plt.figure()
    plt.grid()
    plt.scatter(min_amps, min_spreads)
    plt.xlim(bounds[0])
    plt.ylim(bounds[1])
    plt.xlabel('amp')
    plt.ylabel('spread')
    plt.title('AMPT')
    for x, y, e in zip(min_amps, min_spreads, energies):
        plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
    fig_mins.tight_layout()
    # fig3_data.update({stat: [min_amps, min_spreads]})

    fig1_all, axes1 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig1_all.canvas.manager.set_window_title(f'Raw 2D Chi2')
    axes1 = np.ravel(axes1)
    for i, energy in enumerate(fig1_data):
        cbar = axes1[i].pcolormesh(fig1_data[energy]['X'], fig1_data[energy]['Y'], fig1_data[energy]['Z'],
                                   cmap='jet', vmin=z_min, vmax=z_max)
        axes1[i].scatter(fig1_data[energy]['x'], fig1_data[energy]['y'], color='black', marker='o',
                         s=(72. / fig1_all.dpi) ** 2)
        axes1[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes1[3].set_xlabel('amp')
    axes1[4].set_xlabel('amp')
    axes1[5].set_xlabel('amp')
    axes1[0].set_ylabel('spread')
    axes1[3].set_ylabel('spread')
    # fig1_all.colorbar(cbar)
    fig1_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)
    # fig1_all.tight_layout()

    fig2_all, axes2 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig2_all.canvas.manager.set_window_title(f'Interpolated 2D Chi2')
    axes2 = np.ravel(axes2)
    for i, energy in enumerate(fig2_data):
        cbar = axes2[i].pcolormesh(fig2_data[energy]['X'], fig2_data[energy]['Y'], fig2_data[energy]['Z'],
                                   cmap='jet',
                                   vmin=z_min, vmax=z_max)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='white', marker='^', s=60)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='black', marker='^', s=25, label='Basin Min')
        axes2[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes2[3].set_xlabel('amp')
    axes2[4].set_xlabel('amp')
    axes2[5].set_xlabel('amp')
    axes2[0].set_ylabel('spread')
    axes2[3].set_ylabel('spread')
    # fig2_all.colorbar(cbar)
    fig2_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)

    # fig_mins = plt.figure()
    # plt.grid()
    # plt.xlim(bounds[0])
    # plt.ylim(bounds[1])
    # plt.xlabel('amp')
    # plt.ylabel('spread')
    # plt.title('AMPT')
    # for stat, stat_data in fig3_data.items():
    #     min_amps, min_spreads = stat_data
    #     plt.scatter(min_amps, min_spreads, label=stat)
    #     for x, y, e in zip(min_amps, min_spreads, energies):
    #         plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
    # plt.legend()
    # fig_mins.tight_layout()

    plt.show()


def chi2_sd_slices():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # chi_df_name = 'ampt_new_chi2_sum_dist.csv'
    data_set = 'bes'
    chi_df_name = f'{data_set}_chi2_sum_dist.csv'
    energies = [7, 11, 19, 27, 39, 62]
    x0 = [0.1, 1.0]
    bounds = ((0.0001, 0.99), (0.0, 5.0))
    # chi_df_name = 'chi2_dist_test.csv'
    df_all = pd.read_csv(base_path + chi_df_name)
    for spread in pd.unique(df_all['spread']):
        print(f'spread {spread}, amps: {sorted(list(pd.unique(df_all[df_all["spread"] == spread]["amp"])))}')
    exclude_spreads = []

    df_all = df_all[~df_all['spread'].isin(exclude_spreads)]
    # df_all = df_all[(df_all['spread'] != 0) & (df_all['spread'] != 4)]
    # df = df.sort_values(by=['spread', 'amp'])

    # a = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
    # print(a(1, 2))
    # print(a([1, 3]))
    # print(a.max)
    # print(a(1, 3))
    # return

    fig1_data = {}
    fig2_data = {}
    z_min, z_max = 100, 0
    min_amps = []
    min_spreads = []
    energy_spreads, amp_mins_per_spread = [], []
    for energy in energies:
        print(energy)
        df = df_all[df_all['energy'] == energy]
        df = df.sort_values(by=['spread', 'amp'])

        x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], np.log10(df['chi2_sum'])]]

        # fig1 = plt.figure()
        # fig1.canvas.manager.set_window_title(f'{energy} GeV Raw 2D')
        x_unq, y_unq = np.unique(x), np.unique(y)
        X, Y = np.meshgrid(get_edges(x_unq), get_edges(y_unq))
        Z = []
        for spread in np.unique(df['spread']):
            Z.append(list(np.log10(df[df['spread'] == spread]['chi2_sum'])))
        # pcm1 = plt.pcolormesh(X, Y, Z, cmap='jet')
        # plt.scatter(x, y, color='black', marker='o', s=(72./fig1.dpi)**2)
        # # plt.scatter(x, y, color='white', s=40)
        # # plt.scatter(x, y, c=z, s=25, cmap='jet')
        # plt.xlabel('amp')
        # plt.ylabel('spread')
        # fig1.colorbar(pcm1)
        # fig1.tight_layout()

        fig1_data.update({energy: {'X': X, 'Y': Y, 'Z': Z, 'x': x, 'y': y}})

        f = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
        f_jac = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']), jac=True)

        # tcks = interp.bisplrep(x, y, z, kx=5, ky=5)
        # f = lambda x_ev, y_ev: interp.bisplev(x_ev, y_ev, tcks)

        # f = interp.interp2d(x, y, z)

        min_res = minimize(lambda x_coords: f(*x_coords), x0, bounds=bounds)
        print(min_res)
        bas_bnds = BasinBounds(xmax=[max(x), max(y)], xmin=[min(x), min(y)])
        min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': bounds}
        bas_res = basinhopping(lambda x_coords: f_jac(*x_coords), x0, minimizer_kwargs=min_kwargs,
                               accept_test=bas_bnds, niter=2000, stepsize=1, T=0.2)
        print(bas_res)
        min_amps.append(bas_res.x[0])
        min_spreads.append(bas_res.x[1])

        fig2 = plt.figure()
        fig2.canvas.manager.set_window_title(f'{energy} GeV Interpolate 2D')
        x_intp = np.linspace(min(x), max(x), 200)
        y_intp = np.linspace(min(y), max(y), 200)
        # x_intp = np.linspace(-0.2, 0.4, 200)
        # y_intp = np.linspace(-1.5, 6, 200)
        X_intp, Y_intp = np.meshgrid(get_edges(x_intp), get_edges(y_intp))
        Z_intp = []
        for yi in y_intp:
            Zi_intp = []
            for xi in x_intp:
                Zi_intp.append(f(xi, yi))
            Z_intp.append(Zi_intp)
        pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
        plt.scatter(*x0, color='white', marker='s', s=60)
        plt.scatter(*x0, color='black', marker='s', s=25, label='Initial')
        plt.scatter(*min_res.x, color='white', marker='*', s=60)
        plt.scatter(*min_res.x, color='black', marker='*', s=25, label='Local Min')
        plt.scatter(*bas_res.x, color='white', marker='^', s=60)
        plt.scatter(*bas_res.x, color='black', marker='^', s=25, label='Basin Min')

        plt.axhline(bounds[1][0], ls='--', color='black')
        plt.axhline(bounds[1][1], ls='--', color='black')
        plt.axvline(bounds[0][0], ls='--', color='black')
        plt.axvline(bounds[0][1], ls='--', color='black')

        # plt.axhline(0.2, ls=':', color='black')
        # plt.axhline(4.5, ls=':', color='black')
        # plt.axvline(0.0, ls=':', color='black')
        # plt.axvline(0.2, ls=':', color='black')
        # amps = ['025', '035', '045', '175', '225', '25']
        # for amp in amps:
        #     plt.axvline(float('0.' + amp), color='black', ls=':')
        # spreads = ['225', '275', '325', '375', '5']
        # for spread in spreads:
        #     plt.axhline(float('0.' + spread) * 10, color='black', ls=':')
        # plt.axhline(float('0.' + spread) * 10, color='black', ls=':', label='Simulation in Progress')

        plt.xlabel('amp')
        plt.ylabel('spread')
        plt.legend()
        fig2.colorbar(pcm2)
        fig2.tight_layout()
        plt.subplots_adjust(left=0.084, right=1, bottom=0.096, top=0.986)

        fig2_data.update({energy: {'X': X_intp, 'Y': Y_intp, 'Z': Z_intp, 'basin_min': bas_res.x}})

        for z in [Z, Z_intp]:
            z_min = np.min(z) if np.min(z) < z_min else z_min
            z_max = np.max(z) if np.max(z) > z_max else z_max

        fig3 = plt.figure()
        fig3.canvas.manager.set_window_title(f'{energy} GeV Vs Amp')
        spreads = np.unique(df['spread'])
        color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
        amp_mins_per_spread.append([])
        energy_spreads.append(spreads)
        for spread in spreads:
            c = next(color)
            df_s = df[df['spread'] == spread]
            plt.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), color=c, label=f'spread {spread}')
            zs = []
            amps = np.linspace(min(df_s['amp']), max(df_s['amp']), 10000)
            for amp in amps:
                zs.append(f(amp, spread))
            plt.plot(amps, zs, color=c)
            f_1d_spread = interp.interp1d(df_s['amp'], np.log10(df_s['chi2_sum']), kind='cubic')
            # print(df_s['amp'])
            # min_amp = minimize(f_1d_spread, 0.5, bounds=(bounds[0],))
            bas_bnds = BasinBounds(xmax=bounds[0][-1], xmin=bounds[0][0])
            # min_kwargs = {'method': 'L-BFGS-B', 'bounds': (bounds[0],)}
            # bas_min_amp = basinhopping(f_1d_spread, 0.5, minimizer_kwargs=min_kwargs,
            #                            accept_test=bas_bnds, niter=10, stepsize=1, T=1.2)
            min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': (bounds[0],)}

            def f_jac_slice(amp_i):
                r = f_jac(amp_i, spread)
                return r[0], r[1][0]

            bas_min_amp = basinhopping(f_jac_slice, 0.5, minimizer_kwargs=min_kwargs,
                                       accept_test=bas_bnds, niter=30, stepsize=0.4, T=0.2)
            amp_mins_per_spread[-1].append(bas_min_amp.x[0])

            # plt.plot(amps, f_1d_spread(amps), color=c, ls='--', alpha=0.6)
            # plt.scatter(bas_min_amp.x[0], f_1d_spread(bas_min_amp.x[0]), color=c, marker='P')
            def f_slice(amp_i):
                return f(amp_i, spread)

            plt.plot(amps, f_slice(amps), color=c, ls='--', alpha=0.6)
            plt.scatter(bas_min_amp.x[0], f_slice(bas_min_amp.x[0]), color=c, marker='P')
        plt.axhline(bas_res.fun, ls='--', color='gray')
        plt.xlabel('amp')
        plt.ylabel('chi2_sum')
        plt.legend(bbox_to_anchor=(1.05, 1))
        fig3.tight_layout()

        fig4 = plt.figure()
        fig4.canvas.manager.set_window_title(f'{energy} GeV Vs Spread')
        amps = np.unique(df['amp'])
        color = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
        for amp in amps:
            c = next(color)
            df_a = df[df['amp'] == amp]
            plt.scatter(df_a['spread'], np.log10(df_a['chi2_sum']), color=c, label=f'amp {amp}')
            zs = []
            spreads = np.linspace(min(df_a['spread']), max(df_a['spread']), 10000)
            for spread in spreads:
                zs.append(f(amp, spread))
            plt.plot(spreads, zs, color=c)
            f_1d_amp = interp.interp1d(df_a['spread'], np.log10(df_a['chi2_sum']), kind='cubic')
            plt.plot(spreads, f_1d_amp(spreads), color=c, ls='--', alpha=0.6)
        plt.axhline(bas_res.fun, ls='--', color='gray')
        plt.xlabel('spread')
        plt.ylabel('chi2_sum')
        plt.legend(bbox_to_anchor=(1.05, 1))
        fig4.tight_layout()

    fig_mins = plt.figure()
    fig_mins.canvas.manager.set_window_title(f'Mins in Amp-Spread Space')
    plt.grid()
    # plt.scatter(min_amps, min_spreads)
    # plt.xlim(bounds[0])
    # plt.ylim(bounds[1])
    plt.xlim((0, 0.05))
    plt.ylim((0, 2.5))
    plt.xlabel('amp')
    plt.ylabel('spread')
    plt.title(data_set)
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies))))
    for x, y, am, es, e in zip(min_amps, min_spreads, amp_mins_per_spread, energy_spreads, energies):
        c = next(color)
        plt.scatter(x, y, marker='+', color=c)
        plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center', color=c)
        plt.plot(am, es, color=c, alpha=0.7, label=f'{e}GeV')

    avgs_per_spread = []
    for spread_index in range(len(energy_spreads[0])):
        e_sum = 0
        for energy_mins in amp_mins_per_spread:
            e_sum += energy_mins[spread_index]
        avgs_per_spread.append(e_sum / len(energies))

    plt.plot(avgs_per_spread, energy_spreads[0], ls='--', color='black', label='Average')

    fig_mins.legend()
    fig_mins.tight_layout()
    # fig3_data.update({stat: [min_amps, min_spreads]})

    fig_mins_over_avg = plt.figure()
    fig_mins_over_avg.canvas.manager.set_window_title('Amp Mins divided by Average')
    plt.grid()
    plt.xlabel('spread')
    plt.ylabel('minimum amplitude / energy averaged minimum amplitude')
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies))))
    for energy_mins, energy, spreads in zip(amp_mins_per_spread, energies, energy_spreads):
        c = next(color)
        plt.plot(spreads, np.array(energy_mins) / np.array(avgs_per_spread), color=c, label=f'{energy}GeV')
    plt.axhline(1, color='black', ls='--')
    fig_mins_over_avg.legend()
    fig_mins_over_avg.tight_layout()

    fig_mins_vs_e = plt.figure()
    fig_mins_vs_e.canvas.manager.set_window_title(f'Amp Mins vs Energy')
    plt.grid()
    # plt.scatter(min_amps, min_spreads)
    plt.xlim((0, 80))
    plt.ylim((0, 0.025))
    plt.xlabel('energy')
    plt.ylabel('amp_min')
    plt.title(data_set)
    spreads = energy_spreads[0]  # Assume all energies have all spreads
    max_spread = 2.5
    color = iter(cm.rainbow(np.linspace(0, 1, len([x for x in spreads if x <= max_spread]))))
    inv_amp_mins = np.array(amp_mins_per_spread).T
    for spread, min_amps in zip(spreads, inv_amp_mins):
        if spread > max_spread:
            continue
        c = next(color)
        plt.plot(energies, min_amps, marker='o', color=c, alpha=0.7, label=f'spread {spread:.2f}')
    fig_mins_vs_e.legend()
    fig_mins_vs_e.tight_layout()

    fig1_all, axes1 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig1_all.canvas.manager.set_window_title(f'Raw 2D Chi2')
    axes1 = np.ravel(axes1)
    for i, energy in enumerate(fig1_data):
        cbar = axes1[i].pcolormesh(fig1_data[energy]['X'], fig1_data[energy]['Y'], fig1_data[energy]['Z'],
                                   cmap='jet', vmin=z_min, vmax=z_max)
        axes1[i].scatter(fig1_data[energy]['x'], fig1_data[energy]['y'], color='black', marker='o',
                         s=(72. / fig1_all.dpi) ** 2)
        axes1[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes1[3].set_xlabel('amp')
    axes1[4].set_xlabel('amp')
    axes1[5].set_xlabel('amp')
    axes1[0].set_ylabel('spread')
    axes1[3].set_ylabel('spread')
    # fig1_all.colorbar(cbar)
    fig1_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)
    # fig1_all.tight_layout()

    fig2_all, axes2 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig2_all.canvas.manager.set_window_title(f'Interpolated 2D Chi2')
    axes2 = np.ravel(axes2)
    for i, energy in enumerate(fig2_data):
        cbar = axes2[i].pcolormesh(fig2_data[energy]['X'], fig2_data[energy]['Y'], fig2_data[energy]['Z'],
                                   cmap='jet',
                                   vmin=z_min, vmax=z_max)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='white', marker='^', s=60)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='black', marker='^', s=25, label='Basin Min')
        axes2[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes2[3].set_xlabel('amp')
    axes2[4].set_xlabel('amp')
    axes2[5].set_xlabel('amp')
    axes2[0].set_ylabel('spread')
    axes2[3].set_ylabel('spread')
    # fig2_all.colorbar(cbar)
    fig2_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)

    # fig_mins = plt.figure()
    # plt.grid()
    # plt.xlim(bounds[0])
    # plt.ylim(bounds[1])
    # plt.xlabel('amp')
    # plt.ylabel('spread')
    # plt.title('AMPT')
    # for stat, stat_data in fig3_data.items():
    #     min_amps, min_spreads = stat_data
    #     plt.scatter(min_amps, min_spreads, label=stat)
    #     for x, y, e in zip(min_amps, min_spreads, energies):
    #         plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center')
    # plt.legend()
    # fig_mins.tight_layout()

    plt.show()


def chi2_sd_slices_bs():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    chi_df_name = 'chi2_sum_dist_ampt_new_bs.csv'
    data_set = 'ampt_new'
    # chi_df_name = f'{data_set}_chi2_sum_dist.csv'
    energies = [7, 11, 19, 27, 39, 62]
    x0 = [0.1, 1.0]
    bounds = ((0.0001, 0.99), (0.0, 5.0))
    threads = 15
    # chi_df_name = 'chi2_dist_test.csv'
    df_all = pd.read_csv(base_path + chi_df_name)
    for spread in pd.unique(df_all['spread']):
        print(f'spread {spread}, amps: {sorted(list(pd.unique(df_all[df_all["spread"] == spread]["amp"])))}')
    exclude_spreads = []

    df_all = df_all[~df_all['spread'].isin(exclude_spreads)]
    # df_all = df_all[(df_all['spread'] != 0) & (df_all['spread'] != 4)]
    # df = df.sort_values(by=['spread', 'amp'])

    # a = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_sum']))
    # print(a(1, 2))
    # print(a([1, 3]))
    # print(a.max)
    # print(a(1, 3))
    # return

    fig1_data = {}
    fig2_data = {}
    z_min, z_max = 100, 0
    min_amps = []
    min_spreads = []
    energy_spreads, amp_mins_per_spread, amp_min_sds_per_spread = [], [], []
    for energy in energies:
        print(f'\nStarting {energy}GeV')
        df = df_all[df_all['energy'] == energy]
        df = df.sort_values(by=['spread', 'amp'])

        x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], df['chi2_avg']]]

        # fig1 = plt.figure()
        # fig1.canvas.manager.set_window_title(f'{energy} GeV Raw 2D')
        x_unq, y_unq = np.unique(x), np.unique(y)
        X, Y = np.meshgrid(get_edges(x_unq), get_edges(y_unq))
        Z = []
        for spread in np.unique(df['spread']):
            Z.append(list(np.log10(df[df['spread'] == spread]['chi2_avg'])))
        # pcm1 = plt.pcolormesh(X, Y, Z, cmap='jet')
        # plt.scatter(x, y, color='black', marker='o', s=(72./fig1.dpi)**2)
        # # plt.scatter(x, y, color='white', s=40)
        # # plt.scatter(x, y, c=z, s=25, cmap='jet')
        # plt.xlabel('amp')
        # plt.ylabel('spread')
        # fig1.colorbar(pcm1)
        # fig1.tight_layout()

        fig1_data.update({energy: {'X': X, 'Y': Y, 'Z': Z, 'x': x, 'y': y}})

        print('Starting 2D interpolations...')
        f = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_avg']))
        f_jac = MyInterp2D(df['amp'], df['spread'], np.log10(df['chi2_avg']), jac=True)
        bs_cols = [col for col in df if 'chi2_avg_bs' in col]
        f_jac_bss = []
        for bs_col in bs_cols:
            f_jac_bss.append(MyInterp2D(df['amp'], df['spread'], np.log10(df[bs_col]), jac=True))
        print('2D interpolations finished')

        # tcks = interp.bisplrep(x, y, z, kx=5, ky=5)
        # f = lambda x_ev, y_ev: interp.bisplev(x_ev, y_ev, tcks)

        # f = interp.interp2d(x, y, z)

        min_res = minimize(lambda x_coords: f(*x_coords), x0, bounds=bounds)
        # print(min_res)
        bas_bnds = BasinBounds(xmax=[max(x), max(y)], xmin=[min(x), min(y)])
        min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': bounds}
        bas_res = basinhopping(lambda x_coords: f_jac(*x_coords), x0, minimizer_kwargs=min_kwargs,
                               accept_test=bas_bnds, niter=2000, stepsize=1, T=0.2)
        # print(bas_res)
        min_amps.append(bas_res.x[0])
        min_spreads.append(bas_res.x[1])

        fig2 = plt.figure()
        fig2.canvas.manager.set_window_title(f'{energy} GeV Interpolate 2D')
        x_intp = np.linspace(min(x), max(x), 200)
        y_intp = np.linspace(min(y), max(y), 200)
        # x_intp = np.linspace(-0.2, 0.4, 200)
        # y_intp = np.linspace(-1.5, 6, 200)
        X_intp, Y_intp = np.meshgrid(get_edges(x_intp), get_edges(y_intp))
        Z_intp = []
        for yi in y_intp:
            Zi_intp = []
            for xi in x_intp:
                Zi_intp.append(f(xi, yi))
            Z_intp.append(Zi_intp)
        pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
        plt.scatter(*x0, color='white', marker='s', s=60)
        plt.scatter(*x0, color='black', marker='s', s=25, label='Initial')
        plt.scatter(*min_res.x, color='white', marker='*', s=60)
        plt.scatter(*min_res.x, color='black', marker='*', s=25, label='Local Min')
        plt.scatter(*bas_res.x, color='white', marker='^', s=60)
        plt.scatter(*bas_res.x, color='black', marker='^', s=25, label='Basin Min')

        plt.axhline(bounds[1][0], ls='--', color='black')
        plt.axhline(bounds[1][1], ls='--', color='black')
        plt.axvline(bounds[0][0], ls='--', color='black')
        plt.axvline(bounds[0][1], ls='--', color='black')

        # plt.axhline(0.2, ls=':', color='black')
        # plt.axhline(4.5, ls=':', color='black')
        # plt.axvline(0.0, ls=':', color='black')
        # plt.axvline(0.2, ls=':', color='black')
        # amps = ['025', '035', '045', '175', '225', '25']
        # for amp in amps:
        #     plt.axvline(float('0.' + amp), color='black', ls=':')
        # spreads = ['225', '275', '325', '375', '5']
        # for spread in spreads:
        #     plt.axhline(float('0.' + spread) * 10, color='black', ls=':')
        # plt.axhline(float('0.' + spread) * 10, color='black', ls=':', label='Simulation in Progress')

        plt.xlabel('amp')
        plt.ylabel('spread')
        plt.legend()
        fig2.colorbar(pcm2)
        fig2.tight_layout()
        plt.subplots_adjust(left=0.084, right=1, bottom=0.096, top=0.986)

        fig2_data.update({energy: {'X': X_intp, 'Y': Y_intp, 'Z': Z_intp, 'basin_min': bas_res.x}})

        for z in [Z, Z_intp]:
            z_min = np.min(z) if np.min(z) < z_min else z_min
            z_max = np.max(z) if np.max(z) > z_max else z_max

        fig3 = plt.figure()
        fig3.canvas.manager.set_window_title(f'{energy} GeV Vs Amp')
        spreads = np.unique(df['spread'])
        color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
        amp_mins_per_spread.append([])
        # amp_min_sds_per_spread.append([])
        energy_spreads.append(spreads)
        bs_jobs = []
        for spread in spreads:
            c = next(color)
            df_s = df[df['spread'] == spread]
            plt.scatter(df_s['amp'], np.log10(df_s['chi2_avg']), color=c, label=f'spread {spread}')
            zs = []
            amps = np.linspace(df_s['amp'].min(), df_s['amp'].max(), 10000)
            for amp in amps:
                zs.append(f(amp, spread))
            plt.plot(amps, zs, color=c)
            f_1d_spread = interp.interp1d(df_s['amp'], np.log10(df_s['chi2_avg']), kind='cubic')
            # print(df_s['amp'])
            # min_amp = minimize(f_1d_spread, 0.5, bounds=(bounds[0],))
            bas_bnds = BasinBounds(xmax=bounds[0][-1], xmin=bounds[0][0])
            # min_kwargs = {'method': 'L-BFGS-B', 'bounds': (bounds[0],)}
            # bas_min_amp = basinhopping(f_1d_spread, 0.5, minimizer_kwargs=min_kwargs,
            #                            accept_test=bas_bnds, niter=10, stepsize=1, T=1.2)
            min_kwargs = {'method': 'L-BFGS-B', 'jac': True, 'bounds': (bounds[0],)}

            def f_jac_slice(amp_i):
                r = f_jac(amp_i, spread)
                return r[0], r[1][0]

            bas_min_amp = basinhopping(f_jac_slice, df_s['amp'][df_s['chi2_avg'].idxmin()], minimizer_kwargs=min_kwargs,
                                       accept_test=bas_bnds, niter=30, stepsize=0.4, T=0.2)
            amp_mins_per_spread[-1].append(bas_min_amp.x[0])

            bs_jobs.append((f_jac_bss, bs_cols, df_s, min_kwargs, bas_bnds, spread))
            jobs = []
            for f_jac_bs, bs_col in zip(f_jac_bss, bs_cols):
                # def f_jac_slice_bs(amp_i):
                #     r = f_jac_bs(amp_i, spread)
                #     return r[0], r[1][0]
                jobs.append((f_jac_bs, df_s['amp'][df_s[bs_col].idxmin()], min_kwargs, bas_bnds, spread))

            # print('Getting bootstrap minima')
            # bs_mins = []
            # with Pool(threads) as pool:
            #     for df_entry in tqdm.tqdm(pool.istarmap(basin_min, jobs), total=len(jobs)):
            #         bs_mins.append(df_entry)
            #
            #     # bas_min_amp_bs = basinhopping(f_jac_slice_bs, df_s['amp'][df_s[bs_col].idxmin()], stepsize=0.4, T=0.2,
            #     #                               minimizer_kwargs=min_kwargs, accept_test=bas_bnds, niter=30)
            #     # bs_mins.append(bas_min_amp_bs.x[0])
            # amp_min_sds_per_spread[-1].append(np.std(bs_mins))

            # plt.plot(amps, f_1d_spread(amps), color=c, ls='--', alpha=0.6)
            # plt.scatter(bas_min_amp.x[0], f_1d_spread(bas_min_amp.x[0]), color=c, marker='P')
            def f_slice(amp_i):
                return f(amp_i, spread)

            plt.plot(amps, f_slice(amps), color=c, ls='--', alpha=0.6)
            plt.scatter(bas_min_amp.x[0], f_slice(bas_min_amp.x[0]), color=c, marker='P')
        plt.axhline(bas_res.fun, ls='--', color='gray')
        plt.xlabel('amp')
        plt.ylabel('chi2_avg')
        plt.legend(bbox_to_anchor=(1.05, 1))
        fig3.tight_layout()

        print('Getting bootstrap minima')
        bs_spreads, bs_mins = [], []
        with Pool(threads) as pool:
            for bs_spread, bs_min in tqdm.tqdm(pool.istarmap(basin_min2, bs_jobs), total=len(bs_jobs)):
                bs_spreads.append(bs_spread)
                bs_mins.append(bs_min)
        bs_spreads, bs_mins = zip(*sorted(zip(bs_spreads, bs_mins)))
        amp_min_sds_per_spread.append(bs_mins)

        fig4 = plt.figure()
        fig4.canvas.manager.set_window_title(f'{energy} GeV Vs Spread')
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

    fig_mins = plt.figure()
    fig_mins.canvas.manager.set_window_title(f'Mins in Amp-Spread Space')
    plt.grid()
    # plt.scatter(min_amps, min_spreads)
    # plt.xlim(bounds[0])
    # plt.ylim(bounds[1])
    plt.xlim((0, 0.05))
    plt.ylim((0, 2.5))
    plt.xlabel('amp')
    plt.ylabel('spread')
    plt.title(data_set)
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies))))
    for x, y, am, ame, es, e in \
            zip(min_amps, min_spreads, amp_mins_per_spread, amp_min_sds_per_spread, energy_spreads, energies):
        c = next(color)
        am, ame = np.array(am), np.array(ame)
        plt.scatter(x, y, marker='+', color=c)
        plt.annotate(f'{e}', (x, y), textcoords='offset points', xytext=(0, 10), ha='center', color=c)
        plt.plot(am, es, color=c, alpha=0.7, label=f'{e}GeV')
        plt.fill_betweenx(es, am - ame, am + ame, alpha=0.2, color=c)

    avgs_per_spread = []
    for spread_index in range(len(energy_spreads[0])):
        e_sum = 0
        for energy_mins in amp_mins_per_spread:
            e_sum += energy_mins[spread_index]
        avgs_per_spread.append(e_sum / len(energies))

    plt.plot(avgs_per_spread, energy_spreads[0], ls='--', color='black', label='Average')

    fig_mins.legend()
    fig_mins.tight_layout()
    # fig3_data.update({stat: [min_amps, min_spreads]})

    fig_mins_over_avg = plt.figure()
    fig_mins_over_avg.canvas.manager.set_window_title('Amp Mins divided by Average')
    plt.grid()
    plt.xlabel('spread')
    plt.ylabel('minimum amplitude / energy averaged minimum amplitude')
    color = iter(cm.rainbow(np.linspace(0, 1, len(energies))))
    for energy_mins, energy_mins_sd, energy, spreads in \
            zip(amp_mins_per_spread, amp_min_sds_per_spread, energies, energy_spreads):
        c = next(color)
        plt.plot(spreads, np.array(energy_mins) / np.array(avgs_per_spread), color=c, label=f'{energy}GeV')
        plt.fill_between(spreads, (np.array(energy_mins) - energy_mins_sd) / np.array(avgs_per_spread),
                         (np.array(energy_mins) + energy_mins_sd) / np.array(avgs_per_spread), color=c, alpha=0.2)
    plt.axhline(1, color='black', ls='--')
    fig_mins_over_avg.legend()
    fig_mins_over_avg.tight_layout()

    fig_mins_vs_e = plt.figure()
    fig_mins_vs_e.canvas.manager.set_window_title(f'Amp Mins vs Energy')
    plt.grid()
    # plt.scatter(min_amps, min_spreads)
    plt.xlim((0, 80))
    plt.ylim((0, 0.025))
    plt.xlabel('energy')
    plt.ylabel('amp_min')
    plt.title(data_set)
    spreads = energy_spreads[0]  # Assume all energies have all spreads
    max_spread = 2.5
    color = iter(cm.rainbow(np.linspace(0, 1, len([x for x in spreads if x <= max_spread]))))
    inv_amp_mins = np.array(amp_mins_per_spread).T
    for spread, min_amps in zip(spreads, inv_amp_mins):
        if spread > max_spread:
            continue
        c = next(color)
        plt.plot(energies, min_amps, marker='o', color=c, alpha=0.7, label=f'spread {spread:.2f}')
    fig_mins_vs_e.legend()
    fig_mins_vs_e.tight_layout()

    fig1_all, axes1 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig1_all.canvas.manager.set_window_title(f'Raw 2D Chi2')
    axes1 = np.ravel(axes1)
    for i, energy in enumerate(fig1_data):
        cbar = axes1[i].pcolormesh(fig1_data[energy]['X'], fig1_data[energy]['Y'], fig1_data[energy]['Z'],
                                   cmap='jet', vmin=z_min, vmax=z_max)
        axes1[i].scatter(fig1_data[energy]['x'], fig1_data[energy]['y'], color='black', marker='o',
                         s=(72. / fig1_all.dpi) ** 2)
        axes1[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes1[3].set_xlabel('amp')
    axes1[4].set_xlabel('amp')
    axes1[5].set_xlabel('amp')
    axes1[0].set_ylabel('spread')
    axes1[3].set_ylabel('spread')
    # fig1_all.colorbar(cbar)
    fig1_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)
    # fig1_all.tight_layout()

    fig2_all, axes2 = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 8))
    fig2_all.canvas.manager.set_window_title(f'Interpolated 2D Chi2')
    axes2 = np.ravel(axes2)
    for i, energy in enumerate(fig2_data):
        cbar = axes2[i].pcolormesh(fig2_data[energy]['X'], fig2_data[energy]['Y'], fig2_data[energy]['Z'],
                                   cmap='jet',
                                   vmin=z_min, vmax=z_max)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='white', marker='^', s=60)
        axes2[i].scatter(*fig2_data[energy]['basin_min'], color='black', marker='^', s=25, label='Basin Min')
        axes2[i].text(0.3, 1, f'{energy}GeV', color='white')
    axes2[3].set_xlabel('amp')
    axes2[4].set_xlabel('amp')
    axes2[5].set_xlabel('amp')
    axes2[0].set_ylabel('spread')
    axes2[3].set_ylabel('spread')
    # fig2_all.colorbar(cbar)
    fig2_all.subplots_adjust(top=0.99, bottom=0.06, left=0.04, right=0.99, wspace=0.04, hspace=0.03)

    plt.show()


def basin_min(f, x0, min_kwargs, bas_bnds, spread, niter=30, stepsize=0.4, T=0.2):
    def f_jac_slice_bs(amp_i):
        r = f(amp_i, spread)
        return r[0], r[1][0]
    bas_min_amp_bs = basinhopping(f_jac_slice_bs, x0, stepsize=stepsize, T=T,
                                  minimizer_kwargs=min_kwargs, accept_test=bas_bnds, niter=niter)
    return bas_min_amp_bs.x[0]


def basin_min2(f_jac_bss, bs_cols, df_s, min_kwargs, bas_bnds, spread, niter=30, stepsize=0.4, T=0.2):
    bs_mins = []
    for f_jac_bs, bs_col in zip(f_jac_bss, bs_cols):
        def f_jac_slice_bs(amp_i):
            r = f_jac_bs(amp_i, spread)
            return r[0], r[1][0]
        bas_min_amp_bs = basinhopping(f_jac_slice_bs, df_s['amp'][df_s[bs_col].idxmin()], stepsize=0.4, T=0.2,
                                      minimizer_kwargs=min_kwargs, accept_test=bas_bnds, niter=30)
        bs_mins.append(bas_min_amp_bs.x[0])
    return spread, np.std(bs_mins)
    # amp_min_sds_per_spread[-1].append(np.std(bs_mins))


def gaus2d_test():
    x = np.linspace(0, 10, 50)
    y = np.linspace(-6, -1, 10)
    xv, yv = np.meshgrid(x, y)
    xx = np.ravel(xv)
    yy = np.ravel(yv)

    mu = np.array([2, -2])
    sigma = np.array([[0.5, 0], [0, 1]])

    # print(np.column_stack((xx, yy)))

    z = np.array([gaus(xyi, mu, sigma) for xyi in np.column_stack((xx, yy))])
    # print(z)
    print(z.shape)
    print(z.reshape((len(x), len(y))))

    plt.figure()
    plt.title('Original')
    plt.imshow(z.reshape((len(y), len(x))), extent=[min(x), max(x), min(y), max(y)], aspect='auto', origin='lower',
               cmap='jet')

    print(gaus(np.array([2, -3]), mu, sigma))

    tcks = interp.bisplrep(xx, yy, z, kx=1, ky=1)

    # f = interp2d(df_chi['amp'], df_chi['spread'], df_chi['chi2_sum'], kind='cubic')
    x1 = np.linspace(min(x), max(x), 10)
    y1 = np.linspace(min(y), max(y), 5)
    xx1, yy1 = np.meshgrid(x1, y1)
    # z = f(x, y)
    z1 = interp.bisplev(x1, y1, tcks)
    print(z1)
    print(z1.T)
    print(z1.shape)
    print(z1.T.shape)

    plt.figure()
    plt.title('Interpolated')
    plt.imshow(z1.T, extent=[min(x1), max(x1), min(y1), max(y1)], aspect='auto', origin='lower', cmap='jet')

    plt.show()


def gaus(x, mu, sigma):
    return np.exp(-0.5 * np.dot(np.dot((x - mu), np.linalg.inv(sigma)), (x - mu).T))


def scipy_doc_test_edit():
    x = np.arange(0, 4, 0.25)
    y = np.array([0, 2, 3, 3.5, 3.75, 3.875, 3.9375, 4])
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.pi * X / 2) * np.exp(Y / 2)

    x2 = np.arange(0, 4, 0.05)
    y2 = np.arange(0, 4, 0.05)
    print('x: ', x)
    print('y: ', y)
    print('Z: ', Z)
    f = interp.interp2d(x, y, Z, kind='linear')
    Z2 = f(x2, y2)

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].pcolormesh(X, Y, Z, cmap='jet')

    X2, Y2 = np.meshgrid(x2, y2)
    ax[1].pcolormesh(X2, Y2, Z2, cmap='jet')

    fig2, ax2 = plt.subplots()
    x2_index = dict((val, key) for key, val in dict(enumerate(x2)).items())
    ax2.grid()
    color = iter(cm.rainbow(np.linspace(0, 1, len(x))))
    for x_index, x_val in enumerate(x):
        c = next(color)
        ax2.scatter(y, Z.T[x_index], label=f'x: {x_val}', color=c)
        z_xi = []
        for yi in y:
            # print(xi, y_val, f(xi, y_val)[0])
            z_xi.append(f(x_val, yi)[0])
        # print(f(x, np.full(len(x), y_val)))
        ax2.plot(y, z_xi, color=c)
    ax2.legend()
    ax2.set_xlabel('y')

    fig3, ax3 = plt.subplots()
    x2_index = dict((val, key) for key, val in dict(enumerate(x2)).items())
    ax3.grid()
    color = iter(cm.rainbow(np.linspace(0, 1, len(y))))
    for y_index, y_val in enumerate(y):
        c = next(color)
        ax3.scatter(x, Z[y_index], label=f'y: {y_val}', color=c)
        print(f'y={y_val}')
        z_yi = []
        for xi in x:
            print(xi, y_val, f(xi, y_val)[0])
            z_yi.append(f(xi, y_val)[0])
        # print(f(x, np.full(len(x), y_val)))
        ax3.plot(x, z_yi, color=c)
    ax3.legend()
    ax3.set_xlabel('y')
    print(x)
    print(x2)

    plt.show()


def scipy_doc_test():
    x = np.linspace(0, 4, 13)
    y = np.array([0, 2, 3, 3.5, 3.75, 3.875, 3.9375, 4])
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.pi * X / 2) * np.exp(Y / 2)

    x2 = np.linspace(0, 4, 65)
    y2 = np.linspace(0, 4, 65)
    f = interp.interp2d(x, y, Z, kind='cubic')
    Z2 = f(x2, y2)

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].pcolormesh(X, Y, Z)

    X2, Y2 = np.meshgrid(x2, y2)
    ax[1].pcolormesh(X2, Y2, Z2)

    plt.show()


if __name__ == '__main__':
    main()
