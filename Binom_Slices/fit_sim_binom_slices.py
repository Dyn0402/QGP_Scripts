#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 25 1:30 PM 2022
Created in PyCharm
Created as QGP_Scripts/fit_sim_binom_slices.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from multiprocessing import Pool
import tqdm


def main():
    pars = set_pars()
    df = get_df(pars)
    get_chi2s(df, pars)
    print('donzo')


def set_pars():
    pars = {'thrads': 15,
            'base_path': 'D:/Research/Results/Azimuth_Analysis/',
            'df_name': 'binom_slice_sds_cent8.csv',
            'chi_df_name': 'chi_df_ampt_sds_cent8.csv',
            'energies': [62]}

    return pars


def get_df(pars):
    df = pd.read_csv(pars['base_path'] + pars['df_name'])
    df = df.dropna()
    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    return df


def get_chi2s(df, pars):
    df = pre_filter_df(df)
    jobs = []
    for i in range(10):  # Placeholder for jobs
        jobs.append(i)

    chi2_sets = []
    with Pool(pars['threads']) as pool:
        for chi2_set in tqdm.tqdm(pool.istarmap(chi2_vs_protons, jobs), total=len(jobs)):
            chi2_sets.extend(chi2_set)


def pre_filter_df(df):
    pass


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


if __name__ == '__main__':
    main()
