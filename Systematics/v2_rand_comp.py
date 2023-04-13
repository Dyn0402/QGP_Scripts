#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 13 11:36 AM 2023
Created in PyCharm
Created as QGP_Scripts/v2_rand_comp

@author: Dylan Neff, Dylan
"""

import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def main():
    energies = [7]
    data_set = {
        'run 1': 'F:/Research/Data/default_rand_sys_test/'
                 'rapid05_resample_norotate_strefnoseed_nomix_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/',
        'run 2': 'F:/Research/Data/default_rand_sys_test/'
                 'rapid05_resample_norotate_strefnoseed_nomix_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_1/',
        'run 3': 'F:/Research/Data/default_rand_sys_test/'
                 'rapid05_resample_norotate_strefnoseed_nomix_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_2/',
        'run 4': 'F:/Research/Data/default_rand_sys_test/'
                 'rapid05_resample_norotate_strefnoseed_nomix_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_3/',
    }
    data_set_vals = [1, 2, 3, 4]
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    data_set_vals = dict(zip(data_set.keys(), data_set_vals))

    df = read_flow_data(data_set, energies, cent_map)

    plot_sys_data(df, energies, cent_map, data_set_vals, 1, 2)

    print('donzo')


def read_flow_data(data_sets, energies, cent_map):
    df = []
    for name, directory in data_sets.items():
        print(name, directory)
        for energy in energies:
            v2_file_path = f'{directory}{energy}GeV/v2.txt'
            if os.path.exists(v2_file_path):
                with open(v2_file_path, 'r') as file:
                    lines = file.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        df.append({'data_set': name, 'energy': energy,
                                   'cent_bin': int(line[0]), 'cent_label': cent_map[int(line[0])],
                                   'v2_val': float(line[1].split()[0]), 'v2_err': float(line[1].split()[1]),
                                   'res_val': float(line[2].split()[0]), 'res_err': float(line[2].split()[1])})
            else:
                print(f'{v2_file_path} does not exist!')

    df = pd.DataFrame(df)
    df = df[df['cent_label'] != '80-90%']

    return df


def plot_sys_data(df, energies, cent_map, data_set_vals, def_val, sys_val, sys_name=None, xlabel=None):
    # Plots vs sys variable
    # pdf = PdfPages(pdf_path) if pdf_path is not None else None
    # Plot energies together
    # for cent_bin in sorted(pd.unique(df['cent_bin']), reverse=True):
    #     df_cent = df[df['cent_bin'] == cent_bin]
    #     fig, ax = plt.subplots(dpi=144)
    #     ax.grid()
    #     sys_name = 'Systematic' if sys_name is None else sys_name
    #     title = f'{sys_name} {cent_map[cent_bin]}'
    #     ax.set_title(title)
    #     ax.set_xlabel('Systematic' if xlabel is None else xlabel)
    #     ax.set_ylabel('v2')
    #     for energy in energies:
    #         # print(f'{sys_name} {energy}GeV {cent_map[cent_bin]}')
    #         # if sys_name == 'nHitsFit':
    #         #     print(df_cent)
    #         df_e = df_cent[df_cent['energy'] == energy]
    #         df_e = df_e.assign(data_set_vals=[data_set_vals[data_set] for data_set in df_e['data_set']])
    #         df_e = df_e.sort_values('data_set_vals')
    #         ebar = ax.errorbar(df_e['data_set_vals'], df_e['v2_val'], yerr=df_e['v2_err'], ls='-', marker='o',
    #                            alpha=0.7, label=f'{energy} GeV')
    #         def_val_v2, def_err = df_e[df_e['data_set_vals'] == def_val].iloc[0][['v2_val', 'v2_err']]
    #         var_val_v2, var_err = df_e[df_e['data_set_vals'] == sys_val].iloc[0][['v2_val', 'v2_err']]
    #         barlow = (def_val_v2 - var_val_v2) ** 2 - abs(def_err ** 2 - var_err ** 2)
    #         barlow = 0 if barlow < 0 else np.sqrt(barlow / 12.0)
    #         print(f'{energy}GeV {cent_map[cent_bin]} barlow: {barlow}')
    #         ax.errorbar(def_val, def_val_v2, yerr=barlow, elinewidth=8, alpha=0.3, color=ebar[0].get_color())
    #     ax.legend()
    #     fig.canvas.manager.set_window_title(title)
    #     fig.tight_layout()
    #     # if pdf:
    #     #     pdf.savefig(fig)

    for cent_bin in sorted(pd.unique(df['cent_bin']), reverse=True):
        df_cent = df[df['cent_bin'] == cent_bin]
        for energy in energies:
            fig, ax = plt.subplots(dpi=144)
            ax.grid()
            sys_name = 'Systematic' if sys_name is None else sys_name
            title = f'{sys_name} {energy}GeV {cent_map[cent_bin]}'
            ax.set_title(title)
            ax.set_xlabel('Systematic' if xlabel is None else xlabel)
            ax.set_ylabel('v2')
            df_e = df_cent[df_cent['energy'] == energy]
            df_e = df_e.assign(data_set_vals=[data_set_vals[data_set] for data_set in df_e['data_set']])
            df_e = df_e.sort_values('data_set_vals')
            ebar = ax.errorbar(df_e['data_set_vals'], df_e['v2_val'], yerr=df_e['v2_err'], ls='-', marker='o',
                               alpha=0.7)
            def_val_v2, def_err = df_e[df_e['data_set_vals'] == def_val].iloc[0][['v2_val', 'v2_err']]
            var_val_v2, var_err = df_e[df_e['data_set_vals'] == sys_val].iloc[0][['v2_val', 'v2_err']]
            barlow = (def_val_v2 - var_val_v2) ** 2 - abs(def_err ** 2 - var_err ** 2)
            barlow = 0 if barlow < 0 else np.sqrt(barlow / 12.0)
            print(f'{energy}GeV {cent_map[cent_bin]} barlow: {barlow}')
            ax.errorbar(def_val, def_val_v2, yerr=barlow, elinewidth=8, alpha=0.3, color=ebar[0].get_color())
            std_err = np.std(df_e['v2_val'])
            ax.errorbar(sys_val, var_val_v2, yerr=std_err, elinewidth=8, alpha=0.3, color=ebar[0].get_color())
            fig.canvas.manager.set_window_title(title)
            fig.tight_layout()
    #         if pdf:
    #             pdf.savefig(fig)

    # if pdf:
    #     pdf.close()

    plt.show()


if __name__ == '__main__':
    main()
