#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 21 8:32 PM 2023
Created in PyCharm
Created as QGP_Scripts/bes_v2_plot

@author: Dylan Neff, Dylan
"""
import os.path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

from analyze_binom_slices import read_flow_values


def main():
    # plot_star_models()
    plot_star_sys()

    print('donzo')


def plot_star_models():
    data_sets = {'BES1': 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                         'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
                 'AMPT': 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/'
                         'Ampt_rapid05_resample_norotate_epbins1_0/',
                 'CF': 'F:/Research/Data_CF/default_resample_epbins1/CF_rapid05_resample_norotate_epbins1_0/',
                 'CFEV': 'F:/Research/Data_CFEV/default_resample_epbins1/CFEV_rapid05_resample_norotate_epbins1_0/',
                 'CFEVb342': 'F:/Research/Data_CFEVb342/default_resample_epbins1/'
                             'CFEVb342_rapid05_resample_norotate_epbins1_0/',
                 }
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}
    energies = [7, 11, 19, 27, 39, 62]

    df = read_flow_data(data_sets, energies, cent_map)
    plot_flow_data(df, data_sets, energies, cent_map)


def plot_star_sys():
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}
    # energies = [7, 11, 19, 27, 39, 62]
    energies = [7, 11, 19, 27, 39, 62]
    sys_names = ['DCA', 'nHitsFit', 'tof mass^2 range', 'nsigma proton range']
    x_labels = ['DCA (cm)', 'Minimum nHitsFit', 'mass^2 range (GeV)', '|nsigmaproton| max']
    data_sets_vals = [
        [1, 0.5, 0.8, 1.2, 1.5],
        [15, 20, 25],
        [2, 4, 6, 8, 10],
        [1.5, 1.8, 2.0, 2.2, 2.5]
    ]
    def_vals = [1, 20, 6, 2]
    sys_vals = [0.5, 15, 2, 1.5]
    data_sets = [
        {
            'dca = 1.0 cm': 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                            'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
            'dca = 0.5 cm': 'F:/Research/Data/v2_sys/'
                            'rapid05_resample_norotate_dca05_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'dca = 0.8 cm': 'F:/Research/Data/v2_sys/'
                            'rapid05_resample_norotate_dca08_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'dca = 1.2 cm': 'F:/Research/Data/v2_sys/'
                            'rapid05_resample_norotate_dca12_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'dca = 1.5 cm': 'F:/Research/Data/v2_sys/'
                            'rapid05_resample_norotate_dca15_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
        },
        {
            'nHitsFit = 20': 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                             'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
            'nHitsFit = 15': 'F:/Research/Data/v2_sys/'
                             'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit15_epbins1_calcv2_qaonly_0/',
            'nHitsFit = 25': 'F:/Research/Data/v2_sys/'
                             'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit25_epbins1_calcv2_qaonly_0/',
        },
        {
            'm^2 range = 0.6 GeV':
                'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
            'm^2 range = 0.2 GeV':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r2_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'm^2 range = 0.4 GeV':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r4_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'm^2 range = 0.8 GeV':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r8_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            'm^2 range = 1.0 GeV':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r10_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
        },
        {
            '|nsigmaproton| max = 2.0':
                'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
            '|nsigmaproton| max = 1.5':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx075_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            '|nsigmaproton| max = 1.8':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx09_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            '|nsigmaproton| max = 2.2':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx11_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
            '|nsigmaproton| max = 2.5':
                'F:/Research/Data/v2_sys/'
                'rapid05_resample_norotate_dca1_nsprx125_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_0/',
        }
    ]

    for sys_name, x_label, data_set, data_set_vals, def_val, sys_val in \
            zip(sys_names, x_labels, data_sets, data_sets_vals, def_vals, sys_vals):
        pdf_path = f'F:/Research/Results/BES_v2_Systematics/BES1_proton_v2_systematic_plots_{sys_name}.pdf'
        data_set_vals = dict(zip(data_set.keys(), data_set_vals))

        df = read_flow_data(data_set, energies, cent_map)
        # plot_flow_data(df, data_sets, energies, cent_map)
        plot_sys_data(df, energies, cent_map, data_set_vals, def_val, sys_val, sys_name, x_label, pdf_path)


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


def plot_flow_data(df, data_sets, energies, cent_map):
    # Plots vs centrality
    for data_set in data_sets.keys():
        df_set = df[df['data_set'] == data_set]
        if len(pd.unique(df_set['cent_bin'])) <= 1:
            continue
        fig_v2_vs_cent, ax_v2_vs_cent = plt.subplots(dpi=144)
        plt.xticks(rotation=30)
        ax_v2_vs_cent.grid()
        ax_v2_vs_cent.axhline(0, color='black')
        ax_v2_vs_cent.set_title(f'{data_set} v2')
        fig_res_vs_cent, ax_res_vs_cent = plt.subplots(dpi=144)
        plt.xticks(rotation=30)
        ax_res_vs_cent.grid()
        ax_res_vs_cent.axhline(0, color='black')
        ax_res_vs_cent.set_title(f'{data_set} Resolution')
        for energy in energies:
            df_e = df_set[df_set['energy'] == energy]
            df_e = df_e.sort_values('cent_bin', ascending=False)
            ax_v2_vs_cent.errorbar(df_e['cent_label'], df_e['v2_val'], df_e['v2_err'], ls='-', marker='o',
                                   alpha=0.8, label=f'{energy}GeV')
            ax_res_vs_cent.errorbar(df_e['cent_label'], df_e['res_val'], df_e['res_err'], ls='-', marker='o',
                                    alpha=0.8, label=f'{energy}GeV')
        ax_v2_vs_cent.legend()
        ax_res_vs_cent.legend()
        ax_v2_vs_cent.set_xlabel('Centrality')
        ax_res_vs_cent.set_xlabel('Centrality')
        ax_v2_vs_cent.set_ylabel('v2')
        ax_res_vs_cent.set_ylabel('Event Plane Resolution')
        ax_v2_vs_cent.set_ylim(-0.005, 0.17)
        ax_res_vs_cent.set_ylim(-0.02, 0.53)
        fig_v2_vs_cent.canvas.manager.set_window_title(f'{data_set} v2')
        fig_res_vs_cent.canvas.manager.set_window_title(f'{data_set} Resolution')
        fig_v2_vs_cent.tight_layout()
        fig_res_vs_cent.tight_layout()

    # Plots vs energy
    for cent in cent_map.keys():
        if cent not in pd.unique(df['cent_bin']):
            continue
        df_cent = df[df['cent_bin'] == cent]
        fig_v2_vs_energy, ax_v2_vs_energy = plt.subplots(dpi=144)
        plt.xticks(rotation=30)
        ax_v2_vs_energy.grid()
        ax_v2_vs_energy.axhline(0, color='black')
        ax_v2_vs_energy.set_title(f'{cent_map[cent]} Centrality v2')
        fig_res_vs_energy, ax_res_vs_energy = plt.subplots(dpi=144)
        plt.xticks(rotation=30)
        ax_res_vs_energy.grid()
        ax_res_vs_energy.axhline(0, color='black')
        ax_res_vs_energy.set_title(f'{cent_map[cent]} Centrality Resolution')
        for data_set in data_sets.keys():
            df_set = df_cent[df_cent['data_set'] == data_set]
            df_set = df_set.sort_values('energy')
            ax_v2_vs_energy.errorbar(df_set['energy'], df_set['v2_val'], df_set['v2_err'], ls='-', marker='o',
                                     alpha=0.8, label=data_set)
            ax_res_vs_energy.errorbar(df_set['energy'], df_set['res_val'], df_set['res_err'], ls='-', marker='o',
                                      alpha=0.8, label=data_set)
        ax_v2_vs_energy.legend()
        ax_res_vs_energy.legend()
        ax_v2_vs_energy.set_xlabel('Energy (GeV)')
        ax_res_vs_energy.set_xlabel('Energy (GeV)')
        ax_v2_vs_energy.set_ylabel('v2')
        ax_res_vs_energy.set_ylabel('Event Plane Resolution')
        fig_v2_vs_energy.canvas.manager.set_window_title(f'{cent_map[cent]} Centrality v2')
        fig_res_vs_energy.canvas.manager.set_window_title(f'{cent_map[cent]} Centrality Resolution')
        fig_v2_vs_energy.tight_layout()
        fig_res_vs_energy.tight_layout()

    plt.show()


def plot_sys_data(df, energies, cent_map, data_set_vals, def_val, sys_val, sys_name=None, xlabel=None, pdf_path=None):
    # Plots vs sys variable
    pdf = PdfPages(pdf_path) if pdf_path is not None else None
    # Plot energies together
    for cent_bin in sorted(pd.unique(df['cent_bin']), reverse=True):
        df_cent = df[df['cent_bin'] == cent_bin]
        fig, ax = plt.subplots(dpi=144)
        ax.grid()
        sys_name = 'Systematic' if sys_name is None else sys_name
        title = f'{sys_name} {cent_map[cent_bin]}'
        ax.set_title(title)
        ax.set_xlabel('Systematic' if xlabel is None else xlabel)
        ax.set_ylabel('v2')
        for energy in energies:
            # print(f'{sys_name} {energy}GeV {cent_map[cent_bin]}')
            # if sys_name == 'nHitsFit':
            #     print(df_cent)
            df_e = df_cent[df_cent['energy'] == energy]
            df_e = df_e.assign(data_set_vals=[data_set_vals[data_set] for data_set in df_e['data_set']])
            df_e = df_e.sort_values('data_set_vals')
            ebar = ax.errorbar(df_e['data_set_vals'], df_e['v2_val'], yerr=df_e['v2_err'], ls='-', marker='o',
                               alpha=0.7, label=f'{energy} GeV')
            def_val_v2, def_err = df_e[df_e['data_set_vals'] == def_val].iloc[0][['v2_val', 'v2_err']]
            var_val_v2, var_err = df_e[df_e['data_set_vals'] == sys_val].iloc[0][['v2_val', 'v2_err']]
            barlow = (def_val_v2 - var_val_v2) ** 2 - abs(def_err ** 2 - var_err ** 2)
            barlow = 0 if barlow < 0 else np.sqrt(barlow / 12.0)
            print(f'{energy}GeV {cent_map[cent_bin]} barlow: {barlow}')
            ax.errorbar(def_val, def_val_v2, yerr=barlow, elinewidth=8, alpha=0.3, color=ebar[0].get_color())
        ax.legend()
        fig.canvas.manager.set_window_title(title)
        fig.tight_layout()
        if pdf:
            pdf.savefig(fig)

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
            # print(f'{energy}GeV {cent_map[cent_bin]} barlow: {barlow}')
            ax.errorbar(def_val, def_val_v2, yerr=barlow, elinewidth=8, alpha=0.3, color=ebar[0].get_color())
            fig.canvas.manager.set_window_title(title)
            fig.tight_layout()
            if pdf:
                pdf.savefig(fig)

    if pdf:
        pdf.close()

    # plt.show()


if __name__ == '__main__':
    main()
