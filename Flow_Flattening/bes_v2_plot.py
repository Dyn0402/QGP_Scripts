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
import pandas as pd

from analyze_binom_slices import read_flow_values


def main():
    data_sets = {'BES1': 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/'
                         'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/',
                 'AMPT': 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/'
                         'Ampt_rapid05_resample_norotate_epbins1_0/',
                 'CF': 'F:/Research/Data_CF/default_resample_epbins1/CF_rapid05_resample_norotate_epbins1_0/',
                 'CFEV': 'F:/Research/Data_CFEV/default_resample_epbins1/CFEV_rapid05_resample_norotate_epbins1_0/',
                 # 'CFEVb342': 'F:/Research/Data_CFEVb342/default_resample_epbins1/'
                 #              'CFEVb342_rapid05_resample_norotate_epbins1_0/',
                 }
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}
    energies = [7, 11, 19, 27, 39, 62]

    # energy_figs, energy_axs = list(zip(*[plt.subplots(dpi=144) for energy in energies]))
    # energy_figs, energy_axs = [dict(zip(energies, x)) for x in [energy_figs, energy_axs]]
    # for ax in energy_axs.values():
    #     ax.grid()
    #     ax.axhline(0, color='black')
    #     ax.set_title('')

    df = []
    for name, directory in data_sets.items():
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

    df = pd.DataFrame(df)
    df = df[df['cent_label'] != '80-90%']

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

    print('donzo')


if __name__ == '__main__':
    main()
