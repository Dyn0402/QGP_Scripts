#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 06 9:10 PM 2022
Created in PyCharm
Created as QGP_Scripts/analysis_step_visualizations

@author: Dylan Neff, Dyn04
"""

import numpy as np
import matplotlib.pyplot as plt

from compare_dists import get_norm_dists


def main():
    # bs_visual()
    # full_bs_visual()
    full_bs_div_visual()
    print('donzo')


def bs_visual():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 60

    total_proton = 21

    rmm_bs_plot_frac = 0

    data_set = (base_path, 'default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0',
                energy, cent, div, total_proton, 'Data', 'Data_Mix', 'bes')
    sim_par = ('0125', '08')
    sim_set = (base_path, f'flat80_anticlmulti_spread{sim_par[1]}_amp{sim_par[0]}_resample',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
               62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = data_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist, raw_bs_tp_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist, mix_bs_tp_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    for rm_name, color, def_dists, bs_dists in \
            (('raw', 'blue', raw_tp_dist, raw_bs_tp_dists), ('mix', 'green', mix_tp_dist, mix_bs_tp_dists)):
        fig_raw, ax_raw = plt.subplots()
        ax_raw.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]), color='black', alpha=0.2,
                    label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists[1:]:
            ax_raw.plot(range(bs.size), bs / np.sum(bs), color='black', alpha=0.2)
        ax_raw.plot(range(def_dists.size), def_dists / np.sum(def_dists), alpha=1, color=color,
                    label=rm_name.capitalize())
        ax_raw.set_title(f'{rm_name.capitalize()} {name.capitalize()} {energy}GeV, {total_proton} Protons, '
                         f'{div} divs')
        ax_raw.axhline(0, ls='--', color='black', zorder=0)
        ax_raw.set_xlabel('Protons in Bin')
        ax_raw.legend()
        fig_raw.tight_layout()
        fig_raw.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} '
                                                f'{total_proton} Protons, {div} divs')

        fig_raw_diff, ax_raw_diff = plt.subplots()
        ax_raw_diff.axhline(0, ls='--', color=color)
        ax_raw_diff.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) -
                         def_dists / np.sum(def_dists), color='black',
                         alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_diff.plot(range(bs.size), bs / np.sum(bs) - def_dists / np.sum(def_dists), color='black',
                             alpha=0.2)
        ax_raw_diff.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_diff.set_xlabel('Protons in Bin')
        ax_raw_diff.legend()
        fig_raw_diff.tight_layout()
        fig_raw_diff.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                     f'Protons, {div} divs')

        fig_raw_div, ax_raw_div = plt.subplots()
        ax_raw_div.axhline(1, ls='--', color=color)
        ax_raw_div.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) /
                        (def_dists / np.sum(def_dists)), color='black',
                        alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_div.plot(range(bs.size), bs / np.sum(bs) / (def_dists / np.sum(def_dists)), color='black',
                            alpha=0.2)
        ax_raw_div.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_div.set_xlabel('Protons in Bin')
        ax_raw_div.legend()
        fig_raw_div.tight_layout()
        fig_raw_div.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    raw_sub_mix_def = raw_tp_dist / np.sum(raw_tp_dist) - mix_tp_dist / np.sum(mix_tp_dist)
    fig_raw_sub_mix, ax_raw_sub_mix = plt.subplots()
    ax_raw_sub_mix.axhline(0, ls='--', color='black')
    ax_raw_sub_mix.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='Raw - Mix', zorder=4)
    raw_sub_mix_bss = [raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_tp_dists for mix in mix_bs_tp_dists]
    raw_sub_mix_bss = np.array(raw_sub_mix_bss)
    bs_label = False
    for raw_sub_mix_bs in raw_sub_mix_bss:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix.plot(range(raw_sub_mix_bs.size), raw_sub_mix_bs, alpha=0.3, color='black')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    # ax_raw_sub_mix.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist))
    ax_raw_sub_mix.legend()
    ax_raw_sub_mix.set_title(f'Raw - Mix {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_sub_mix.tight_layout()
    fig_raw_sub_mix.canvas.manager.set_window_title(f'Raw - Mix {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    fig_raw_sub_mix_sd, ax_raw_sub_mix_sd = plt.subplots()
    ax_raw_sub_mix_sd.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sd.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='Raw - Mix')
    raw_sub_mix_sd = np.std(raw_sub_mix_bss, axis=0)
    # n_bands = 100
    # for sigma in np.linspace(0.0, 10, n_bands):
    #     ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - sigma * raw_sub_mix_sd,
    #                                    raw_sub_mix_def + sigma * raw_sub_mix_sd, color='red', alpha=1.0 / n_bands)
    ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - raw_sub_mix_sd,
                                   raw_sub_mix_def + raw_sub_mix_sd, color='red', alpha=0.4)
    # ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - 2 * raw_sub_mix_sd,
    #                                raw_sub_mix_def + 2 * raw_sub_mix_sd, color='red', alpha=0.4)
    ax_raw_sub_mix_sd.legend()
    ax_raw_sub_mix_sd.set_title(f'Raw - Mix SD {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_sub_mix_sd.tight_layout()
    fig_raw_sub_mix_sd.canvas.manager.set_window_title(f'Raw - Mix SD {name.capitalize()} {total_proton} '
                                                       f'Protons, {div} divs')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = sim_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist_sim, raw_bs_tp_dists_sim = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist_sim, mix_bs_tp_dists_sim = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    raw_sub_mix_def_sim = raw_tp_dist_sim / np.sum(raw_tp_dist_sim) - mix_tp_dist_sim / np.sum(mix_tp_dist_sim)
    raw_sub_mix_bss_sim = [raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_tp_dists_sim
                           for mix in mix_bs_tp_dists_sim]
    raw_sub_mix_sd_sim = np.std(np.array(raw_sub_mix_bss_sim), axis=0)

    fig_raw_sub_mix_sim_comp, ax_raw_sub_mix_sim_comp = plt.subplots()
    ax_raw_sub_mix_sim_comp.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='BES Raw - Mix')
    ax_raw_sub_mix_sim_comp.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - raw_sub_mix_sd,
                                         raw_sub_mix_def + raw_sub_mix_sd, color='red', alpha=0.4)
    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def_sim.size), raw_sub_mix_def_sim, color='purple',
                                 label='Sim Raw - Mix')
    ax_raw_sub_mix_sim_comp.fill_between(range(raw_tp_dist_sim.size), raw_sub_mix_def_sim - raw_sub_mix_sd_sim,
                                         raw_sub_mix_def_sim + raw_sub_mix_sd_sim, color='purple', alpha=0.4)
    ax_raw_sub_mix_sim_comp.legend()
    ax_raw_sub_mix_sim_comp.set_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_sub_mix_sim_comp.tight_layout()
    fig_raw_sub_mix_sim_comp.canvas.manager.set_window_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    diff = raw_sub_mix_def - raw_sub_mix_def_sim
    diff_sd = np.sqrt(raw_sub_mix_sd ** 2 + raw_sub_mix_sd_sim ** 2)
    fig_raw_sub_mix_sim_diff, ax_raw_sub_mix_sim_diff = plt.subplots()
    ax_raw_sub_mix_sim_diff.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_diff.plot(range(diff.size), diff, color='magenta', marker='o', label='BES - Sim')
    ax_raw_sub_mix_sim_diff.fill_between(range(diff_sd.size), diff - diff_sd, diff + diff_sd, alpha=0.4,
                                         color='magenta')
    ax_raw_sub_mix_sim_diff.legend()
    ax_raw_sub_mix_sim_diff.set_title(f'Raw - Mix BES Sim Diff {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_sub_mix_sim_diff.tight_layout()
    fig_raw_sub_mix_sim_diff.canvas.manager.set_window_title(f'Raw - Mix BES Sim Diff {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    chi2 = (raw_sub_mix_def - raw_sub_mix_def_sim) ** 2 / (raw_sub_mix_sd ** 2 + raw_sub_mix_sd_sim ** 2)
    fig_raw_sub_mix_sim_chi, ax_raw_sub_mix_sim_chi = plt.subplots()
    ax_raw_sub_mix_sim_chi.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_chi.plot(range(chi2.size), chi2, color='salmon', marker='o', label='BES - Sim Chi2')
    ax_raw_sub_mix_sim_chi.legend()
    ax_raw_sub_mix_sim_chi.set_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} {total_proton} '
                                     f'Protons, {div} divs')
    fig_raw_sub_mix_sim_chi.tight_layout()
    fig_raw_sub_mix_sim_chi.canvas.manager.set_window_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} '
                                                            f'{total_proton} Protons, {div} divs')

    plt.show()


def full_bs_visual():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 60

    total_proton = 21

    rmm_bs_plot_frac = 0.001

    data_set = (base_path, 'default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0',
                energy, cent, div, total_proton, 'Data', 'Data_Mix', 'bes')
    # sim_par = ('0125', '08')
    sim_par = ('0', '1')
    sim_set = (base_path, f'flat80_anticlmulti_spread{sim_par[1]}_amp{sim_par[0]}_resample',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
               62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = data_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist, raw_bs_tp_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist, mix_bs_tp_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    for rm_name, color, def_dists, bs_dists in \
            (('raw', 'blue', raw_tp_dist, raw_bs_tp_dists), ('mix', 'green', mix_tp_dist, mix_bs_tp_dists)):
        fig_raw, ax_raw = plt.subplots()
        ax_raw.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]), color='black', alpha=0.2,
                    label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists[1:]:
            ax_raw.plot(range(bs.size), bs / np.sum(bs), color='black', alpha=0.2)
        ax_raw.plot(range(def_dists.size), def_dists / np.sum(def_dists), alpha=1, color=color,
                    label=rm_name.capitalize())
        ax_raw.set_title(f'{rm_name.capitalize()} {name.capitalize()} {energy}GeV, {total_proton} Protons, '
                         f'{div} divs')
        ax_raw.axhline(0, ls='--', color='black', zorder=0)
        ax_raw.set_xlabel('Protons in Bin')
        ax_raw.legend()
        fig_raw.tight_layout()
        fig_raw.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} '
                                                f'{total_proton} Protons, {div} divs')

        fig_raw_diff, ax_raw_diff = plt.subplots()
        ax_raw_diff.axhline(0, ls='--', color=color)
        ax_raw_diff.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) -
                         def_dists / np.sum(def_dists), color='black',
                         alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_diff.plot(range(bs.size), bs / np.sum(bs) - def_dists / np.sum(def_dists), color='black',
                             alpha=0.2)
        ax_raw_diff.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_diff.set_xlabel('Protons in Bin')
        ax_raw_diff.legend()
        fig_raw_diff.tight_layout()
        fig_raw_diff.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                     f'Protons, {div} divs')

        fig_raw_div, ax_raw_div = plt.subplots()
        ax_raw_div.axhline(1, ls='--', color=color)
        ax_raw_div.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) /
                        (def_dists / np.sum(def_dists)), color='black',
                        alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_div.plot(range(bs.size), bs / np.sum(bs) / (def_dists / np.sum(def_dists)), color='black',
                            alpha=0.2)
        ax_raw_div.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_div.set_xlabel('Protons in Bin')
        ax_raw_div.legend()
        fig_raw_div.tight_layout()
        fig_raw_div.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    raw_sub_mix_def = raw_tp_dist / np.sum(raw_tp_dist) - mix_tp_dist / np.sum(mix_tp_dist)
    fig_raw_sub_mix, ax_raw_sub_mix = plt.subplots()
    ax_raw_sub_mix.axhline(0, ls='--', color='black')
    ax_raw_sub_mix.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='Raw - Mix', zorder=4)
    raw_sub_mix_bss = [raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_tp_dists for mix in mix_bs_tp_dists]
    raw_sub_mix_bss = np.array(raw_sub_mix_bss)
    bs_label = False
    for raw_sub_mix_bs in raw_sub_mix_bss:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix.plot(range(raw_sub_mix_bs.size), raw_sub_mix_bs, alpha=0.3, color='black')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    # ax_raw_sub_mix.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist))
    ax_raw_sub_mix.legend()
    ax_raw_sub_mix.set_title(f'Raw - Mix {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_sub_mix.tight_layout()
    fig_raw_sub_mix.canvas.manager.set_window_title(f'Raw - Mix {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    fig_raw_sub_mix_sd, ax_raw_sub_mix_sd = plt.subplots()
    ax_raw_sub_mix_sd.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sd.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='Raw - Mix')
    raw_sub_mix_sd = np.std(raw_sub_mix_bss, axis=0)
    # n_bands = 100
    # for sigma in np.linspace(0.0, 10, n_bands):
    #     ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - sigma * raw_sub_mix_sd,
    #                                    raw_sub_mix_def + sigma * raw_sub_mix_sd, color='red', alpha=1.0 / n_bands)
    ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - raw_sub_mix_sd,
                                   raw_sub_mix_def + raw_sub_mix_sd, color='red', alpha=0.4)
    # ax_raw_sub_mix_sd.fill_between(range(raw_tp_dist.size), raw_sub_mix_def - 2 * raw_sub_mix_sd,
    #                                raw_sub_mix_def + 2 * raw_sub_mix_sd, color='red', alpha=0.4)
    ax_raw_sub_mix_sd.legend()
    ax_raw_sub_mix_sd.set_title(f'Raw - Mix SD {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_sub_mix_sd.tight_layout()
    fig_raw_sub_mix_sd.canvas.manager.set_window_title(f'Raw - Mix SD {name.capitalize()} {total_proton} '
                                                       f'Protons, {div} divs')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = sim_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist_sim, raw_bs_tp_dists_sim = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist_sim, mix_bs_tp_dists_sim = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    raw_sub_mix_def_sim = raw_tp_dist_sim / np.sum(raw_tp_dist_sim) - mix_tp_dist_sim / np.sum(mix_tp_dist_sim)
    raw_sub_mix_bss_sim = [raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_tp_dists_sim
                           for mix in mix_bs_tp_dists_sim]
    raw_sub_mix_sd_sim = np.std(np.array(raw_sub_mix_bss_sim), axis=0)

    fig_raw_sub_mix_sim_comp, ax_raw_sub_mix_sim_comp = plt.subplots()
    ax_raw_sub_mix_sim_comp.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='white', lw=5, zorder=4)
    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def.size), raw_sub_mix_def, color='red', label='BES Raw - Mix',
                                 zorder=4)
    bs_label = False
    for raw_sub_mix_bs in raw_sub_mix_bss:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_bs.size), raw_sub_mix_bs, alpha=0.3, color='red')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True

    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def_sim.size), raw_sub_mix_def_sim, color='white', lw=5, zorder=4)
    ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_def_sim.size), raw_sub_mix_def_sim, color='purple',
                                 label='Sim Raw - Mix', zorder=4)
    bs_label = False
    for raw_sub_mix_bs in raw_sub_mix_bss_sim:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix_sim_comp.plot(range(raw_sub_mix_bs.size), raw_sub_mix_bs, alpha=0.3, color='purple')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_sub_mix_sim_comp.legend()
    ax_raw_sub_mix_sim_comp.set_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_sub_mix_sim_comp.tight_layout()
    fig_raw_sub_mix_sim_comp.canvas.manager.set_window_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    # Plot difference between bes and sim
    diff_def = raw_sub_mix_def - raw_sub_mix_def_sim
    fig_raw_sub_mix_sim_diff, ax_raw_sub_mix_sim_diff = plt.subplots()
    ax_raw_sub_mix_sim_diff.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_diff.plot(range(diff_def.size), diff_def, color='white', lw=5, zorder=4)
    ax_raw_sub_mix_sim_diff.plot(range(diff_def.size), diff_def, color='magenta', marker='o', label='BES - Sim',
                                 zorder=4)
    for bes, sim in zip(raw_sub_mix_bss, raw_sub_mix_bss_sim):
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix_sim_diff.plot(range(bes.size), bes - sim, alpha=0.3, color='magenta')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_sub_mix_sim_diff.legend()
    ax_raw_sub_mix_sim_diff.set_title(f'Raw - Mix BES Sim Diff {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_sub_mix_sim_diff.tight_layout()
    fig_raw_sub_mix_sim_diff.canvas.manager.set_window_title(f'Raw - Mix BES Sim Diff {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    diff2_def = (raw_sub_mix_def - raw_sub_mix_def_sim) ** 2
    fig_raw_sub_mix_sim_diff2, ax_raw_sub_mix_sim_diff2 = plt.subplots()
    ax_raw_sub_mix_sim_diff2.axhline(0, ls='--', color='black')
    ax_raw_sub_mix_sim_diff2.plot(range(diff2_def.size), diff2_def, color='white', lw=5, zorder=4)
    ax_raw_sub_mix_sim_diff2.plot(range(diff2_def.size), diff2_def, color='magenta', marker='o', label='BES - Sim',
                                  zorder=4)
    for bes, sim in zip(raw_sub_mix_bss, raw_sub_mix_bss_sim):
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_sub_mix_sim_diff2.plot(range(bes.size), (bes - sim) ** 2, alpha=0.3, color='magenta')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_sub_mix_sim_diff2.legend()
    ax_raw_sub_mix_sim_diff2.set_title(f'Raw - Mix BES Sim Diff^2 {name.capitalize()} {total_proton} '
                                       f'Protons, {div} divs')
    fig_raw_sub_mix_sim_diff2.tight_layout()
    fig_raw_sub_mix_sim_diff2.canvas.manager.set_window_title(f'Raw - Mix BES Sim Diff^2 {name.capitalize()} '
                                                              f'{total_proton} Protons, {div} divs')

    # chi2 = (raw_sub_mix_def - raw_sub_mix_def_sim) ** 2 / (raw_sub_mix_sd ** 2 + raw_sub_mix_sd_sim ** 2)
    # fig_raw_sub_mix_sim_chi, ax_raw_sub_mix_sim_chi = plt.subplots()
    # ax_raw_sub_mix_sim_chi.axhline(0, ls='--', color='black')
    # ax_raw_sub_mix_sim_chi.plot(range(chi2.size), chi2, color='salmon', marker='o', label='BES - Sim Chi2')
    # ax_raw_sub_mix_sim_chi.legend()
    # ax_raw_sub_mix_sim_chi.set_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} {total_proton} '
    #                                  f'Protons, {div} divs')
    # fig_raw_sub_mix_sim_chi.tight_layout()
    # fig_raw_sub_mix_sim_chi.canvas.manager.set_window_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} '
    #                                                         f'{total_proton} Protons, {div} divs')

    fig_weights, ax_weights = plt.subplots()
    plot_set = (('bes_raw', 'blue', raw_tp_dist), ('bes_mix', 'green', mix_tp_dist),
                ('sim_raw', 'cyan', raw_tp_dist_sim), ('sim_mix', 'lime', mix_tp_dist_sim))
    for name, color, def_dist in plot_set:
        ax_weights.plot(range(def_dist.size), def_dist / np.sum(def_dist), color=color, label=name, alpha=0.8)
    weights = np.sum([x[2] for x in plot_set], axis=0) / len(plot_set)
    ax_weights.plot(range(weights.size), weights, color='red', label='Average')
    ax_weights.legend()
    ax_weights.set_title(f'Raw & Mix BES & Sim {name.capitalize()} {total_proton} '
                         f'Protons, {div} divs')
    fig_weights.tight_layout()
    fig_weights.canvas.manager.set_window_title(f'Raw & Mix BES & Sim {name.capitalize()} {total_proton} '
                                                f'Protons, {div} divs')

    plt.show()


def full_bs_div_visual():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 60

    total_proton = 21

    rmm_bs_plot_frac = 0.001

    data_set = (base_path, 'default_resample', 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0',
                energy, cent, div, total_proton, 'Data', 'Data_Mix', 'bes')
    # sim_par = ('0125', '08')
    sim_par = ('0', '1')
    sim_set = (base_path, f'flat80_anticlmulti_spread{sim_par[1]}_amp{sim_par[0]}_resample',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
               62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
               f'Sim_spread{sim_par[1]}_amp{sim_par[0]}')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = data_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist, raw_bs_tp_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist, mix_bs_tp_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    for rm_name, color, def_dists, bs_dists in \
            (('raw', 'blue', raw_tp_dist, raw_bs_tp_dists), ('mix', 'green', mix_tp_dist, mix_bs_tp_dists)):
        fig_raw, ax_raw = plt.subplots()
        ax_raw.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]), color='black', alpha=0.2,
                    label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists[1:]:
            ax_raw.plot(range(bs.size), bs / np.sum(bs), color='black', alpha=0.2)
        ax_raw.plot(range(def_dists.size), def_dists / np.sum(def_dists), alpha=1, color=color,
                    label=rm_name.capitalize())
        ax_raw.set_title(f'{rm_name.capitalize()} {name.capitalize()} {energy}GeV, {total_proton} Protons, '
                         f'{div} divs')
        ax_raw.axhline(0, ls='--', color='black', zorder=0)
        ax_raw.set_xlabel('Protons in Bin')
        ax_raw.legend()
        fig_raw.tight_layout()
        fig_raw.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} '
                                                f'{total_proton} Protons, {div} divs')

        fig_raw_diff, ax_raw_diff = plt.subplots()
        ax_raw_diff.axhline(0, ls='--', color=color)
        ax_raw_diff.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) -
                         def_dists / np.sum(def_dists), color='black',
                         alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_diff.plot(range(bs.size), bs / np.sum(bs) - def_dists / np.sum(def_dists), color='black',
                             alpha=0.2)
        ax_raw_diff.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_diff.set_xlabel('Protons in Bin')
        ax_raw_diff.legend()
        fig_raw_diff.tight_layout()
        fig_raw_diff.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                     f'Protons, {div} divs')

        fig_raw_div, ax_raw_div = plt.subplots()
        ax_raw_div.axhline(1, ls='--', color=color)
        ax_raw_div.plot(range(bs_dists[0].size), bs_dists[0] / np.sum(bs_dists[0]) /
                        (def_dists / np.sum(def_dists)), color='black',
                        alpha=0.2, label=f'{len(bs_dists)} Bootstraps')
        for bs in bs_dists:
            ax_raw_div.plot(range(bs.size), bs / np.sum(bs) / (def_dists / np.sum(def_dists)), color='black',
                            alpha=0.2)
        ax_raw_div.set_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} Protons, {div} divs')
        ax_raw_div.set_xlabel('Protons in Bin')
        ax_raw_div.legend()
        fig_raw_div.tight_layout()
        fig_raw_div.canvas.manager.set_window_title(f'{rm_name.capitalize()} {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    raw_div_mix_def = (raw_tp_dist / np.sum(raw_tp_dist)) / (mix_tp_dist / np.sum(mix_tp_dist))
    fig_raw_div_mix, ax_raw_div_mix = plt.subplots()
    ax_raw_div_mix.axhline(0, ls='--', color='black')
    ax_raw_div_mix.plot(range(raw_div_mix_def.size), raw_div_mix_def, color='red', label='Raw - Mix', zorder=4)
    raw_div_mix_bss = [(raw / np.sum(raw)) / (mix / np.sum(mix)) for raw in raw_bs_tp_dists for mix in mix_bs_tp_dists]
    raw_div_mix_bss = np.array(raw_div_mix_bss)
    bs_label = False
    for raw_div_mix_bs in raw_div_mix_bss:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_div_mix.plot(range(raw_div_mix_bs.size), raw_div_mix_bs, alpha=0.3, color='black')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    # ax_raw_div_mix.plot(range(raw_tp_dist.size), raw_tp_dist / np.sum(raw_tp_dist))
    ax_raw_div_mix.legend()
    ax_raw_div_mix.set_title(f'Raw - Mix {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_div_mix.tight_layout()
    fig_raw_div_mix.canvas.manager.set_window_title(f'Raw - Mix {name.capitalize()} {total_proton} '
                                                    f'Protons, {div} divs')

    fig_raw_div_mix_sd, ax_raw_div_mix_sd = plt.subplots()
    ax_raw_div_mix_sd.axhline(0, ls='--', color='black')
    ax_raw_div_mix_sd.plot(range(raw_div_mix_def.size), raw_div_mix_def, color='red', label='Raw - Mix')
    raw_div_mix_sd = np.std(raw_div_mix_bss, axis=0)
    # n_bands = 100
    # for sigma in np.linspace(0.0, 10, n_bands):
    #     ax_raw_div_mix_sd.fill_between(range(raw_tp_dist.size), raw_div_mix_def - sigma * raw_div_mix_sd,
    #                                    raw_div_mix_def + sigma * raw_div_mix_sd, color='red', alpha=1.0 / n_bands)
    ax_raw_div_mix_sd.fill_between(range(raw_tp_dist.size), raw_div_mix_def - raw_div_mix_sd,
                                   raw_div_mix_def + raw_div_mix_sd, color='red', alpha=0.4)
    # ax_raw_div_mix_sd.fill_between(range(raw_tp_dist.size), raw_div_mix_def - 2 * raw_div_mix_sd,
    #                                raw_div_mix_def + 2 * raw_div_mix_sd, color='red', alpha=0.4)
    ax_raw_div_mix_sd.legend()
    ax_raw_div_mix_sd.set_title(f'Raw / Mix SD {name.capitalize()} {total_proton} Protons, {div} divs')
    fig_raw_div_mix_sd.tight_layout()
    fig_raw_div_mix_sd.canvas.manager.set_window_title(f'Raw / Mix SD {name.capitalize()} {total_proton} '
                                                       f'Protons, {div} divs')

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = sim_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_tp_dist_sim, raw_bs_tp_dists_sim = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_tp_dist_sim, mix_bs_tp_dists_sim = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    raw_div_mix_def_sim = raw_tp_dist_sim / np.sum(raw_tp_dist_sim) - mix_tp_dist_sim / np.sum(mix_tp_dist_sim)
    raw_div_mix_bss_sim = [raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_tp_dists_sim
                           for mix in mix_bs_tp_dists_sim]
    raw_div_mix_sd_sim = np.std(np.array(raw_div_mix_bss_sim), axis=0)

    fig_raw_div_mix_sim_comp, ax_raw_div_mix_sim_comp = plt.subplots()
    ax_raw_div_mix_sim_comp.axhline(0, ls='--', color='black')
    ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_def.size), raw_div_mix_def, color='white', lw=5, zorder=4)
    ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_def.size), raw_div_mix_def, color='red', label='BES Raw - Mix',
                                 zorder=4)
    bs_label = False
    for raw_div_mix_bs in raw_div_mix_bss:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_bs.size), raw_div_mix_bs, alpha=0.3, color='red')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True

    ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_def_sim.size), raw_div_mix_def_sim, color='white', lw=5, zorder=4)
    ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_def_sim.size), raw_div_mix_def_sim, color='purple',
                                 label='Sim Raw - Mix', zorder=4)
    bs_label = False
    for raw_div_mix_bs in raw_div_mix_bss_sim:
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_div_mix_sim_comp.plot(range(raw_div_mix_bs.size), raw_div_mix_bs, alpha=0.3, color='purple')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_div_mix_sim_comp.legend()
    ax_raw_div_mix_sim_comp.set_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_div_mix_sim_comp.tight_layout()
    fig_raw_div_mix_sim_comp.canvas.manager.set_window_title(f'Raw - Mix BES Sim Comparison {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    # Plot difference between bes and sim
    diff_def = raw_div_mix_def - raw_div_mix_def_sim
    fig_raw_div_mix_sim_diff, ax_raw_div_mix_sim_diff = plt.subplots()
    ax_raw_div_mix_sim_diff.axhline(0, ls='--', color='black')
    ax_raw_div_mix_sim_diff.plot(range(diff_def.size), diff_def, color='white', lw=5, zorder=4)
    ax_raw_div_mix_sim_diff.plot(range(diff_def.size), diff_def, color='magenta', marker='o', label='BES - Sim',
                                 zorder=4)
    for bes, sim in zip(raw_div_mix_bss, raw_div_mix_bss_sim):
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_div_mix_sim_diff.plot(range(bes.size), bes - sim, alpha=0.3, color='magenta')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_div_mix_sim_diff.legend()
    ax_raw_div_mix_sim_diff.set_title(f'Raw - Mix BES Sim Diff {name.capitalize()} {total_proton} '
                                      f'Protons, {div} divs')
    fig_raw_div_mix_sim_diff.tight_layout()
    fig_raw_div_mix_sim_diff.canvas.manager.set_window_title(f'Raw - Mix BES Sim Diff {name.capitalize()} '
                                                             f'{total_proton} Protons, {div} divs')

    diff2_def = (raw_div_mix_def - raw_div_mix_def_sim) ** 2
    fig_raw_div_mix_sim_diff2, ax_raw_div_mix_sim_diff2 = plt.subplots()
    ax_raw_div_mix_sim_diff2.axhline(0, ls='--', color='black')
    ax_raw_div_mix_sim_diff2.plot(range(diff2_def.size), diff2_def, color='white', lw=5, zorder=4)
    ax_raw_div_mix_sim_diff2.plot(range(diff2_def.size), diff2_def, color='magenta', marker='o', label='BES - Sim',
                                  zorder=4)
    for bes, sim in zip(raw_div_mix_bss, raw_div_mix_bss_sim):
        if np.random.random() > rmm_bs_plot_frac:  # Only plot x percent
            continue
        bs, = ax_raw_div_mix_sim_diff2.plot(range(bes.size), (bes - sim) ** 2, alpha=0.3, color='magenta')
        if not bs_label:
            bs.set_label(f'{rmm_bs_plot_frac} * $250^2$ bootstrap samples')
            bs_label = True
    ax_raw_div_mix_sim_diff2.legend()
    ax_raw_div_mix_sim_diff2.set_title(f'Raw - Mix BES Sim Diff^2 {name.capitalize()} {total_proton} '
                                       f'Protons, {div} divs')
    fig_raw_div_mix_sim_diff2.tight_layout()
    fig_raw_div_mix_sim_diff2.canvas.manager.set_window_title(f'Raw - Mix BES Sim Diff^2 {name.capitalize()} '
                                                              f'{total_proton} Protons, {div} divs')

    # chi2 = (raw_div_mix_def - raw_div_mix_def_sim) ** 2 / (raw_div_mix_sd ** 2 + raw_div_mix_sd_sim ** 2)
    # fig_raw_div_mix_sim_chi, ax_raw_div_mix_sim_chi = plt.subplots()
    # ax_raw_div_mix_sim_chi.axhline(0, ls='--', color='black')
    # ax_raw_div_mix_sim_chi.plot(range(chi2.size), chi2, color='salmon', marker='o', label='BES - Sim Chi2')
    # ax_raw_div_mix_sim_chi.legend()
    # ax_raw_div_mix_sim_chi.set_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} {total_proton} '
    #                                  f'Protons, {div} divs')
    # fig_raw_div_mix_sim_chi.tight_layout()
    # fig_raw_div_mix_sim_chi.canvas.manager.set_window_title(f'Raw - Mix BES Sim Chi2 {name.capitalize()} '
    #                                                         f'{total_proton} Protons, {div} divs')

    fig_weights, ax_weights = plt.subplots()
    plot_set = (('bes_raw', 'blue', raw_tp_dist), ('bes_mix', 'green', mix_tp_dist),
                ('sim_raw', 'cyan', raw_tp_dist_sim), ('sim_mix', 'lime', mix_tp_dist_sim))
    for name, color, def_dist in plot_set:
        ax_weights.plot(range(def_dist.size), def_dist / np.sum(def_dist), color=color, label=name, alpha=0.8)
    weights = np.sum([x[2] for x in plot_set], axis=0) / len(plot_set)
    ax_weights.plot(range(weights.size), weights, color='red', label='Average')
    ax_weights.legend()
    ax_weights.set_title(f'Raw & Mix BES & Sim {name.capitalize()} {total_proton} '
                         f'Protons, {div} divs')
    fig_weights.tight_layout()
    fig_weights.canvas.manager.set_window_title(f'Raw & Mix BES & Sim {name.capitalize()} {total_proton} '
                                                f'Protons, {div} divs')

    plt.show()


if __name__ == '__main__':
    main()
