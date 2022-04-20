#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 12 2:08 PM 2022
Created in PyCharm
Created as QGP_Scripts/sim_sim_comp.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, leastsq
from scipy.optimize import minimize

from compare_dists import get_norm_dists
from calc_binom_slices import find_sim_sets, get_name_amp_spread


def main():
    test()
    print('donzo')


def test():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 60

    total_proton = 21

    rmm_bs_plot_frac = 0

    comp_sim_par = ('15', '1')
    comp_sim_set = (base_path, f'flat80_anticlmulti_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}_resample',
                    f'Sim_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
                    62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
                    f'Sim_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}')

    sim_pars = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'], True)
    sim_amps = pd.unique(sim_pars['amp'])
    print(sim_amps)

    sim_sets = []
    spread = '1'
    for amp in sim_amps:
        if amp == comp_sim_par[0] and spread == comp_sim_par[1]:
            continue
        name = f'Sim_spread{spread}_amp{amp}'
        amp, spread = get_name_amp_spread(name)
        sim_sets.append((base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                         f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                         62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
                         name, amp, spread))

    sim_sets = sorted(sim_sets, key = lambda x: x[-2])  # Sort on amp value

    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = comp_sim_set
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_dist_comp, raw_bs_dists_comp = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_dist_comp, mix_bs_dists_comp = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)
    raw_sub_mix_comp = raw_dist_comp / np.sum(raw_dist_comp) - mix_dist_comp / np.sum(mix_dist_comp)
    raw_sub_mix_bss_comp = np.array([raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_dists_comp
                                     for mix in mix_bs_dists_comp])
    raw_sub_mix_sd_comp = np.std(raw_sub_mix_bss_comp)
    raw_sub_mix_bss_slim_comp = np.array([raw / np.sum(raw) - mix / np.sum(mix) for raw, mix in
                                          zip(raw_bs_dists_comp, mix_bs_dists_comp)])

    amps, diff2s, chi2s = [], [], []
    costs_colors = {'Diff2': 'blue', 'Chi2': 'green'}
    costs = {cost: [] for cost in costs_colors.keys()}
    costs_bss = {cost: [] for cost in costs_colors.keys()}
    # i = 0
    for sim_set in sim_sets:
        # i += 1
        # if i > 5:
        #     break
        base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name = sim_set
        print(name)
        file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
        path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
        raw_dist, raw_bs_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
        mix_dist, mix_bs_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)
        raw_sub_mix = raw_dist / np.sum(raw_dist) - mix_dist / np.sum(mix_dist)
        raw_sub_mix_bss = np.array([raw / np.sum(raw) - mix / np.sum(mix) for raw in raw_bs_dists
                                    for mix in mix_bs_dists])
        raw_sub_mix_sd = np.std(raw_sub_mix_bss)
        raw_sub_mix_bss_slim = np.array([raw / np.sum(raw) - mix / np.sum(mix) for raw, mix in
                                         zip(raw_bs_dists, mix_bs_dists)])

        amp, spread = get_name_amp_spread(name)

        diff = raw_sub_mix_comp - raw_sub_mix
        diff2 = diff ** 2
        err2 = (raw_sub_mix_sd_comp + raw_sub_mix_sd) ** 2
        chi2 = diff2 / err2

        diff2_bss = []
        chi2_bss = []
        for i in range(len(raw_sub_mix_bss_slim_comp)):
            diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
            diff2_bs = diff_bs**2
            diff2_bss.append(np.sum(diff2_bs))
            chi2_bs = diff2_bs / err2  # Placeholder probably?
            chi2_bss.append(np.sum(chi2_bs))

        diff2 = np.sum(diff2)
        chi2 = np.sum(chi2)

        # fig_diff, ax_diff = plt.subplots()
        # ax_diff.axhline(0, color='black', ls='--')
        # ax_diff.set_title(f'Diff Amp {amp}')
        # ax_diff.plot(range(diff.size), diff)
        # fig_diff.tight_layout()
        # fig_diff.canvas.manager.set_window_title(f'Diff Amp {amp}')

        # fig_diff2, ax_diff2 = plt.subplots()
        # ax_diff2.axhline(0, color='black', ls='--')
        # ax_diff2.set_title(f'Diff2 Amp {amp}')
        # ax_diff2.plot(range(diff2.size), diff2)
        # fig_diff2.tight_layout()
        # fig_diff2.canvas.manager.set_window_title(f'Diff2 Amp {amp}')

        amps.append(amp)
        costs['Diff2'].append(np.sum(diff2))
        costs_bss['Diff2'].append(diff2_bss)
        costs['Chi2'].append(np.sum(chi2))
        costs_bss['Chi2'].append()

    # amps = [-amp for amp in amps] + amps  # Mirror
    # diff2s = diff2s + diff2s
    # chi2s = chi2s + chi2s
    # for i in range(len(amps)):  # Remove duplicate zero
    #     if amps[i] == 0:
    #         amps.pop(i)
    #         diff2s.pop(i)
    #         chi2s.pop(i)
    #         break
    # amps, *costs_vals = zip(*sorted(zip(amps, *costs.values())))
    # costs = {key: val for key, val in zip(costs.keys(), costs_vals)}

    costs_norm, costs_interp, costs_min_res = [{key: [] for key in costs} for i in range(3)]
    for cost in costs.keys():
        for cost_i in [costs[cost], *costs_bss[cost]]:
            costs_norm[cost].append((np.array(val) - min(val)) / (max(val) - min(val)) for key, val in cost_i.items())

    costs_norm = {key: (np.array(val) - min(val)) / (max(val) - min(val)) for key, val in costs.items()}
    costs_interp = {key: interp1d(amps, val, 'cubic', bounds_error=True, fill_value=max(val)) for key, val in
                    costs_norm.items()}
    costs_min_res = {key: minimize_scalar(val, method='bounded', bounds=(0.001, 0.999)) for key, val in
                     costs_interp.items()}
    for cost, min_res in costs_min_res.items():
        print(f'\n{cost} Minimum: \n{min_res}')
    # amps, diff2s, chi2s = zip(*sorted(zip(amps, diff2s, chi2s)))
    # diff2s_norm = (np.array(diff2s) - min(diff2s)) / (max(diff2s) - min(diff2s))  # Normalize
    # chi2s_norm = (np.array(chi2s) - min(chi2s)) / (max(chi2s) - min(chi2s))
    #
    # diff2s_interp = interp1d(amps, diff2s_norm, 'cubic', bounds_error=True, fill_value=max(diff2s_norm))
    # chi2s_interp = interp1d(amps, chi2s_norm, 'cubic', bounds_error=True, fill_value=max(chi2s_norm))
    #
    # diff2s_min_res = minimize_scalar(diff2s_interp, method='bounded', bounds=(0.001, 0.999))
    # chi2s_min_res = minimize_scalar(chi2s_interp, method='bounded', bounds=(0.001, 0.999))
    #
    # # diff2s_lsq_res = leastsq(diff2s_interp, x0=diff2s_min_res.x)
    # # chi2s_lsq_res = leastsq(chi2s_interp, x0=chi2s_min_res.x)
    #
    # print(f'\nDiff Minimum: \n{diff2s_min_res}')
    # # print(np.sqrt(max(1, abs(diff2s_min_res.fun)) * ftol * diff2s_min_res.hess_inv(np.zeros(1))[0])))
    # print(f'\nChi Minimum: \n{chi2s_min_res}')

    # print(f'Diff Err: {np.sqrt(np.diag(diff2s_lsq_res.cov_x))}')

    true_amp = get_name_amp_spread(f'Sim_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}')[0]

    xs = np.linspace(min(amps), max(amps), 10000)

    fig_cost_amp, ax_cost_amp = plt.subplots()
    ax_cost_amp.axvline(true_amp, color='black', ls=':', label='True')
    ax_cost_amp.axhline(0, color='black', ls='--')
    for cost in costs.keys():
        ax_cost_amp.plot(xs, costs_interp[cost](xs), color=costs_colors[cost])
        ax_cost_amp.scatter(amps, costs_norm[cost], color=costs_colors[cost], label=cost)
        ax_cost_amp.axvline(costs_min_res[cost].x, color=costs_colors[cost], label=f'{cost} Minimum')

    # ax_cost_amp.plot(xs, diff2s_interp(xs), color='blue')
    # ax_cost_amp.plot(xs, chi2s_interp(xs), color='green')
    # ax_cost_amp.scatter(amps, diff2s_norm, color='blue', label='Diff2')
    # ax_cost_amp.scatter(amps, chi2s_norm, color='green', label='Chi2')
    # ax_cost_amp.axvline(diff2s_min_res.x, color='blue', label='Diff2 Minimum')
    # ax_cost_amp.axvline(chi2s_min_res.x, color='green', label='Chi2 Minimum')
    ax_cost_amp.set_title(f'Cost vs Amp')
    ax_cost_amp.set_ylabel('Cost')
    ax_cost_amp.set_xlabel('Amp')
    ax_cost_amp.legend()
    fig_cost_amp.tight_layout()
    fig_cost_amp.canvas.manager.set_window_title(f'Cost vs Amp')

    plt.show()


if __name__ == '__main__':
    main()
