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
    sim_self_compare()
    print('donzo')


def sim_self_compare():
    base_path = 'D:/Research/'
    energy = 62
    cent = 8
    div = 120

    total_proton = 21

    edge_buffer = 4  # Number of points to skip on either extreme to avoid boundaries
    rmm_bs_plot_frac = 0

    # comp_sim_par = ('15', '1')
    # comp_sim_set = (base_path, f'flat80_anticlmulti_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}_resample',
    #                 f'Sim_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}_flat80_anticlmulti_norotate_resample_0',
    #                 62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix',
    #                 f'Sim_spread{comp_sim_par[1]}_amp{comp_sim_par[0]}')

    sim_pars = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'], True)
    sim_amps = pd.unique(sim_pars['amp'])
    print(sim_amps)

    sim_sets = {}
    # sim_names, comp_names, amp_fs, amp_fs_comp = [], [], [], []
    sim_names, amp_fs = [], []
    spread = '1'
    for amp in sim_amps:
        name = f'Sim_spread{spread}_amp{amp}'
        amp_f, spread_f = get_name_amp_spread(name)
        sim_pars = (base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                    f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                    62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix', name, amp_f, spread_f)
        sim_sets.update({name: calc_set(sim_pars)})
        amp_fs.append(amp_f)
        sim_names.append(name)
        # Filter any unwanted comparison sets here
        # amp_fs_comp.append(amp_f)
        # comp_names.append(name)

    amp_fs, sim_names = zip(*sorted(zip(amp_fs, sim_names)))
    # amp_fs_comp, comp_names = zip(*sorted(zip(amp_fs_comp, comp_names)))
    amp_fs_comp, comp_names = amp_fs[edge_buffer:-edge_buffer], sim_names[edge_buffer:-edge_buffer]  # skip edges

    # comp_amp, comp_cost_def, comp_cost_sd = [], [], []
    costs_colors = {'Diff2': 'blue', 'Chi2': 'green', 'Diff2_Div': 'red', 'Chi2_Div': 'purple'}
    comp_amp = []
    comp_cost_min_def, comp_cost_min_sd = [{cost: [] for cost in costs_colors.keys()} for i in range(2)]
    for comp_name in comp_names:
        print(f'Comparison Set: {comp_name}')
        raw_sub_mix_comp, raw_sub_mix_sd_comp, raw_sub_mix_bss_slim_comp, raw_div_mix_comp, raw_div_mix_sd_comp, \
            raw_div_mix_bss_slim_comp = sim_sets[comp_name]

        true_amp = get_name_amp_spread(comp_name)[0]

        amps, diff2s, chi2s = [], [], []
        plot_comp = true_amp == 0.1
        costs = {cost: [] for cost in costs_colors.keys()}
        costs_bss = {cost: [] for cost in costs_colors.keys()}

        for sim_name in sim_names:
            if sim_name == comp_name:
                continue  # Don't compare to self.

            raw_sub_mix, raw_sub_mix_sd, raw_sub_mix_bss_slim, raw_div_mix, raw_div_mix_sd, raw_div_mix_bss_slim = \
                sim_sets[sim_name]
            amp = get_name_amp_spread(sim_name)[0]

            diff = raw_sub_mix_comp - raw_sub_mix
            diff2 = diff ** 2
            err2 = (raw_sub_mix_sd_comp + raw_sub_mix_sd) ** 2
            chi2 = diff2 / err2

            diff_div = raw_div_mix_comp - raw_div_mix
            diff2_div = diff_div ** 2
            err2_div = (raw_div_mix_sd_comp + raw_div_mix_sd) ** 2
            chi2_div = diff2_div / err2_div

            diff2_bss, chi2_bss, diff2_div_bss, chi2_div_bss = [], [], [], []
            for i in range(len(raw_sub_mix_bss_slim_comp)):
                diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
                diff2_bs = diff_bs ** 2
                diff2_bss.append(np.sum(diff2_bs))

                chi2_bs = diff2_bs / err2  # Placeholder probably?
                chi2_bss.append(np.sum(chi2_bs))

                diff_div_bs = raw_div_mix_bss_slim_comp[i] - raw_div_mix_bss_slim[i]
                diff2_div_bs = diff_div_bs ** 2
                diff2_div_bss.append(np.nansum(diff2_div_bs))

                chi2_div_bs = diff2_div_bs / err2_div  # Placeholder probably?
                chi2_div_bss.append(np.nansum(chi2_div_bs))

            diff2 = np.sum(diff2)
            chi2 = np.sum(chi2)
            diff2_div = np.nansum(diff2_div)
            chi2_div = np.nansum(chi2_div)

            amps.append(amp)
            costs['Diff2'].append(diff2)
            costs_bss['Diff2'].append(diff2_bss)
            costs['Chi2'].append(chi2)
            costs_bss['Chi2'].append(chi2_bss)
            costs['Diff2_Div'].append(diff2_div)
            costs_bss['Diff2_Div'].append(diff2_div_bss)
            costs['Chi2_Div'].append(chi2_div)
            costs_bss['Chi2_Div'].append(chi2_div_bss)

        costs_norm, costs_interp, costs_min_res = [{key: [] for key in costs} for i in range(3)]
        # costs_min, costs_sd = [{key: [] for key in costs} for i in range(2)]
        for cost in costs.keys():
            costs_bss[cost] = list(np.array(costs_bss[cost]).T)  # Group by amps instead of bs #
            for cost_i in [costs[cost], *costs_bss[cost]]:
                cost_norm = (np.array(cost_i) - min(cost_i)) / (max(cost_i) - min(cost_i))
                cost_interp = interp1d(amps, cost_norm, 'cubic', bounds_error=True, fill_value=max(cost_norm))
                # print(cost_i)
                # cost_interp = interp1d(amps, cost_i, 'cubic', bounds_error=True, fill_value=max(cost_i))
                cost_min_res = minimize_scalar(cost_interp, method='bounded', bounds=(min(amps), max(amps)))

                costs_norm[cost].append(cost_norm)
                costs_interp[cost].append(cost_interp)
                costs_min_res[cost].append(cost_min_res)
            # costs_min[cost] = costs_min_res[cost][0]
            # costs_sd[cost] = np.std([min_res.x for min_res in costs_min_res[cost][1:]])
            comp_cost_min_def[cost].append(costs_min_res[cost][0].x)
            comp_cost_min_sd[cost].append(np.std([min_res.x for min_res in costs_min_res[cost][1:]]))

        # comp_cost_def.append(costs_min)
        # comp_cost_sd.append(costs_sd)
        comp_amp.append(true_amp)

        if plot_comp:
            plot_comp_set(amps, costs_norm, costs_interp, costs_min_res, costs_colors, true_amp)

    fig_min_comp_amp, ax_min_comp_amp = plt.subplots()
    ax_min_comp_amp.axhline(0, ls='--', color='gray', label='True')
    comp_amp = np.array(comp_amp)
    num_costs = len(comp_cost_min_def)
    space = np.min(np.diff(comp_amp)) / num_costs / 3
    for i, cost in enumerate(comp_cost_min_def.keys()):
        dev_from_true = np.array(comp_cost_min_def[cost]) - comp_amp
        ax_min_comp_amp.errorbar(comp_amp - space * (i - num_costs / 2 + 0.5), dev_from_true, color=costs_colors[cost],
                                 yerr=comp_cost_min_sd[cost], ls='none', marker='o', label=cost)
    ax_min_comp_amp.set_xlabel('Comparison Amplitude')
    ax_min_comp_amp.set_ylabel('Minimum Value of Cost Interpolation Minus True Value')
    ax_min_comp_amp.set_title('Cost Minima vs Comparison Set Amplitude')
    ax_min_comp_amp.legend()
    fig_min_comp_amp.canvas.manager.set_window_title('Cost Minima vs Comparison Set Amplitude')

    fig_diff_sd_cor, ax_diff_sd_cor = plt.subplots()
    for cost in comp_cost_min_def.keys():
        dev_from_true = np.array(comp_cost_min_def[cost]) - comp_amp
        ax_diff_sd_cor.scatter(abs(dev_from_true), comp_cost_min_sd[cost], color=costs_colors[cost], label=cost)
    ax_diff_sd_cor.set_xlabel('|Deviation| from True Value')
    ax_diff_sd_cor.set_ylabel('Minimum Uncertainty')
    ax_diff_sd_cor.set_title('Correlation Uncertainty vs |Deviation| from True')
    ax_diff_sd_cor.legend()
    fig_diff_sd_cor.tight_layout()
    fig_diff_sd_cor.canvas.manager.set_window_title('Correlation Uncertainty vs |Deviation| from True')

    plt.show()


def calc_set(set_pars):
    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name, amp, spread \
        = set_pars
    print(f'Calculate Set: {name}')
    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_dist, raw_bs_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_dist, mix_bs_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)
    raw_sum, mix_sum = np.sum(raw_dist), np.sum(mix_dist)
    raw_bs_sums, mix_bs_sums = [np.sum(raw) for raw in raw_bs_dists], [np.sum(mix) for mix in mix_bs_dists]
    raw_sub_mix = raw_dist / raw_sum - mix_dist / mix_sum
    raw_sub_mix_bss = np.array([raw / raw_sum_bs - mix / mix_sum_bs for raw, raw_sum_bs in
                                zip(raw_bs_dists, raw_bs_sums) for mix, mix_sum_bs in zip(mix_bs_dists, mix_bs_sums)])
    raw_sub_mix_sd = np.std(raw_sub_mix_bss)
    raw_sub_mix_bss_slim = np.array([raw / raw_sum_bs - mix / mix_sum_bs for raw, raw_sum_bs, mix, mix_sum_bs in
                                     zip(raw_bs_dists, raw_bs_sums, mix_bs_dists, mix_bs_sums)])
    raw_div_mix = (raw_dist / raw_sum) / (mix_dist / mix_sum)
    raw_div_mix_bss = np.array([(raw / raw_sum_bs) / (mix / mix_sum_bs) for raw, raw_sum_bs in
                                zip(raw_bs_dists, raw_bs_sums) for mix, mix_sum_bs in zip(mix_bs_dists, mix_bs_sums)])
    raw_div_mix_sd = np.std(raw_div_mix_bss)
    raw_div_mix_bss_slim = np.array([(raw / raw_sum_bs) / (mix / mix_sum_bs) for raw, raw_sum_bs, mix, mix_sum_bs in
                                     zip(raw_bs_dists, raw_bs_sums, mix_bs_dists, mix_bs_sums)])

    # set_dists = {'raw_sub_mix': raw_sub_mix, 'raw_sub_mix_sd': raw_sub_mix_sd,
    #              'raw_sub_mix_bss_slim': raw_sub_mix_bss_slim}
    set_dists = [raw_sub_mix, raw_sub_mix_sd, raw_sub_mix_bss_slim, raw_div_mix, raw_div_mix_sd, raw_div_mix_bss_slim]

    return set_dists


def plot_comp_set(amps, costs_norm, costs_interp, costs_min_res, costs_colors, true_amp, show=False):
    xs = np.linspace(min(amps), max(amps), 10000)

    fig_cost_amp, ax_cost_amp = plt.subplots()
    ax_cost_amp.axvline(true_amp, color='black', ls=':', label='True')
    ax_cost_amp.axhline(0, color='black', ls='--')
    for cost in costs_norm.keys():
        ax_cost_amp.plot(xs, costs_interp[cost][0](xs), color=costs_colors[cost])
        ax_cost_amp.scatter(amps, costs_norm[cost][0], color=costs_colors[cost], label=cost)
        ax_cost_amp.axvline(costs_min_res[cost][0].x, color=costs_colors[cost], label=f'{cost} Minimum')
        for i in range(1, len(costs_norm[cost])):
            ax_cost_amp.plot(xs, costs_interp[cost][i](xs), color=costs_colors[cost], alpha=0.2)
            ax_cost_amp.axvline(costs_min_res[cost][i].x, color=costs_colors[cost], alpha=0.2)

    ax_cost_amp.set_title(f'Cost vs Amp')
    ax_cost_amp.set_ylabel('Cost')
    ax_cost_amp.set_xlabel('Amp')
    ax_cost_amp.legend(loc='upper right')
    fig_cost_amp.tight_layout()
    fig_cost_amp.canvas.manager.set_window_title(f'Cost vs Amp')

    fig_cost_min, ax_cost_min = plt.subplots()
    for cost in costs_norm.keys():
        bs_minima = [min_res.x for min_res in costs_min_res[cost][1:]]
        ax_cost_min.hist(bs_minima, histtype='step', color=costs_colors[cost], label=f'{cost} bootstraps')
        ax_cost_min.axvline(costs_min_res[cost][0].x, ls='--', color=costs_colors[cost], label=f'{cost} default')
        bs_minima_sd = np.std(bs_minima)
        ax_cost_min.axvspan(costs_min_res[cost][0].x - bs_minima_sd, costs_min_res[cost][0].x + bs_minima_sd,
                            color=costs_colors[cost], alpha=0.4, label=f'{cost} BS 1σ')
    ax_cost_min.axvline(true_amp, ls='--', color='black', label='True')
    ax_cost_min.set_xlabel('Amp')
    ax_cost_min.set_title('Cost Amp Minima Bootstrap distributions')
    ax_cost_min.legend()
    fig_cost_min.tight_layout()
    fig_cost_min.canvas.manager.set_window_title(f'Cost Amp Minima')

    if show:
        plt.show()


if __name__ == '__main__':
    main()
