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
from scipy import stats

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from compare_dists import get_norm_dists, get_norm_dists_all, get_total_proton_list
from calc_binom_slices import find_sim_sets, get_name_amp_spread


def main():
    sim_self_compare_single()
    # sim_self_compare()
    print('donzo')


def sim_self_compare_single():
    base_path = 'F:/Research/'
    cent = 8
    div = 120

    total_proton = 21

    edge_buffer = 4  # Number of points to skip on either extreme to avoid boundaries

    sim_pars = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'], True)
    sim_amps = pd.unique(sim_pars['amp'])
    print(sim_amps)

    sim_sets = {}
    sim_names, amp_fs = [], []
    spread = '1'
    for amp in sim_amps:
        name = f'Sim_spread{spread}_amp{amp}'
        amp_f, spread_f = get_name_amp_spread(name)
        sim_pars = (base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                    f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                    62, cent, div, total_proton, 'Data_Sim', 'Data_Sim_Mix', name, amp_f, spread_f)
        sim_sets.update({name: get_set(sim_pars)})
        amp_fs.append(amp_f)
        sim_names.append(name)
        # Filter any unwanted comparison sets here
        # amp_fs_comp.append(amp_f)
        # comp_names.append(name)

    amp_fs, sim_names = zip(*sorted(zip(amp_fs, sim_names)))
    amp_fs_comp, comp_names = amp_fs[edge_buffer:-edge_buffer], sim_names[edge_buffer:-edge_buffer]  # skip edges

    costs = {
        'Diff2': {'func': diff2_cost, 'color': 'blue', 'min_def': [], 'min_sd': []},
        'Chi2': {'func': chi2_cost, 'color': 'green', 'min_def': [], 'min_sd': []},
        'Diff2_Weight': {'func': diff2_cut_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        'Chi2_Weight': {'func': chi2_cut_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
        # 'Diff2_Div': {'func': diff2_div_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        # 'Chi2_Div': {'func': chi2_div_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
        # 'Diff2_Div_CompDiv': {'func': diff2_div_compdiv_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        # 'Chi2_Div_CompDiv': {'func': chi2_div_compdiv_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
    }
    costs_colors = {cost: costs[cost]['color'] for cost in costs}

    comp_amp = []
    for comp_name in comp_names:
        print(f'Comparison Set: {comp_name}')

        true_amp = get_name_amp_spread(comp_name)[0]

        amps = []
        plot_comp = true_amp == 0.1
        costs_def, costs_bss = [{cost: [] for cost in costs.keys()} for i in range(2)]

        for sim_name in sim_names:
            if sim_name == comp_name:
                continue  # Don't compare to self.

            amps.append(get_name_amp_spread(sim_name)[0])

            for cost in costs:
                cost_def, cost_bss = costs[cost]['func'](sim_sets[sim_name], sim_sets[comp_name], False)
                costs_def[cost].append(cost_def)
                costs_bss[cost].append(cost_bss)

        costs_norm, costs_interp, costs_min_res = [{key: [] for key in costs} for i in range(3)]
        for cost in costs:
            costs_bss[cost] = list(np.array(costs_bss[cost]).T)  # Group by amps instead of bs #
            for cost_i in [costs_def[cost], *costs_bss[cost]]:
                cost_norm = (np.array(cost_i) - min(cost_i)) / (max(cost_i) - min(cost_i))
                cost_interp = interp1d(amps, cost_norm, 'cubic', bounds_error=True, fill_value=max(cost_norm))
                # print(cost_i)
                # cost_interp = interp1d(amps, cost_i, 'cubic', bounds_error=True, fill_value=max(cost_i))
                cost_min_res = minimize_scalar(cost_interp, method='bounded', bounds=(min(amps), max(amps)))

                costs_norm[cost].append(cost_norm)
                costs_interp[cost].append(cost_interp)
                costs_min_res[cost].append(cost_min_res)
            costs[cost]['min_def'].append(costs_min_res[cost][0].x)
            costs[cost]['min_sd'].append(np.std([min_res.x for min_res in costs_min_res[cost][1:]]))

        comp_amp.append(true_amp)

        if plot_comp:
            plot_comp_set(amps, costs_norm, costs_interp, costs_min_res, costs_colors, true_amp)

    fig_min_comp_amp, ax_min_comp_amp = plt.subplots()
    ax_min_comp_amp.axhline(0, ls='--', color='gray', label='True')
    comp_amp = np.array(comp_amp)
    num_costs = len(costs)
    space = np.min(np.diff(comp_amp)) / num_costs / 3
    for i, cost in enumerate(costs):
        dev_from_true = np.array(costs[cost]['min_def']) - comp_amp
        ax_min_comp_amp.errorbar(comp_amp - space * (i - num_costs / 2 + 0.5), dev_from_true, color=costs_colors[cost],
                                 yerr=costs[cost]['min_sd'], ls='none', marker='o', label=cost)
    ax_min_comp_amp.set_xlabel('Comparison Amplitude')
    ax_min_comp_amp.set_ylabel('Minimum Value of Cost Interpolation Minus True Value')
    ax_min_comp_amp.set_title('Cost Minima vs Comparison Set Amplitude')
    ax_min_comp_amp.legend()
    fig_min_comp_amp.canvas.manager.set_window_title('Cost Minima vs Comparison Set Amplitude')

    fig_diff_sd_cor, ax_diff_sd_cor = plt.subplots()
    for cost in costs:
        dev_from_true = np.array(costs[cost]['min_def']) - comp_amp
        ax_diff_sd_cor.scatter(abs(dev_from_true), costs[cost]['min_sd'], color=costs_colors[cost], label=cost)
    ax_diff_sd_cor.set_xlabel('|Deviation| from True Value')
    ax_diff_sd_cor.set_ylabel('Minimum Uncertainty')
    ax_diff_sd_cor.set_title('Correlation Uncertainty vs |Deviation| from True')
    ax_diff_sd_cor.legend()
    fig_diff_sd_cor.tight_layout()
    fig_diff_sd_cor.canvas.manager.set_window_title('Correlation Uncertainty vs |Deviation| from True')

    for cost in costs:
        sigma_deviation = (np.array(costs[cost]['min_def']) - comp_amp) / np.array(costs[cost]['min_sd'])
        fig_sigmas, ax_sigmas = plt.subplots()
        ax_sigmas.set_xlabel('Sigmas from True')
        ax_sigmas.set_title(f'{cost} Number of Sigmas from True Value')
        ax_sigmas.hist(sigma_deviation, density=True, label=cost)
        norm_xs = np.linspace(*ax_sigmas.get_xlim(), 1000)
        norm_dist = stats.norm()
        ax_sigmas.plot(norm_xs, norm_dist.pdf(norm_xs), color='red', label='Standard Normal')
        ax_sigmas.legend()
        fig_sigmas.tight_layout()
        fig_sigmas.canvas.manager.set_window_title(f'{cost} Number of Sigmas from True Value')

    fig_dev_sum, ax_dev_sum = plt.subplots()
    ax_dev_sum.set_xlabel('Cost Function')
    ax_dev_sum.set_ylabel('RMS Deviation from True')
    cost_labels, dev_true_rmss = [], []
    for cost in costs:
        dev_true_rms = np.sqrt(np.mean((np.array(costs[cost]['min_def']) - comp_amp) ** 2))
        cost_labels.append(cost)
        dev_true_rmss.append(dev_true_rms)
    ax_dev_sum.axhline(0, ls='--', color='black')
    ax_dev_sum.scatter(cost_labels, dev_true_rmss)
    ax_dev_sum.grid()
    fig_dev_sum.tight_layout()
    fig_dev_sum.canvas.manager.set_window_title(f'{cost} RMS Deviation from True')

    plt.show()


def sim_self_compare():
    base_path = 'F:/Research/'
    cent = 8
    # [60, 72, 89, 90, 120, 180, 240, 270, 288, 300, 356]
    divs = [60, 72, 90, 120, 180, 240, 270, 288, 300]  # [60, 72, 90, 120, 180, 240, 270, 288, 300]
    threads = 16
    min_spread, max_spread = 0.5, 1  # 0.1, 2  # 0.5, 0.65

    edge_buffer = 4  # Number of points to skip on either extreme to avoid boundaries

    sim_pars = find_sim_sets(f'{base_path}Data_Sim/', ['flat80', 'anticlmulti', 'resample'], ['test'], True)
    sim_amps, sim_spreads = pd.unique(sim_pars['amp']), pd.unique(sim_pars['spread'])
    print(sim_amps)

    jobs = []
    for amp in sim_amps:
        for spread in sim_spreads:
            name = f'Sim_spread{spread}_amp{amp}'
            amp_f, spread_f = get_name_amp_spread(name)
            if not min_spread < spread_f < max_spread:
                continue
            sim_pars = (base_path, f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
                        f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_0',
                        62, cent, divs, 'Data_Sim', 'Data_Sim_Mix', amp_f, spread_f)
            jobs.append(sim_pars)

    sim_sets = {}
    print('Getting sim set ')
    with Pool(threads) as pool:
        for amp_f, spread_f, set_dists in tqdm.tqdm(pool.istarmap(get_sets, jobs), total=len(jobs)):
            if amp_f not in sim_sets:
                sim_sets.update({amp_f: {}})
            sim_sets[amp_f].update({spread_f: set_dists})

    comp_amps = sorted(sim_sets.keys())[edge_buffer:-edge_buffer]

    costs = {
        'Diff2': {'func': diff2_cost, 'color': 'blue', 'min_def': [], 'min_sd': []},
        'Chi2': {'func': chi2_cost, 'color': 'green', 'min_def': [], 'min_sd': []},
        'Diff2_Weight': {'func': diff2_cut_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        'Chi2_Weight': {'func': chi2_cut_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
        # 'Diff2_Div': {'func': diff2_div_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        # 'Chi2_Div': {'func': chi2_div_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
        # 'Diff2_Div_CompDiv': {'func': diff2_div_compdiv_cost, 'color': 'red', 'min_def': [], 'min_sd': []},
        # 'Chi2_Div_CompDiv': {'func': chi2_div_compdiv_cost, 'color': 'purple', 'min_def': [], 'min_sd': []},
    }

    cost_res = {cost: {'min_def': [], 'min_def_bs': [], 'min_sd': []} for cost in costs.keys()}
    true_amps = []

    jobs = [(sim_sets, costs, comp_amp) for comp_amp in comp_amps]
    print('Compare sets')
    with Pool(threads) as pool:
        for comp_amp, res in tqdm.tqdm(pool.istarmap(get_comp, jobs), total=len(jobs)):
            true_amps.append(comp_amp)
            for cost in costs:
                for key in cost_res[cost].keys():
                    cost_res[cost][key].append(res[cost][key])

    for cost in costs:
        sigma_deviation = (np.array(cost_res[cost]['min_def']) - true_amps) / np.array(cost_res[cost]['min_sd'])
        fig_sigmas, ax_sigmas = plt.subplots()
        ax_sigmas.set_xlabel('Sigmas from True')
        ax_sigmas.set_title(f'{cost} Number of Sigmas from True Value')
        bins = max(10, len(sigma_deviation) / 100)
        ax_sigmas.hist(sigma_deviation, density=True, bins=bins, label=cost)
        norm_xs = np.linspace(*ax_sigmas.get_xlim(), 1000)
        norm_dist = stats.norm()
        ax_sigmas.plot(norm_xs, norm_dist.pdf(norm_xs), color='red', label='Standard Normal')
        ax_sigmas.legend()
        fig_sigmas.tight_layout()
        fig_sigmas.canvas.manager.set_window_title(f'{cost} Number of Sigmas from True Value')

    fig_dev_sum, ax_dev_sum = plt.subplots()
    ax_dev_sum.set_xlabel('Cost Function')
    ax_dev_sum.set_ylabel('RMS Deviation from True')
    cost_labels, dev_true_rmss = [], []
    for cost in costs:
        dev_true_rms = np.sqrt(np.mean((np.array(cost_res[cost]['min_def']) - true_amps) ** 2))
        dev_true_rms_sd = np.std(np.sqrt(np.mean((np.array(cost_res[cost]['min_def_bs']).T - true_amps) ** 2, axis=1)))
        cost_labels.append(cost)
        dev_true_rmss.append(dev_true_rms)
    ax_dev_sum.axhline(0, ls='--', color='black')
    ax_dev_sum.errorbar(cost_labels, dev_true_rmss, marker='o', yerr=dev_true_rms_sd, ls='none')
    ax_dev_sum.grid()
    fig_dev_sum.tight_layout()
    fig_dev_sum.canvas.manager.set_window_title(f'{cost} RMS Deviation from True')

    plt.show()


def get_comp(sim_sets, costs, comp_amp):
    all_amps = sorted(list(sim_sets.keys()))

    cost_res = {cost: {'min_def': None, 'min_def_bs': None, 'min_sd': None} for cost in costs.keys()}

    for spread in sim_sets[all_amps[0]].keys():
        for div in sim_sets[all_amps[0]][spread].keys():
            for total_proton in sim_sets[all_amps[0]][spread][div].keys():
                for cost in costs:
                    cost_defs, cost_bss, amps = [], [], []
                    for amp in all_amps:
                        if amp != comp_amp:  # Don't compare to self
                            amps.append(amp)
                            cost_def, cost_bs = costs[cost]['func'](sim_sets[amp][spread][div][total_proton],
                                                                    sim_sets[comp_amp][spread][div][total_proton],
                                                                    False)
                            cost_defs.append(cost_def)
                            cost_bss.append(cost_bs)

                    cost_bss = list(np.array(cost_bss).T)  # Group by amps instead of bs #
                    costs_min_res = []
                    for cost_i in [cost_defs, *cost_bss]:
                        cost_norm = (np.array(cost_i) - min(cost_i)) / (max(cost_i) - min(cost_i))
                        cost_interp = interp1d(amps, cost_norm, 'cubic', bounds_error=True, fill_value=max(cost_norm))
                        # cost_interp = interp1d(amps, cost_i, 'cubic', bounds_error=True, fill_value=max(cost_i))
                        costs_min_res.append(minimize_scalar(cost_interp, method='bounded',
                                                             bounds=(min(amps), max(amps))))
                    cost_res[cost]['min_def'] = costs_min_res[0].x
                    cost_res[cost]['min_def_bs'] = [min_res.x for min_res in costs_min_res[1:]]
                    cost_res[cost]['min_sd'] = np.std(cost_res[cost]['min_def_bs'])

    return comp_amp, cost_res


def get_set(set_pars):
    base_path, set_group, set_name, energy_set, cent, div, total_proton, raw_folder, mix_folder, name, amp, spread \
        = set_pars
    print(f'Calculate Set: {name}')

    file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
    path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
    raw_dist, raw_bs_dists = get_norm_dists(f'{base_path}{raw_folder}/{path_sufx}', total_proton)
    mix_dist, mix_bs_dists = get_norm_dists(f'{base_path}{mix_folder}/{path_sufx}', total_proton)

    return calc_set(raw_dist, mix_dist, raw_bs_dists, mix_bs_dists)


def calc_set_old(raw_dist, mix_dist, raw_bs_dists, mix_bs_dists):
    raw_sub_mix = raw_dist - mix_dist
    raw_sub_mix_bss = np.array([raw - mix for raw in raw_bs_dists for mix in mix_bs_dists])
    raw_sub_mix_sd = np.nanstd(raw_sub_mix_bss, axis=0)
    raw_sub_mix_bss_slim = np.array([raw - mix for raw, mix in zip(raw_bs_dists, mix_bs_dists)])

    raw_div_mix = raw_dist / mix_dist
    raw_div_mix_bss = np.array([raw / mix for raw in raw_bs_dists for mix in mix_bs_dists])
    raw_div_mix_sd = np.nanstd(raw_div_mix_bss, axis=0)
    raw_div_mix_bss_slim = np.array([raw / mix for raw, mix in zip(raw_bs_dists, mix_bs_dists)])

    set_dists = {'raw_sub_mix': raw_sub_mix, 'raw_sub_mix_sd': raw_sub_mix_sd,
                 'raw_sub_mix_bss_slim': raw_sub_mix_bss_slim, 'raw_div_mix': raw_div_mix,
                 'raw_div_mix_sd': raw_div_mix_sd, 'raw_div_mix_bss_slim': raw_div_mix_bss_slim,
                 'raw_norm_dist': np.array(raw_dist), 'mix_norm_dist': np.array(mix_dist)}

    return set_dists


def calc_set(raw_dist, mix_dist, raw_bss, mix_bss):
    raw_bss, mix_bss = np.array(raw_bss), np.array(mix_bss)

    raw_sub_mix = raw_dist - mix_dist
    raw_sub_mix_bss = (raw_bss[..., None] - mix_bss.T).transpose(0, 2, 1).reshape(-1, len(raw_bss[0]))

    raw_sub_mix_sd = np.nanstd(raw_sub_mix_bss)
    raw_sub_mix_bss_slim = raw_bss - mix_bss

    raw_div_mix = raw_dist / mix_dist
    raw_div_mix_bss = (raw_bss[..., None] / mix_bss.T).transpose(0, 2, 1).reshape(-1, len(raw_bss[0]))
    raw_div_mix_sd = np.nanstd(raw_div_mix_bss)
    raw_div_mix_bss_slim = raw_bss / mix_bss

    set_dists = {'raw_sub_mix': raw_sub_mix, 'raw_sub_mix_sd': raw_sub_mix_sd,
                 'raw_sub_mix_bss_slim': raw_sub_mix_bss_slim, 'raw_div_mix': raw_div_mix,
                 'raw_div_mix_sd': raw_div_mix_sd, 'raw_div_mix_bss_slim': raw_div_mix_bss_slim,
                 'raw_norm_dist': np.array(raw_dist), 'mix_norm_dist': np.array(mix_dist)}

    return set_dists


def get_sets(base_path, set_group, set_name, energy_set, cent, divs, raw_folder, mix_folder, amp_f, spread_f):
    # base_path, set_group, set_name, energy_set, cent, divs, raw_folder, mix_folder, amp_f, spread_f = set_pars
    # print(f'Calculate Set: {name}')

    np.seterr(divide='ignore', invalid='ignore')  # Needs to be done in each process separately
    set_dists = {}
    for div in divs:
        file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'
        path_sufx = f'{set_group}/{set_name}/{energy_set}GeV/{file_name}'
        raw_dists, raw_bs_distss = get_norm_dists_all(f'{base_path}{raw_folder}/{path_sufx}')
        mix_dists, mix_bs_distss = get_norm_dists_all(f'{base_path}{mix_folder}/{path_sufx}')

        total_protons = get_total_proton_list(raw_dists, mix_dists, raw_bs_distss, mix_bs_distss)

        set_dists.update({div: {}})
        for total_proton in total_protons:
            raw_dist, mix_dist, raw_bs_dists, mix_bs_dists = [x[total_proton] for x in (raw_dists, mix_dists,
                                                                                        raw_bs_distss, mix_bs_distss)]
            set_dist = calc_set(raw_dist, mix_dist, raw_bs_dists, mix_bs_dists)
            # set_dists.append(set_dist)
            set_dists[div].update({total_proton: set_dist})

    return amp_f, spread_f, set_dists


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

    for cost in costs_norm.keys():
        for i, min_res in enumerate(costs_min_res[cost]):
            if abs(min_res.x - true_amp) > 0.01:
                fig_bad, ax_bad = plt.subplots()
                ax_bad.plot(xs, costs_interp[cost][i](xs), color=costs_colors[cost])
                ax_bad.scatter(amps, costs_norm[cost][i], color=costs_colors[cost], label=cost)
                ax_cost_amp.axvline(costs_min_res[cost][i].x, color=costs_colors[cost], label=f'{cost} Minimum')

    if show:
        plt.show()


def diff2_cost(set_data, comp_set_data, plot=False):
    raw_sub_mix_comp, raw_sub_mix = comp_set_data['raw_sub_mix'], set_data['raw_sub_mix']
    raw_sub_mix_bss_slim_comp = comp_set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_bss_slim = set_data['raw_sub_mix_bss_slim']

    diff = raw_sub_mix_comp - raw_sub_mix
    diff2 = diff ** 2

    diff2_bss = []
    for i in range(len(raw_sub_mix_bss_slim_comp)):
        diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
        diff2_bs = diff_bs ** 2
        diff2_bss.append(np.nansum(diff2_bs))

    # diff2 = np.nansum(diff2) / np.sum(~np.isnan(diff2))
    diff2 = np.nansum(diff2)

    return diff2, diff2_bss


def chi2_cost(set_data, comp_set_data, plot=False):
    raw_sub_mix_comp, raw_sub_mix = comp_set_data['raw_sub_mix'], set_data['raw_sub_mix']
    raw_sub_mix_bss_slim_comp = comp_set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_bss_slim = set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_sd_comp, raw_sub_mix_sd = comp_set_data['raw_sub_mix_sd'], set_data['raw_sub_mix_sd']

    diff = raw_sub_mix_comp - raw_sub_mix
    diff2 = diff ** 2
    err2 = (raw_sub_mix_sd_comp + raw_sub_mix_sd) ** 2
    chi2 = diff2 / err2
    # print(f'diff2: {diff2}')
    # print(f'err2: {err2}')
    # print(f'chi2: {chi2}')

    chi2_bss = []
    for i in range(len(raw_sub_mix_bss_slim_comp)):
        diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
        diff2_bs = diff_bs ** 2

        chi2_bs = diff2_bs / err2  # Placeholder probably?
        chi2_bss.append(np.nansum(chi2_bs))

    # chi2 = np.nansum(chi2) / np.sum(~np.isnan(chi2))
    chi2 = np.nansum(chi2)

    return chi2, chi2_bss


def diff2_cut_cost(set_data, comp_set_data, plot=False):
    raw_sub_mix_comp, raw_sub_mix = comp_set_data['raw_sub_mix'], set_data['raw_sub_mix']
    raw_sub_mix_bss_slim_comp = comp_set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_bss_slim = set_data['raw_sub_mix_bss_slim']
    raw_dist_comp, raw_dist = comp_set_data['raw_norm_dist'], set_data['raw_norm_dist']
    mix_dist_comp, mix_dist = comp_set_data['mix_norm_dist'], set_data['mix_norm_dist']

    diff = raw_sub_mix_comp - raw_sub_mix
    diff2 = diff ** 2 * (raw_dist + raw_dist_comp + mix_dist + mix_dist_comp) / 4

    diff2_bss = []
    for i in range(len(raw_sub_mix_bss_slim_comp)):
        diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
        diff2_bs = diff_bs ** 2 * (raw_dist + raw_dist_comp + mix_dist + mix_dist_comp) / 4
        diff2_bss.append(np.nansum(diff2_bs))

    # diff2 = np.nansum(diff2) / np.sum(~np.isnan(diff2))
    diff2 = np.nansum(diff2)

    return diff2, diff2_bss


def chi2_cut_cost(set_data, comp_set_data, plot=False):
    raw_sub_mix_comp, raw_sub_mix = comp_set_data['raw_sub_mix'], set_data['raw_sub_mix']
    raw_sub_mix_bss_slim_comp = comp_set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_bss_slim = set_data['raw_sub_mix_bss_slim']
    raw_sub_mix_sd_comp, raw_sub_mix_sd = comp_set_data['raw_sub_mix_sd'], set_data['raw_sub_mix_sd']
    raw_dist_comp, raw_dist = comp_set_data['raw_norm_dist'], set_data['raw_norm_dist']
    mix_dist_comp, mix_dist = comp_set_data['mix_norm_dist'], set_data['mix_norm_dist']

    diff = raw_sub_mix_comp - raw_sub_mix
    diff2 = diff ** 2
    err2 = (raw_sub_mix_sd_comp + raw_sub_mix_sd) ** 2
    chi2 = diff2 / err2 * (raw_dist + raw_dist_comp + mix_dist + mix_dist_comp) / 4

    chi2_bss = []
    for i in range(len(raw_sub_mix_bss_slim_comp)):
        diff_bs = raw_sub_mix_bss_slim_comp[i] - raw_sub_mix_bss_slim[i]
        diff2_bs = diff_bs ** 2

        chi2_bs = diff2_bs / err2 * (raw_dist + raw_dist_comp + mix_dist + mix_dist_comp) / 4  # Placeholder probably?
        chi2_bss.append(np.nansum(chi2_bs))

    # chi2 = np.nansum(chi2) / np.sum(~np.isnan(chi2))
    chi2 = np.nansum(chi2)

    return chi2, chi2_bss


def diff2_div_cost(set_data, comp_set_data, plot=False):
    raw_div_mix_comp, raw_div_mix = comp_set_data['raw_div_mix'], set_data['raw_div_mix']
    raw_div_mix_bss_slim_comp = comp_set_data['raw_div_mix_bss_slim']
    raw_div_mix_bss_slim = set_data['raw_div_mix_bss_slim']

    diff_div = raw_div_mix_comp - raw_div_mix
    diff2_div = diff_div ** 2

    diff2_div_bss = []
    for i in range(len(raw_div_mix_bss_slim_comp)):
        diff_div_bs = raw_div_mix_bss_slim_comp[i] - raw_div_mix_bss_slim[i]
        diff2_div_bs = diff_div_bs ** 2
        diff2_div_bss.append(np.nansum(diff2_div_bs))

    diff2_div = np.nansum(diff2_div) / np.sum(~np.isnan(diff2_div))

    return diff2_div, diff2_div_bss


def chi2_div_cost(set_data, comp_set_data, plot=False):
    raw_div_mix_comp, raw_div_mix = comp_set_data['raw_div_mix'], set_data['raw_div_mix']
    raw_div_mix_bss_slim_comp = comp_set_data['raw_div_mix_bss_slim']
    raw_div_mix_bss_slim = set_data['raw_div_mix_bss_slim']
    raw_div_mix_sd_comp, raw_div_mix_sd = comp_set_data['raw_div_mix_sd'], set_data['raw_div_mix_sd']

    diff_div = raw_div_mix_comp - raw_div_mix
    diff2_div = diff_div ** 2
    err2_div = (raw_div_mix_sd_comp + raw_div_mix_sd) ** 2
    chi2_div = diff2_div / err2_div

    chi2_div_bss = []
    for i in range(len(raw_div_mix_bss_slim_comp)):
        diff_div_bs = raw_div_mix_bss_slim_comp[i] - raw_div_mix_bss_slim[i]
        diff2_div_bs = diff_div_bs ** 2

        chi2_div_bs = diff2_div_bs / err2_div  # Placeholder probably?
        chi2_div_bss.append(np.nansum(chi2_div_bs))

    chi2_div = np.nansum(chi2_div) / np.sum(~np.isnan(chi2_div))

    if plot:
        fig_comp, ax_comp = plt.subplots()
        ax_comp.plot(raw_div_mix_comp, label='comp set')
        ax_comp.plot(raw_div_mix, label='set')
        plt.show()

    return chi2_div, chi2_div_bss


def diff2_div_compdiv_cost(set_data, comp_set_data, plot=False):
    raw_div_mix_comp, raw_div_mix = comp_set_data['raw_div_mix'], set_data['raw_div_mix']
    raw_div_mix_bss_slim_comp = comp_set_data['raw_div_mix_bss_slim']
    raw_div_mix_bss_slim = set_data['raw_div_mix_bss_slim']

    diff_div = raw_div_mix_comp - raw_div_mix
    diff2_div = diff_div ** 2

    diff2_div_bss = []
    for i in range(len(raw_div_mix_bss_slim_comp)):
        diff_div_bs = raw_div_mix_bss_slim_comp[i] / raw_div_mix_bss_slim[i] - 1
        diff2_div_bs = diff_div_bs ** 2
        diff2_div_bss.append(np.nansum(diff2_div_bs))

    diff2_div = np.nansum(diff2_div) / np.sum(~np.isnan(diff2_div))

    return diff2_div, diff2_div_bss


def chi2_div_compdiv_cost(set_data, comp_set_data, plot=False):
    raw_div_mix_comp, raw_div_mix = comp_set_data['raw_div_mix'], set_data['raw_div_mix']
    raw_div_mix_bss_slim_comp = comp_set_data['raw_div_mix_bss_slim']
    raw_div_mix_bss_slim = set_data['raw_div_mix_bss_slim']
    raw_div_mix_sd_comp, raw_div_mix_sd = comp_set_data['raw_div_mix_sd'], set_data['raw_div_mix_sd']

    diff_div = raw_div_mix_comp - raw_div_mix
    diff2_div = diff_div ** 2
    err2_div = (raw_div_mix_sd_comp + raw_div_mix_sd) ** 2
    chi2_div = diff2_div / err2_div

    chi2_div_bss = []
    for i in range(len(raw_div_mix_bss_slim_comp)):
        diff_div_bs = raw_div_mix_bss_slim_comp[i] / raw_div_mix_bss_slim[i] - 1
        diff2_div_bs = diff_div_bs ** 2

        chi2_div_bs = diff2_div_bs / err2_div  # Placeholder probably?
        chi2_div_bss.append(np.nansum(chi2_div_bs))

    chi2_div = np.nansum(chi2_div) / np.sum(~np.isnan(chi2_div))

    if plot:
        fig_comp, ax_comp = plt.subplots()
        ax_comp.plot(raw_div_mix_comp, label='comp set')
        ax_comp.plot(raw_div_mix, label='set')
        plt.show()

    return chi2_div, chi2_div_bss


if __name__ == '__main__':
    main()
