#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 31 6:22 PM 2023
Created in PyCharm
Created as QGP_Scripts/ampt_v2

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import math

import uproot
import awkward as ak
import vector

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from Azimuthal_Correlations.az_corr import get_ampt_ref3_edges
from Sub_Resample.sub_event_resample_algorithm import get_resamples4
from Measure import Measure


def main():
    calculate_v2()
    calculate_v2_cfmodel()
    print('donzo')


def calculate_v2():
    base_path = 'F:/Research/'
    data_path = f'{base_path}AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    ampt_cent_path = f'{base_path}Ampt_Centralities_New_Coalescence/string_melting/'
    out_dir = f'{base_path}Data_Ampt_New_Coal/default_resample_epbins1/' \
              f'Ampt_rapid05_resample_norotate_epbins1_0/'
    # out_dir = None
    out_sufx = '_new'
    plot_out_dir = 'F:/Research/Results/AMPT_New_Coal_QA_Plots/Flow/'
    print_results = False
    # plot_out_dir = None
    energies = [7, 11, 19, 27, 39, 62]
    max_rapid = 0.5
    min_pt = 0.4  # GeV
    max_pt = 2.0  # GeV
    max_p = 3.0  # GeV
    eta_gap = 0.2  # Was 0.2 in first run
    proton_pid = 2212
    pids = [211, 321, -211, -321, 2212, -2212]
    # pids = [2212]
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    # cents = [8]
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%',
                0: '80-90%', -1: '90-100%'}
    calc_quantities = {'v2_ep_sum': 0, 'v2_ep_sum2': 0, 'n_v2': 0, 'ep_res_sum': 0, 'ep_res_sum2': 0, 'n_psi': 0,
                       'vn_cum_events': 0}
    rp_harmonics = [1, 2, 3, 4, 5, 6]
    for n in rp_harmonics:
        calc_quantities.update({f'v{n}_rp_sum': 0, f'v{n}_rp_sum2': 0, f'v{n}_cum_sum': 0, f'v{n}_cum_sum2': 0})
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    threads = 15

    for energy in energies:
        print(f'Starting {energy}GeV')
        file_dir = f'{data_path}{energy}GeV/'
        file_paths = os.listdir(file_dir)[:]

        ref3_edges = get_ampt_ref3_edges(ampt_cent_path, energy)

        jobs = [(f'{file_dir}{path}', read_branches, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
                 calc_quantities, rp_harmonics) for path in file_paths]

        v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
        with Pool(threads) as pool:
            for file_v2_data in tqdm.tqdm(pool.istarmap(read_file, jobs), total=len(jobs)):
                for pid in pids:
                    for cent in cents:
                        for quant in v2_data[pid][cent]:
                            v2_data[pid][cent][quant] += file_v2_data[pid][cent][quant]

        for pid in pids:
            v2_avgs, v2_sds, v2_avg_err, cent_labels = [], [], [], []
            vn_rp_avgs, vn_rp_sds, vn_rp_avg_err, vn_cum_avgs, vn_cum_sds, vn_cum_avg_err = {}, {}, {}, {}, {}, {}
            for n in rp_harmonics:
                vn_rp_avgs.update({n: []})
                vn_rp_sds.update({n: []})
                vn_rp_avg_err.update({n: []})
                vn_cum_avgs.update({n: []})
                vn_cum_sds.update({n: []})
                vn_cum_avg_err.update({n: []})
            reses, res_err = [], []
            for cent in cents:
                data = v2_data[pid][cent]
                if data['n_v2'] > 0 and data['n_psi'] > 0 and data['vn_cum_events'] > 0:
                    res = np.sqrt(data['ep_res_sum'] / data['n_psi'])
                    reses.append(res)
                    res_sum_err = np.sqrt(data['ep_res_sum2'] / data['n_psi'] - res ** 4) / \
                                  np.sqrt(data['n_psi'])
                    res_err.append(res_sum_err / (2 * res))  # Error propagation for sqrt
                    # Assume no error for now, figure out how to deal properly later
                    v2_avgs.append(data['v2_ep_sum'] / data['n_v2'] / res)
                    v2_sds.append(np.sqrt(data['v2_ep_sum2'] / data['n_v2'] - v2_avgs[-1] ** 2) / res)
                    v2_avg_err.append(v2_sds[-1] / np.sqrt(data['n_v2']))
                    for n in rp_harmonics:
                        vn_rp_avgs[n].append(data[f'v{n}_rp_sum'] / data['n_v2'])
                        vn_rp_sds[n].append(
                            np.sqrt(data[f'v{n}_rp_sum2'] / data['n_v2'] - vn_rp_avgs[n][-1] ** 2))
                        vn_rp_avg_err[n].append(vn_rp_sds[n][-1] / np.sqrt(data['n_v2']))

                        vn_cum = data[f'v{n}_cum_sum'] / data['vn_cum_events']
                        vn_cum_sd = np.sqrt(data[f'v{n}_cum_sum2'] / data['vn_cum_events'] - vn_cum ** 2)
                        vn_cum_err = vn_cum_sd / np.sqrt(data['vn_cum_events'])
                        vn_cum = np.sign(vn_cum) * np.sqrt(np.abs(Measure(vn_cum, vn_cum_err)))
                        vn_cum_avgs[n].append(vn_cum.val)
                        vn_cum_sds[n].append(np.sqrt(vn_cum_sd))
                        vn_cum_avg_err[n].append(vn_cum.err)
                    cent_labels.append(cent_map[cent])

            if plot_out_dir:
                fig, ax = plt.subplots(dpi=144)
                ax.grid()
                ax.axhline(0, color='black', ls='--')
                ax.errorbar(cent_labels, v2_avgs, yerr=v2_avg_err, ls='none', color='blue', marker='o',
                            label='Event Plane')
                ax.errorbar(cent_labels, vn_rp_avgs[2], yerr=vn_rp_avg_err[2], ls='none', color='red', marker='o',
                            label='Reaction Plane')
                ax.errorbar(cent_labels, vn_cum_avgs[2], yerr=vn_cum_avg_err[2], ls='none', color='green', marker='o',
                            label='Two Particle Correlation')

                ax.set_ylabel('v2')
                ax.set_xlabel('Centrality')
                plt.xticks(rotation=30)
                ax.legend()
                ax.set_title(f'AMPT {energy} GeV PID {pid} Elliptic Flow')
                fig.tight_layout()

                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2{out_sufx}.png')
                ax.set_ylim(bottom=-0.01, top=0.1)
                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2_zoom{out_sufx}.png')

                for n in rp_harmonics:
                    if n == 2:
                        continue
                    fig, ax = plt.subplots(dpi=144)
                    ax.grid()
                    ax.axhline(0, color='black', ls='--')
                    ax.errorbar(cent_labels, vn_rp_avgs[n], yerr=vn_rp_avg_err[n], ls='none', color='red', marker='o',
                                label='Reaction Plane')
                    ax.errorbar(cent_labels, vn_cum_avgs[n], yerr=vn_cum_avg_err[n], ls='none', color='green',
                                marker='o', label='Two Particle Correlation')

                    ax.set_ylabel(f'v{n}')
                    ax.set_xlabel('Centrality')
                    plt.xticks(rotation=30)
                    ax.legend()
                    ax.set_title(f'AMPT {energy} GeV PID {pid} v{n} Flow')
                    fig.tight_layout()

                    fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v{n}{out_sufx}.png')

            if out_dir and pid == proton_pid:
                out_path = f'{out_dir}{energy}GeV/v2{out_sufx}.txt'
                write_vn(out_path, cents, v2_avgs, v2_avg_err, reses, res_err, n=2)
                for n in rp_harmonics:
                    out_path_rp = f'{out_dir}{energy}GeV/v{n}_rp{out_sufx}.txt'
                    write_v2(out_path_rp, cents, vn_rp_avgs[n], vn_rp_avg_err[n], reses, res_err)
                    out_path_cum = f'{out_dir}{energy}GeV/v{n}_cum{out_sufx}.txt'
                    write_v2(out_path_cum, cents, vn_cum_avgs[n], vn_cum_avg_err[n], reses, res_err)
            if print_results:
                print(f'pid {pid}:\n{cents}\nv2s: {[Measure(val, err) for val, err in zip(v2_avgs, v2_avg_err)]}'
                      f'\nresolutions: {[Measure(val, err) for val, err in zip(reses, res_err)]}')


def calculate_v2_cfmodel():
    base_path = 'F:/Research/'
    cf_type = '_EVb342'  # '', '_EV', '_EVb342'
    data_path = f'{base_path}Cooper_Frye{cf_type}_Trees/All_Particles/'
    cf_type = cf_type.strip('_')
    # out_dir = f'{base_path}Data_CF{cf_type}/default_resample_epbins1/CF{cf_type}_rapid05_resample_norotate_epbins1_0/'
    out_dir = None
    out_sufx = '_new'
    # plot_out_dir = f'F:/Research/Results/Cooper_Frye{cf_type}_QA/'
    plot_out_dir = None  # Not really ready for plotting, using AMPT code to plot against centrality
    batch_size = '30000 kB'
    print_results = False
    # energies = [7, 19, 27, 39, 62]
    energies = [62]
    max_rapid = 0.5
    min_pt = 0.4  # GeV
    max_pt = 2.0  # GeV
    max_p = 3.0  # GeV
    eta_gap = 0.2  # Was 0.2 in first run
    proton_pid = 2212
    pids = [211, 321, -211, -321, 2212, -2212]
    # pids = [2212]
    cents = [8]
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%',
                0: '80-90%', -1: '90-100%'}
    calc_quantities = {'v2_ep_sum': 0, 'v2_ep_sum2': 0, 'n_v2': 0, 'ep_res_sum': 0, 'ep_res_sum2': 0, 'n_psi': 0,
                       'vn_cum_events': 0}
    rp_harmonics = [1, 2, 3, 4, 5, 6]
    for n in rp_harmonics:
        calc_quantities.update({f'v{n}_rp_sum': 0, f'v{n}_rp_sum2': 0, f'v{n}_cum_sum': 0, f'v{n}_cum_sum2': 0})
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    threads = 15

    for energy in energies:
        print(f'\nStarting {energy}GeV')
        file_dir = f'{data_path}{energy}GeV/'
        file_paths = [file_dir + file_name for file_name in os.listdir(file_dir)[:]]
        uproot_batches = uproot.iterate(file_paths, read_branches, step_size=batch_size)
        # num_events = ak.count([x for x in uproot.iterate(file_paths, ['refmult3'], step_size=batch_size)])
        test_file = uproot.open(file_paths[0])['tree']
        events_per_file = test_file.num_entries
        events_per_batch = test_file.num_entries_for(batch_size, read_branches)
        num_files = len(file_paths)
        est_total_batches = math.ceil(num_files * events_per_file / events_per_batch * 1.015)
        print(f'{num_files} files, {events_per_file} events/file, {events_per_batch} events/batch, '
              f'{est_total_batches} estimated total batches')
        # num_batches = len([1 for x in uproot_batches])

        ref3_edges = get_cf_ref3_edge()  # Overkill but matches style of AMPT v2 calculation

        def gen():
            for batch in uproot_batches:
                yield (batch, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
                       calc_quantities, rp_harmonics)

        # jobs = [(batch, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
        #          calc_quantities, rp_harmonics) for batch in uproot_batches]

        v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
        with Pool(threads) as pool:
            for file_v2_data in tqdm.tqdm(pool.istarmap(read_batch, gen()), total=est_total_batches):
                for pid in pids:
                    for cent in cents:
                        for quant in v2_data[pid][cent]:
                            v2_data[pid][cent][quant] += file_v2_data[pid][cent][quant]

        for pid in pids:
            v2_avgs, v2_sds, v2_avg_err, cent_labels = [], [], [], []
            vn_rp_avgs, vn_rp_sds, vn_rp_avg_err, vn_cum_avgs, vn_cum_sds, vn_cum_avg_err = {}, {}, {}, {}, {}, {}
            for n in rp_harmonics:
                vn_rp_avgs.update({n: []})
                vn_rp_sds.update({n: []})
                vn_rp_avg_err.update({n: []})
                vn_cum_avgs.update({n: []})
                vn_cum_sds.update({n: []})
                vn_cum_avg_err.update({n: []})
            reses, res_err = [], []
            for cent in cents:
                data = v2_data[pid][cent]
                if data['n_v2'] > 0 and data['n_psi'] > 0 and data['vn_cum_events'] > 0:
                    res = np.sqrt(data['ep_res_sum'] / data['n_psi'])
                    reses.append(res)
                    res_sum_err = np.sqrt(data['ep_res_sum2'] / data['n_psi'] - res ** 4) / np.sqrt(data['n_psi'])
                    res_err.append(res_sum_err / (2 * res))  # Error propagation for sqrt
                    # Assume no error for now, figure out how to deal properly later
                    v2_avgs.append(data['v2_ep_sum'] / data['n_v2'] / res)
                    v2_sds.append(np.sqrt(data['v2_ep_sum2'] / data['n_v2'] - v2_avgs[-1] ** 2) / res)
                    v2_avg_err.append(v2_sds[-1] / np.sqrt(data['n_v2']))
                    for n in rp_harmonics:
                        vn_rp_avgs[n].append(data[f'v{n}_rp_sum'] / data['n_v2'])
                        vn_rp_sds[n].append(
                            np.sqrt(data[f'v{n}_rp_sum2'] / data['n_v2'] - vn_rp_avgs[n][-1] ** 2))
                        vn_rp_avg_err[n].append(vn_rp_sds[n][-1] / np.sqrt(data['n_v2']))

                        vn_cum_avgs[n].append(data[f'v{n}_cum_sum'] / data['vn_cum_events'])
                        vn_cum_sds[n].append(
                            np.sqrt(data[f'v{n}_cum_sum2'] / data['vn_cum_events'] - vn_cum_avgs[n][-1] ** 2))
                        vn_cum_avg_err[n].append(vn_cum_sds[n][-1] / np.sqrt(data['vn_cum_events']))

                    cent_labels.append(cent_map[cent])

            if plot_out_dir:
                fig, ax = plt.subplots(dpi=144)
                ax.grid()
                ax.axhline(0, color='black', ls='--')
                ax.errorbar(cent_labels, v2_avgs, yerr=v2_avg_err, ls='none', color='blue', marker='o',
                            label='Event Plane')
                ax.errorbar(cent_labels, vn_rp_avgs[2], yerr=vn_rp_avg_err[2], ls='none', color='red', marker='o',
                            label='Reaction Plane')
                ax.errorbar(cent_labels, vn_cum_avgs[2], yerr=vn_cum_avg_err[2], ls='none', color='green', marker='o',
                            label='Two Particle Correlation')

                ax.set_ylabel('v2')
                ax.set_xlabel('Centrality')
                plt.xticks(rotation=30)
                ax.legend()
                ax.set_title(f'AMPT {energy} GeV PID {pid} Elliptic Flow')
                fig.tight_layout()

                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2{out_sufx}.png')
                ax.set_ylim(bottom=-0.01, top=0.1)
                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2_zoom{out_sufx}.png')

                for n in rp_harmonics:
                    if n == 2:
                        continue
                    fig, ax = plt.subplots(dpi=144)
                    ax.grid()
                    ax.axhline(0, color='black', ls='--')
                    ax.errorbar(cent_labels, vn_rp_avgs[n], yerr=vn_rp_avg_err[n], ls='none', color='red', marker='o',
                                label='Reaction Plane')
                    ax.errorbar(cent_labels, vn_cum_avgs[n], yerr=vn_cum_avg_err[n], ls='none', color='green',
                                marker='o', label='Two Particle Correlation')

                    ax.set_ylabel(f'v{n}')
                    ax.set_xlabel('Centrality')
                    plt.xticks(rotation=30)
                    ax.legend()
                    ax.set_title(f'AMPT {energy} GeV PID {pid} v{n} Flow')
                    fig.tight_layout()

                    fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v{n}{out_sufx}.png')

            if out_dir and pid == proton_pid:
                out_path = f'{out_dir}{energy}GeV/v2{out_sufx}.txt'
                write_vn(out_path, cents, v2_avgs, v2_avg_err, reses, res_err, n=2)
                for n in rp_harmonics:
                    out_path_rp = f'{out_dir}{energy}GeV/v{n}_rp{out_sufx}.txt'
                    write_v2(out_path_rp, cents, vn_rp_avgs[n], vn_rp_avg_err[n], reses, res_err)
                    out_path_cum = f'{out_dir}{energy}GeV/v{n}_cum{out_sufx}.txt'
                    write_v2(out_path_cum, cents, vn_cum_avgs[n], vn_cum_avg_err[n], reses, res_err)
            if print_results:
                print(f'pid {pid}:\n{cents}\nv2s: {[Measure(val, err) for val, err in zip(v2_avgs, v2_avg_err)]}'
                      f'\nresolutions: {[Measure(val, err) for val, err in zip(reses, res_err)]}')


def read_file(file_path, read_branches, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
              calc_quantities, rp_harmonics):
    vector.register_awkward()
    v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
    # data_nproton = {cent: {n: {'sum': 0, 'sum2': 0, 'n': 0} for n in range(100)} for cent in cents}
    ep_pids = [211, 321, -211, -321, 2212, -2212]  # Just use pions and kaons for now
    proton_mass = 0.09382721  # GeV
    with uproot.open(file_path) as file:
        # for batch in file['tree'].iterate(step_size=1000):
        tracks_all = file['tree'].arrays(read_branches)
        for cent in cents:
            tracks = tracks_all[(tracks_all.refmult3 <= ref3_edges[cent][0]) &
                                (tracks_all.refmult3 > ref3_edges[cent][1])]
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                            with_name='Momentum3D')

            tracks = tracks[(abs(tracks.eta) < 1.) & (tracks.pt > 0.2) & (tracks.pt < 2.0)]
            for pid in pids:
                non_protons = []
                ep_pids_hold = [ep_pid for ep_pid in ep_pids if abs(ep_pid) != abs(pid)]
                for ep_pid in ep_pids_hold:  # Probably a columnar way to do this but dimension is too high for my head
                    non_protons.append(tracks[tracks['pid'] == ep_pid])
                non_protons = ak.concatenate(non_protons, axis=1)
                non_protons_west = non_protons[non_protons.eta < -eta_gap]
                non_protons_east = non_protons[non_protons.eta > eta_gap]

                qx_west = ak.mean(non_protons_west.pt * np.cos(2 * non_protons_west.phi), axis=1)
                qy_west = ak.mean(non_protons_west.pt * np.sin(2 * non_protons_west.phi), axis=1)
                qx_east = ak.mean(non_protons_east.pt * np.cos(2 * non_protons_east.phi), axis=1)
                qy_east = ak.mean(non_protons_east.pt * np.sin(2 * non_protons_east.phi), axis=1)
                psi_east = np.arctan2(qy_east, qx_east) / 2.
                psi_west = np.arctan2(qy_west, qx_west) / 2.
                psi_res = np.cos(2 * (psi_east - psi_west))

                protons = tracks[(tracks['pid'] == pid)]
                protons = protons[(protons.pt > min_pt) & (protons.pt < max_pt) & (protons.p < max_p)]
                if abs(pid) == 2212:
                    proton_rapid = rapidity(protons.px, protons.py, protons.pz, protons.px * 0 + proton_mass)
                    protons = protons[abs(proton_rapid) < max_rapid]
                else:
                    protons = protons[abs(protons.eta) < max_rapid]
                protons_west = protons[protons.eta < 0].phi
                protons_east = protons[protons.eta >= 0].phi
                v2_west = np.cos(2 * (protons_west - psi_east))
                v2_east = np.cos(2 * (protons_east - psi_west))
                v2 = ak.concatenate([v2_east, v2_west], axis=1)

                v2_data[pid][cent]['v2_ep_sum'] = ak.sum(v2)
                v2_data[pid][cent]['v2_ep_sum2'] = ak.sum(v2 ** 2)
                v2_data[pid][cent]['n_v2'] = ak.count(v2)
                v2_data[pid][cent]['ep_res_sum'] = ak.sum(psi_res)
                v2_data[pid][cent]['ep_res_sum2'] = ak.sum(psi_res ** 2)
                v2_data[pid][cent]['n_psi'] = ak.count(psi_res)

                for n in rp_harmonics:
                    vn_rp = np.cos(n * protons.phi)  # AMPT reaction plane always at 0
                    v2_data[pid][cent][f'v{n}_rp_sum'] = ak.sum(vn_rp)
                    v2_data[pid][cent][f'v{n}_rp_sum2'] = ak.sum(vn_rp ** 2)

                proton_combos = ak.combinations(protons.phi, 2)
                combo_a, combo_b = ak.unzip(proton_combos)
                combos_dphi = combo_a - combo_b
                # for event_i in range(2):
                #     print(f'event {event_i}')
                #     print(f'proton_phis: {len(protons.phi[event_i])} {protons.phi[event_i]}')
                #     print(f'proton_combos: {len(proton_combos[event_i])} {proton_combos[event_i]}')
                #     print(f'combos_dphi: {len(combos_dphi[event_i])} {combos_dphi[event_i]}')
                # v2_data[pid][cent]['vn_cum_pairs'] = ak.count(combos_dphi)
                for n in rp_harmonics:
                    vn_cum = np.cos(n * combos_dphi)
                    # v2_data[pid][cent][f'v{n}_cum_sum'] = ak.sum(vn_cum)
                    # v2_data[pid][cent][f'v{n}_cum_sum2'] = ak.sum(vn_cum ** 2)
                    vn_cum_avg = ak.mean(vn_cum, axis=1)
                    v2_data[pid][cent][f'v{n}_cum_sum'] = ak.sum(vn_cum_avg)
                    v2_data[pid][cent][f'v{n}_cum_sum2'] = ak.sum(vn_cum_avg ** 2)
                    # for event_i in range(2):
                    #     print(f'event {event_i}')
                    #     print(f'v{n}_cum: {len(vn_cum[event_i])} {vn_cum[event_i]}')
                    #     print(f'v{n}_cum_avg: {vn_cum_avg[event_i]}')
                v2_data[pid][cent]['vn_cum_events'] = ak.count(vn_cum_avg)

    return v2_data


def read_batch(batch, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap, calc_quantities,
               rp_harmonics):
    vector.register_awkward()
    v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
    # data_nproton = {cent: {n: {'sum': 0, 'sum2': 0, 'n': 0} for n in range(100)} for cent in cents}
    ep_pids = [211, 321, -211, -321, 2212, -2212]  # Just use pions and kaons for now
    proton_mass = 0.09382721  # GeV

    for cent in cents:
        tracks = batch[(batch.refmult3 <= ref3_edges[cent][0]) & (batch.refmult3 > ref3_edges[cent][1])]
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')

        tracks = tracks[(abs(tracks.eta) < 1.) & (tracks.pt > 0.2) & (tracks.pt < 2.0)]
        for pid in pids:
            non_protons = []
            ep_pids_hold = [ep_pid for ep_pid in ep_pids if abs(ep_pid) != abs(pid)]
            for ep_pid in ep_pids_hold:  # Probably a columnar way to do this but dimension is too high for my head
                non_protons.append(tracks[tracks['pid'] == ep_pid])
            non_protons = ak.concatenate(non_protons, axis=1)
            non_protons_west = non_protons[non_protons.eta < -eta_gap]
            non_protons_east = non_protons[non_protons.eta > eta_gap]

            qx_west = ak.mean(non_protons_west.pt * np.cos(2 * non_protons_west.phi), axis=1)
            qy_west = ak.mean(non_protons_west.pt * np.sin(2 * non_protons_west.phi), axis=1)
            qx_east = ak.mean(non_protons_east.pt * np.cos(2 * non_protons_east.phi), axis=1)
            qy_east = ak.mean(non_protons_east.pt * np.sin(2 * non_protons_east.phi), axis=1)
            psi_east = np.arctan2(qy_east, qx_east) / 2.
            psi_west = np.arctan2(qy_west, qx_west) / 2.
            psi_res = np.cos(2 * (psi_east - psi_west))

            protons = tracks[(tracks['pid'] == pid)]
            protons = protons[(protons.pt > min_pt) & (protons.pt < max_pt) & (protons.p < max_p)]
            if abs(pid) == 2212:
                proton_rapid = rapidity(protons.px, protons.py, protons.pz, protons.px * 0 + proton_mass)
                protons = protons[abs(proton_rapid) < max_rapid]
            else:
                protons = protons[abs(protons.eta) < max_rapid]
            protons_west = protons[protons.eta < 0].phi
            protons_east = protons[protons.eta >= 0].phi
            v2_west = np.cos(2 * (protons_west - psi_east))
            v2_east = np.cos(2 * (protons_east - psi_west))
            v2 = ak.concatenate([v2_east, v2_west], axis=1)

            v2_data[pid][cent]['v2_ep_sum'] = ak.sum(v2)
            v2_data[pid][cent]['v2_ep_sum2'] = ak.sum(v2 ** 2)
            v2_data[pid][cent]['n_v2'] = ak.count(v2)
            v2_data[pid][cent]['ep_res_sum'] = ak.sum(psi_res)
            v2_data[pid][cent]['ep_res_sum2'] = ak.sum(psi_res ** 2)
            v2_data[pid][cent]['n_psi'] = ak.count(psi_res)

            for n in rp_harmonics:
                vn_rp = np.cos(n * protons.phi)  # AMPT reaction plane always at 0
                v2_data[pid][cent][f'v{n}_rp_sum'] = ak.sum(vn_rp)
                v2_data[pid][cent][f'v{n}_rp_sum2'] = ak.sum(vn_rp ** 2)

            proton_combos = ak.combinations(protons.phi, 2)
            combo_a, combo_b = ak.unzip(proton_combos)
            combos_dphi = combo_a - combo_b
            # for event_i in range(2):
            #     print(f'event {event_i}')
            #     print(f'proton_phis: {len(protons.phi[event_i])} {protons.phi[event_i]}')
            #     print(f'proton_combos: {len(proton_combos[event_i])} {proton_combos[event_i]}')
            #     print(f'combos_dphi: {len(combos_dphi[event_i])} {combos_dphi[event_i]}')
            # v2_data[pid][cent]['vn_cum_pairs'] = ak.count(combos_dphi)
            for n in rp_harmonics:
                vn_cum = np.cos(n * combos_dphi)
                # vn_cum_avg = ak.mean(vn_cum, axis=1)
                # v2_data[pid][cent][f'v{n}_cum_sum'] = ak.sum(vn_cum)
                # v2_data[pid][cent][f'v{n}_cum_sum2'] = ak.sum(vn_cum ** 2)
                vn_cum_avg = ak.mean(vn_cum, axis=1)
                v2_data[pid][cent][f'v{n}_cum_sum'] = ak.sum(vn_cum_avg)
                v2_data[pid][cent][f'v{n}_cum_sum2'] = ak.sum(vn_cum_avg ** 2)
                # for event_i in range(2):
                #     print(f'event {event_i}')
                #     print(f'v{n}_cum: {len(vn_cum[event_i])} {vn_cum[event_i]}')
                #     print(f'v{n}_cum_avg: {vn_cum_avg[event_i]}')
            v2_data[pid][cent]['vn_cum_events'] = ak.count(vn_cum_avg)

    return v2_data


def write_v2(out_path, cents, v2s, v2_errs, reses, res_errs):
    with open(out_path, 'w') as file:
        file.write('centrality\tv2\tresolution\n')
        for cent, v2, v2_err, res, res_err in zip(cents, v2s, v2_errs, reses, res_errs):
            file.write(f'{cent}\t{v2} {v2_err}\t{res} {res_err}\n')


def write_vn(out_path, cents, vns, vn_errs, reses, res_errs, n=2):
    with open(out_path, 'w') as file:
        file.write(f'centrality\tv{n}\tresolution\n')
        for cent, vn, vn_err, res, res_err in zip(cents, vns, vn_errs, reses, res_errs):
            file.write(f'{cent}\t{vn} {vn_err}\t{res} {res_err}\n')


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


def get_cf_ref3_edge():
    edges_9bin = {8: [1000, 0]}  # MUSIC+FIST model only has hypersurfaces for 0-5%

    return edges_9bin


if __name__ == '__main__':
    main()
