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

import uproot
import awkward as ak
import vector

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from Azimuthal_Correlations.az_corr import get_ampt_ref3_edges
from Sub_Resample.sub_event_resample_algorithm import get_resamples4


def main():
    calculate_v2()
    print('donzo')


def calculate_v2():
    base_path = 'F:/Research/'
    data_path = f'{base_path}AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    ampt_cent_path = f'{base_path}Ampt_Centralities_New_Coalescence/string_melting/'
    out_dir = f'{base_path}Data_Ampt_New_Coal/default_resample_epbins1/' \
              f'Ampt_rapid05_resample_norotate_epbins1_0/'
    out_dir = None
    plot_out_dir = 'F:/Research/Results/AMPT_New_Coal_QA_Plots/Flow/'
    # plot_out_dir = None
    energies = [7, 11, 19, 27, 39, 62]
    max_rapid = 0.5
    min_pt = 0.4  # GeV
    max_pt = 2.0  # GeV
    max_p = 3.0  # GeV
    eta_gap = 0.4  # Was 0.2 in first run
    proton_pid = 2212
    pids = [211, 321, -211, -321, 2212, -2212]
    # pids = [2212]
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    # cents = [8]
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%',
                0: '80-90%', -1: '90-100%'}
    calc_quantities = {'v2_ep_sum': 0, 'v2_ep_sum2': 0, 'v2_rp_sum': 0, 'v2_rp_sum2': 0, 'n_v2': 0,
                       'ep_res_sum': 0, 'ep_res_sum2': 0, 'n_psi': 0}
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    threads = 11

    for energy in energies:
        print(f'Starting {energy}GeV')
        file_dir = f'{data_path}{energy}GeV/'
        file_paths = os.listdir(file_dir)[:]

        ref3_edges = get_ampt_ref3_edges(ampt_cent_path, energy)

        jobs = [(f'{file_dir}{path}', read_branches, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
                 calc_quantities) for path in file_paths]

        v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
        with Pool(threads) as pool:
            for file_v2_data in tqdm.tqdm(pool.istarmap(read_file, jobs), total=len(jobs)):
                for pid in pids:
                    for cent in cents:
                        for quant in v2_data[pid][cent]:
                            v2_data[pid][cent][quant] += file_v2_data[pid][cent][quant]

        for pid in pids:
            v2_avgs, v2_sds, v2_avg_err, v2_rp_avgs, v2_rp_sds, v2_rp_avg_err, cent_labels = [], [], [], [], [], [], []
            reses, res_err = [], []
            for cent in cents:
                data = v2_data[pid][cent]
                if data['n_v2'] <= 0:
                    continue
                res = np.sqrt(data['ep_res_sum'] / data['n_psi'])
                reses.append(res)
                res_sum_err = np.sqrt(data['ep_res_sum2'] / data['n_psi'] - res ** 4) / \
                              np.sqrt(data['n_psi'])
                res_err.append(res_sum_err / (2 * res))  # Error propagation for sqrt
                # Assume no error for now, figure out how to deal properly later
                v2_avgs.append(data['v2_ep_sum'] / data['n_v2'] / res)
                v2_sds.append(np.sqrt(data['v2_ep_sum2'] / data['n_v2'] - v2_avgs[-1] ** 2) / res)
                v2_avg_err.append(v2_sds[-1] / np.sqrt(data['n_v2']))
                v2_rp_avgs.append(data['v2_rp_sum'] / data['n_v2'])
                v2_rp_sds.append(
                    np.sqrt(data['v2_rp_sum2'] / data['n_v2'] - v2_rp_avgs[-1] ** 2) / res)
                v2_rp_avg_err.append(v2_rp_sds[-1] / np.sqrt(data['n_v2']))
                cent_labels.append(cent_map[cent])

            if plot_out_dir:
                fig, ax = plt.subplots(dpi=144)
                ax.grid()
                ax.axhline(0, color='black', ls='--')
                ax.errorbar(cent_labels, v2_avgs, yerr=v2_avg_err, ls='none', color='blue', marker='o',
                            label='Event Plane')
                ax.errorbar(cent_labels, v2_rp_avgs, yerr=v2_rp_avg_err, ls='none', color='red', marker='o',
                            label='Reaction Plane')

                ax.set_ylabel('v2')
                ax.set_xlabel('Centrality')
                plt.xticks(rotation=30)
                ax.legend()
                ax.set_title(f'AMPT {energy} GeV PID {pid} Elliptic Flow')
                fig.tight_layout()

                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2.png')
                ax.set_ylim(bottom=-0.01, top=0.1)
                fig.savefig(f'{plot_out_dir}Ampt_{energy}GeV_pid_{pid}_v2_zoom.png')

            if out_dir and pid == proton_pid:
                out_path = f'{out_dir}{energy}GeV/v2_04gap.txt'
                write_v2(out_path, cents, v2_avgs, v2_avg_err, reses, res_err)
                out_path_rp = f'{out_dir}{energy}GeV/v2_rp_04gap.txt'
                write_v2(out_path_rp, cents, v2_rp_avgs, v2_rp_avg_err, reses, res_err)


def read_file(file_path, read_branches, cents, pids, ref3_edges, max_rapid, min_pt, max_pt, max_p, eta_gap,
              calc_quantities):
    vector.register_awkward()
    v2_data = {pid: {cent: calc_quantities.copy() for cent in cents} for pid in pids}
    # data_nproton = {cent: {n: {'sum': 0, 'sum2': 0, 'n': 0} for n in range(100)} for cent in cents}
    ep_pids = [211, 321, -211, -321, 2212, -2212]  # Just use pions and kaons for now
    proton_mass = 0.09382721  # GeV
    with uproot.open(file_path) as file:
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
                v2_rp = np.cos(2 * protons.phi)  # AMPT reaction plane always at 0
                v2 = ak.concatenate([v2_east, v2_west], axis=1)

                v2_data[pid][cent]['v2_ep_sum'] = ak.sum(v2)
                v2_data[pid][cent]['v2_ep_sum2'] = ak.sum(v2 ** 2)
                v2_data[pid][cent]['v2_rp_sum'] = ak.sum(v2_rp)
                v2_data[pid][cent]['v2_rp_sum2'] = ak.sum(v2_rp ** 2)
                v2_data[pid][cent]['n_v2'] = ak.count(v2)
                v2_data[pid][cent]['ep_res_sum'] = ak.sum(psi_res)
                v2_data[pid][cent]['ep_res_sum2'] = ak.sum(psi_res ** 2)
                v2_data[pid][cent]['n_psi'] = ak.count(psi_res)

    return v2_data


def write_v2(out_path, cents, v2s, v2_errs, reses, res_errs):
    with open(out_path, 'w') as file:
        file.write('centrality\tv2\tresolution\n')
        for cent, v2, v2_err, res, res_err in zip(cents, v2s, v2_errs, reses, res_errs):
            file.write(f'{cent}\t{v2} {v2_err}\t{res} {res_err}\n')


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
