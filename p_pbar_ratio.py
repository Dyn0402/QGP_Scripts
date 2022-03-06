#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 05 8:09 PM 2022
Created in PyCharm
Created as QGP_Scripts/p_pbar_ratio.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

import uproot
import awkward as ak
import vector
import os

from particle import Particle
import hepunits

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats


def main():
    mpl.rcParams.update({'figure.max_open_warning': 0})
    threads = 16
    data_sets = ['ampt_new', 'ampt_old']  # ['ampt_old', 'ampt_new', 'bes']
    energies = [7, 11, 19, 27, 39, 62]
    base_path = 'D:/Research/'
    out_name = 'Ampt_old_new_comp'
    pdf_out_path = f'D:/Research/Results/Presentations/3-5-22_p_pbar/' \
                   f'{out_name}.pdf'
    # pdf_pz_out_path = f'D:/Research/Results/Presentations/2-15-22_Ampt_eta_asym/' \
    #                   f'{out_name}_pz.pdf'
    tree_name = 'tree'
    track_attributes_ampt = ['pid', 'px', 'py', 'pz', 'refmult3']
    track_attributes_bes = ['proton.charge', 'proton.pt', 'proton.eta', 'proton.nsigma', 'proton.beta', 'proton.dca',
                            'proton.nhits_fit', 'refmult3']
    pids = {'p': 2212, 'p_bar': -2212}
    p_edges = np.arange(-0.5, 100.5)
    p_bar_edges = np.arange(-0.5, 25.5)
    bins = {'p': p_edges, 'p_bar': p_bar_edges}
    ref3_cut = {'bes': {7: 271, 11: 344, 19: 449, 27: 491, 39: 523, 62: 572},
                'ampt_old': {7: 240, 11: 351, 19: 423, 27: 461, 39: 513, 62: 580},
                'ampt_new': {7: 241, 11: 355, 19: 428, 27: 467, 39: 519, 62: 587}}
    base_paths = {'bes': f'{base_path}BES1_Trees/',
                  'ampt_old': f'{base_path}AMPT_Old_Trees/min_bias/string_melting/',
                  'ampt_new': f'{base_path}AMPT_Trees/min_bias/string_melting/'}
    cuts = {'refmult3': None, 'rapid_min': -0.5, 'rapid_max': 0.5, 'p_min': 0.1, 'pt_min': 0.4, 'pt_max': 2, 'min_beta': 0,
            'max_beta': 1, 'min_dca': 0, 'max_dca': 1, }

    fig_list = []
    for data_set in data_sets:
        for energy in energies:
            path = f'{base_paths[data_set]}{energy}GeV/'
            cuts['refmult3'] = ref3_cut[data_set][energy]

            p_hist, p_bar_hist = np.zeros(len(p_edges) - 1), np.zeros(len(p_bar_edges) - 1)
            p_pbar_hist = np.zeros((len(p_edges) - 1, len(p_bar_edges) - 1))

            jobs = []
            for file_name in os.listdir(path):
                if 'ampt' in data_set:
                    jobs.append((path + file_name, tree_name, track_attributes_ampt, pids, cuts, bins))
                else:
                    jobs.append((path + file_name, tree_name, track_attributes_bes, pids, cuts, bins))

            print(f'Starting {data_set} {energy}GeV:')
            func = read_ampt_file if 'ampt' in data_set else read_bes_file
            with Pool(threads) as pool:
                for p_hist_i, p_bar_hist_i, p_pbar_hist2d_i in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
                    p_hist += p_hist_i
                    p_bar_hist += p_bar_hist_i
                    p_pbar_hist += p_pbar_hist2d_i

            p_bin_centers = (bins['p'][1:] + bins['p'][:-1]) / 2
            p_bin_widths = bins['p'][1:] - bins['p'][:-1]
            p_bar_bin_centers = (bins['p_bar'][1:] + bins['p_bar'][:-1]) / 2
            p_bar_bin_widths = bins['p_bar'][1:] - bins['p_bar'][:-1]

            fig_p, ax_p = plt.subplots()
            ax_p.grid(zorder=0)
            ax_p.bar(p_bin_centers, p_hist, p_bin_widths, align='center', zorder=3)
            ax_p.set_xlabel('Protons per Event')
            ax_p.set_title(f'{data_set} {energy}GeV Protons')
            stats_p = DistStats(dict(zip(p_bin_centers, p_hist)), debug=False)
            ax_p.text(0.05, 0.88,
                      f'mean: {stats_p.get_mean()}\nstd: {stats_p.get_sd()}\nskewness: {stats_p.get_skewness()}',
                      transform=ax_p.transAxes)
            # ax.set_xlim(fig_x_lim)
            # plt.axvline(0, color='black', ls=':')
            fig_p.tight_layout()
            fig_list.append(fig_p)

            fig_p_bar, ax_p_bar = plt.subplots()
            ax_p_bar.grid(zorder=0)
            ax_p_bar.bar(p_bar_bin_centers, p_bar_hist, p_bar_bin_widths, align='center', zorder=3)
            ax_p_bar.set_xlabel('Anti-Protons per Event')
            ax_p_bar.set_title(f'{data_set} {energy}GeV Anti-Protons')
            stats_p_bar = DistStats(dict(zip(p_bar_bin_centers, p_bar_hist)), debug=False)
            ax_p_bar.text(0.05, 0.88, f'mean: {stats_p_bar.get_mean()}\nstd: {stats_p_bar.get_sd()}\n'
                                      f'skewness: {stats_p_bar.get_skewness()}', transform=ax_p_bar.transAxes)
            # ax.set_xlim(fig_x_lim)
            # plt.axvline(0, color='black', ls=':')
            fig_p_bar.tight_layout()
            fig_list.append(fig_p_bar)

    pdf = PdfPages(pdf_out_path)
    for figure in fig_list:
        # figure.savefig(f'{out_name}_{figure._suptitle.get_text()}.png')
        pdf.savefig(figure)
    pdf.close()

    print('donzo')


def read_ampt_file(root_path, tree_name, track_attributes, pids, cuts, bins):
    vector.register_awkward()
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        tracks = tracks[tracks.refmult3 > cuts['refmult3']]
        mass = tracks['pid'] * 0 + Particle.from_pdgid(pids['p']).mass / hepunits.GeV
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'],
                         'pz': tracks['pz'], 'M': mass}, with_name='Momentum4D')
        tracks = tracks[(tracks.rapidity > cuts['rapid_min']) & (tracks.rapidity < cuts['rapid_max']) &
                        (tracks.mag > cuts['p_min']) & (tracks.pt > cuts['pt_min']) & (tracks.pt < cuts['pt_max'])]
        p_tracks = ak.count(tracks[tracks.pid == pids['p']].pid, axis=1)
        p_bar_tracks = ak.count(tracks[tracks.pid == pids['p_bar']].pid, axis=1)

        p_hist, p_bar_hist, p_pbar_hist2d = get_hists(p_tracks, p_bar_tracks, bins)

    return p_hist, p_bar_hist, p_pbar_hist2d


def read_bes_file(root_path, tree_name, track_attributes, pids, cuts, bins):
    vector.register_awkward()
    # cuts = {'refmult3': None, 'rapid_min': -0.5, 'rapid_max': 0.5, 'p_min': 0.1, 'pt_min': 0.4, 'pt_max': 2, 'min_beta': 0,
    #         'max_beta': 1, 'min_dca': 0, 'max_dca': 1, }
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        print(tracks)
        tracks = tracks[tracks.refmult3 > cuts['refmult3']]
        tracks = tracks[tracks['proton.pt'] > 1]
    #     tracks = tracks[tracks.refmult3 > cuts['refmult3']]
    #     mass = tracks['pid'] * 0 + Particle.from_pdgid(pids['p']).mass / hepunits.GeV
    #     tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'],
    #                      'pz': tracks['pz'], 'M': mass}, with_name='Momentum4D')
    #     tracks = tracks[(tracks.eta < cuts['rapid_max']) & (tracks.mag > cuts['p_min']) & (tracks.pt > cuts['pt_min']) &
    #                     (tracks.pt < cuts['pt_max'])]
    #     p_tracks = ak.count(tracks[tracks.pid == pids['p']].pid, axis=1)
    #     p_bar_tracks = ak.count(tracks[tracks.pid == pids['p_bar']].pid, axis=1)
    #
    #     p_hist, p_bar_hist, p_pbar_hist2d = get_hists(p_tracks, p_bar_tracks, bins)
    #
    # return p_hist, p_bar_hist, p_pbar_hist2d


def get_hists(p_tracks, p_bar_tracks, bins):
    p_hist = np.histogram(p_tracks, bins=bins['p'])[0]
    p_bar_hist = np.histogram(p_bar_tracks, bins=bins['p_bar'])[0]
    p_pbar_hist2d = np.histogram2d(p_tracks, p_bar_tracks, [bins['p'], bins['p_bar']])[0]

    return p_hist, p_bar_hist, p_pbar_hist2d


def rapidity(px, py, pz, m):
    e = np.sqrt(m**2 + px**2 + py**2 + pz**2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
