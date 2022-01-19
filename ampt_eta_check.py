#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 14 6:49 PM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_eta_check.py

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

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats


def main():
    mpl.rcParams.update({'figure.max_open_warning': 0})
    threads = 15
    # energies = [7, 11, 19, 27, 39, 62]
    # sets = [(energy, ['default', 'string_melting']) if energy in [7, 11] else (energy, ['string_melting'])
    #         for energy in energies]
    energies = [7]
    sets = [(energy, ['string_melting']) for energy in energies]
    min_bias_path = f'D:/Research/AMPT_Trees/min_bias/'
    fix_most_cent_path = f'D:/Research/AMPT_Trees/ref_fix_most_central/'
    pdf_out_path = 'D:/Research/Results/Presentations/1-14-22_Ampt_eta/Ampt_eta_dists_7GeV_eta_fix.pdf'
    tree_name = 'tree'
    track_attributes = ['pid', 'px', 'py', 'pz']
    pids = [2212, -2212, 3122, -3122]
    particle_pids = {2212: 'proton', -2212: 'anti-proton', 3122: 'lambda', -3122: 'anti-lambda'}
    bin_edges = np.linspace(-8, 8, 1001)

    fig_list = []
    for energy, ampt_modes in sets:
        for ampt_mode in ampt_modes:
            path = f'{fix_most_cent_path}{ampt_mode}/{energy}GeV/'
            eta_hists = {pid: np.zeros(len(bin_edges) - 1) for pid in pids}

            jobs = []
            for root_name in os.listdir(path):
                jobs.append((path + root_name, tree_name, track_attributes, pids, bin_edges))

            print(f'Starting {energy}GeV {ampt_mode}:')
            with Pool(threads) as pool:
                for file_hist in tqdm.tqdm(pool.istarmap(read_file, jobs), total=len(jobs)):
                    for pid in pids:
                        eta_hists[pid] += file_hist[pid]

            bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
            bin_widths = bin_edges[1:] - bin_edges[:-1]
            for pid in pids:
                fig, ax = plt.subplots()
                ax.grid(zorder=0)
                ax.bar(bin_centers, eta_hists[pid], bin_widths, align='center', zorder=3)
                ax.set_xlabel('Pseudo-Rapidity')
                ax.set_title(f'AMPT {energy}GeV {ampt_mode}-mode {particle_pids[pid]}s')
                stats = DistStats(dict(zip(bin_centers, eta_hists[pid])))
                ax.text(0.55, 0.9, f'skewness: {stats.get_skewness()}', transform=ax.transAxes)
                # plt.axvline(0, color='black', ls=':')
                fig.tight_layout()
                fig_list.append(fig)
            # plt.show()

    pdf = PdfPages(pdf_out_path)
    for figure in fig_list:
        pdf.savefig(figure)
    pdf.close()
    print('donzo')


def read_file(root_path, tree_name, track_attributes, pids, bin_edges):
    vector.register_awkward()
    eta_hists = {}
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        for pid in pids:
            pid_tracks = tracks[tracks.pid == pid]
            eta_hists.update({pid: np.histogram(ak.ravel(pid_tracks.eta), bins=bin_edges)[0]})

    return eta_hists


if __name__ == '__main__':
    main()
