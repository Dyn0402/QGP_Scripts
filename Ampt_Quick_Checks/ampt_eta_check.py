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

from particle import Particle
import hepunits

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats


def main():
    mpl.rcParams.update({'figure.max_open_warning': 0})
    threads = 16
    # energies = [7, 11, 19, 27, 39, 62]
    # sets = [(energy, ['default', 'string_melting']) if energy in [7, 11] else (energy, ['string_melting'])
    #         for energy in energies]
    energies = [7, 11, 19, 27, 39, 62]
    file_type = 'root'
    test_num = 4
    sets = [(energy, ['string_melting']) for energy in energies]
    # min_bias_path = f'D:/Research/AMPT_Trees/min_bias/'
    # fix_most_cent_path = f'D:/Research/AMPT_Trees/ref_fix{test_num}_most_central/'
    files_path = f'F:/Research/AMPT_Trees_Baryon_First/min_bias/'
    out_name = 'Ampt_baryon_first_eta_minbias'
    pdf_out_path = f'F:/Research/Results/Presentations/7-20-22/' \
                   f'{out_name}.pdf'
    pdf_pz_out_path = f'F:/Research/Results/Presentations/7-20-22/' \
                      f'{out_name}_pz.pdf'
    tree_name = 'tree'
    track_attributes = ['pid', 'px', 'py', 'pz']
    pids = [2212, -2212, 3122, -3122]
    particle_pids = {2212: 'proton', -2212: 'anti-proton', 3122: 'lambda', -3122: 'anti-lambda'}
    # particle_masses = {pid: Particle.from_pdgid(pid).mass for pid in pids}
    bin_edges = np.linspace(-8, 8, 501)
    fig_x_lim = (-8, 8)

    fig_list = []
    fig_pz_list = []
    for energy, ampt_modes in sets:
        for ampt_mode in ampt_modes:
            path = f'{files_path}{ampt_mode}/{energy}GeV/'
            eta_hists = {pid: np.zeros(len(bin_edges) - 1) for pid in pids}
            pz_hists = {pid: np.zeros(2) for pid in pids}

            jobs = []
            for file_name in os.listdir(path):
                if file_type == 'root':
                    jobs.append((path + file_name, tree_name, track_attributes, pids, bin_edges))
                elif file_type == 'dat':
                    jobs.append((path + file_name, pids, bin_edges))

            print(f'Starting {energy}GeV {ampt_mode}:')
            func = read_dat_file if file_type == 'dat' else read_file
            with Pool(threads) as pool:
                for eta_hist, pz_hist in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
                    for pid in pids:
                        eta_hists[pid] += eta_hist[pid]
                        pz_hists[pid] += pz_hist[pid]

            bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
            bin_widths = bin_edges[1:] - bin_edges[:-1]
            for pid in pids:
                fig, ax = plt.subplots()
                ax.grid(zorder=0)
                ax.bar(bin_centers, eta_hists[pid], bin_widths, align='center', zorder=3)
                ax.set_xlabel('Pseudorapidity')
                ax.set_title(f'AMPT {energy}GeV {ampt_mode}-mode {particle_pids[pid]}s')
                stats = DistStats(dict(zip(bin_centers, eta_hists[pid])))
                ax.text(0.55, 0.92, f'mean: {stats.get_mean()}\nskewness: {stats.get_skewness()}',
                        transform=ax.transAxes)
                ax.set_xlim(fig_x_lim)
                # plt.axvline(0, color='black', ls=':')
                fig.tight_layout()
                fig_list.append(fig)

                fig_pz, ax_pz = plt.subplots()
                ax_pz.grid(zorder=0)
                ax_pz.bar([-4, 4], pz_hists[pid], [8, 8], align='center', zorder=3)
                ax_pz.set_xlabel('pz')
                ax_pz.set_title(f'AMPT {energy}GeV {ampt_mode}-mode {particle_pids[pid]}s')
                fig_pz.tight_layout()
                fig_pz_list.append(fig_pz)
            # plt.show()

    pdf = PdfPages(pdf_out_path)
    for figure in fig_list:
        # figure.savefig(f'{out_name}_{figure._suptitle.get_text()}.png')
        pdf.savefig(figure)
    pdf.close()
    pdf_pz = PdfPages(pdf_pz_out_path)
    for figure in fig_pz_list:
        pdf_pz.savefig(figure)
    pdf_pz.close()
    print('donzo')


def read_file(root_path, tree_name, track_attributes, pids, bin_edges):
    vector.register_awkward()
    eta_hists = {}
    pz_hists = {}
    with uproot.open(root_path) as file:
        # print(root_path)
        tracks = file[tree_name].arrays(track_attributes)
        for pid in pids:
            pid_tracks = tracks[tracks.pid == pid]
            mass = pid_tracks['pid'] * 0 + Particle.from_pdgid(pid).mass / hepunits.GeV
            pid_tracks = ak.zip({'pid': pid_tracks['pid'], 'px': pid_tracks['px'], 'py': pid_tracks['py'],
                                 'pz': pid_tracks['pz'], 'M': mass}, with_name='Momentum4D')
            # print(ak.ravel(pid_tracks.pz))
            eta_hists.update({pid: np.histogram(ak.to_numpy(ak.ravel(pid_tracks.eta)), bins=bin_edges)[0]})
            # print(ak.ravel(pid_tracks.pz))
            # print(type(ak.ravel(pid_tracks.pz)))
            # print(ak.to_numpy(ak.ravel(pid_tracks.pz)))
            # print(np.histogram(ak.to_numpy(ak.ravel(pid_tracks.pz)), bins=[-8, 0, 8]))
            pz_hists.update({pid: np.histogram(ak.to_numpy(ak.ravel(pid_tracks.pz)), bins=[-8, 0, 8])[0]})

    return eta_hists, pz_hists


def read_dat_file(dat_path, pids, bin_edges):
    rapid_lists = {pid: [] for pid in pids}
    pz_lists = {pid: [] for pid in pids}
    with open(dat_path, 'r') as file:
        lines = file.readlines()
        index = 0
        num_lines = len(lines)
        while index < num_lines:
            n_particles = int(lines[index].strip().split()[2])
            index += 1
            event_particle_count = 0
            while event_particle_count < n_particles:
                if index >= num_lines:
                    break
                pid, px, py, pz, m, x, y, z, t = lines[index].strip().split()
                pid = int(pid)
                px, py, pz, m = [float(x) for x in (px, py, pz, m)]
                if pid in pids:
                    rapid_lists[pid].append(rapidity(px, py, pz, m))
                    pz_lists[pid].append(pz)
                event_particle_count += 1
                index += 1

    rapid_hists = {pid: np.histogram(rapids, bins=bin_edges)[0] for pid, rapids in rapid_lists.items()}
    pz_hists = {pid: np.histogram(pzs, bins=[-8, 0, 8])[0] for pid, pzs in pz_lists.items()}

    return rapid_hists, pz_hists


def rapidity(px, py, pz, m):
    e = np.sqrt(m**2 + px**2 + py**2 + pz**2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
