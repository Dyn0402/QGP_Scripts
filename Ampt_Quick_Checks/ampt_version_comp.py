#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 28 8:01 PM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_version_comp.py

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector

from multiprocessing import Pool
import tqdm
try:
    import istarmap
except ModuleNotFoundError:
    try:
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Analyzer'))
        import istarmap
    except ModuleNotFoundError:
        print('Can\'t find istarmap!')


def main():
    base_path = 'F:/Research/'
    out_dir = 'F:/Research/Results/Presentations/7-20-22/'
    ampt_version_paths = {
        'Baryon_First': 'AMPT_Trees_Baryon_First',
        'New_Coalescence': 'AMPT_Trees',
        'Meson_First': 'AMPT_Trees_Meson_First'
    }
    mid_path = '/min_bias/string_melting/'
    # base_path = '/gpfs01/star/pwg/dneff/data/AMPT/'
    # out_dir = '/star/u/dneff/ampt_version_comp/'
    # ampt_version_paths = {
    #     'Baryon_First': 'min_bias_baryon_first',
    #     'New_Coalescence': 'min_bias',
    #     'Meson_First': 'min_bias_meson_first'
    # }
    # mid_path = '/string_melting/'
    energies = [7, 11, 19, 27, 39, 62]
    # energies = [7]

    num_files = None  # None for all
    threads = 16

    tree_name = 'tree'

    # job_func = read_file_att
    # tree_attribute = 'refmult'
    # job_pars = [tree_name, tree_attribute]
    # # bins = np.linspace(-0.1, 25, 200)  # b
    # # bins = np.arange(-0.5, 920.5, 1)  # ref3
    # bins = np.arange(-0.5, 500.5, 1)  # ref

    job_func = read_file_protons
    eta_max = None
    b_max = 5
    job_pars = [tree_name, ['pid', 'px', 'py', 'pz'], 2212, b_max, eta_max]
    tree_attribute = f'proton_dist_etamax{eta_max}_bmax{b_max}'  # Just used for file name in this case
    bins = np.arange(-0.5, 300.5, 1)

    density = True
    y_log = False

    ampt_version_dists = {}

    for energy in energies:
        fig, ax = plt.subplots()
        ax.set_title(f'{energy}GeV {tree_attribute}')
        fig.canvas.manager.set_window_title(f'{energy}GeV {tree_attribute}')
        for v_name, v_path in ampt_version_paths.items():
            print(f'Starting {energy}GeV {v_name}')
            dist = np.array([])
            path = f'{base_path}{v_path}{mid_path}{energy}GeV/'
            jobs = [(f'{path}{root_name}', *job_pars) for root_name in os.listdir(path)[:num_files]]
            with Pool(threads) as pool:
                for file_att in tqdm.tqdm(pool.istarmap(job_func, jobs), total=len(jobs)):
                    dist = np.append(dist, file_att)
            ampt_version_dists.update({v_name: np.histogram(dist, bins=bins)})
            ax.hist(dist, bins=bins, density=density, histtype='step', alpha=0.8, label=v_name)
        if y_log:
            ax.set_yscale('log')
        ax.set_xlabel(tree_attribute)
        ax.legend()
        fig.tight_layout()

        fig.savefig(f'{out_dir}ampt_{tree_attribute}_version_comp_{energy}GeV.png')
    plt.show()
    print('donzo')


def read_file_att(root_path, tree_name, tree_attribute):
    with uproot.open(root_path) as file:
        file_att = file[tree_name].arrays(tree_attribute, library='np')[tree_attribute]
    return file_att


def read_file_protons(root_path, tree_name, track_attributes, proton_pid=2212, b_max=20, eta_max=0.5):
    vector.register_awkward()
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        imps = file[tree_name].arrays(['imp'])
        tracks = tracks if b_max is None else tracks[imps.imp < b_max]
        tracks = tracks[tracks.pid == proton_pid]
        tracks = ak.zip({'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']}, with_name='Momentum3D')
        tracks = tracks if eta_max is None else tracks[abs(tracks.eta) < eta_max]
        proton_counts = ak.count(tracks.px, axis=1)

    return proton_counts


if __name__ == '__main__':
    main()
