#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 04 11:21 AM 2024
Created in PyCharm
Created as QGP_Scripts/cluster_pcons_combination.py

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

import uproot
import awkward as ak
import vector

from PConsSim import PConsSim
from ClusterSim import ClusterSim
from momentum_conservation_model import plot_vectors

from Analysis_POCs.poc_functions import bin_experiment, bin_experiment_no_bs
from PConsSim import PConsSim, rotate_vector
from DistStats import DistStats
from Measure import Measure
from analyze_binom_slices import vn_divs


def main():
    full_sim()
    # cluster_test()
    print('donzo')


def full_sim():
    out_path = 'C:/Users/Dylan/Desktop/test/'
    # out_path = '/home/dylan/mom_cons_model/'
    ss = np.random.SeedSequence(42)

    n_tracks = [30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 240, 280, 320, 390, 500, 640, 800]
    # n_tracks = [30, 40, 50]
    n_tracks.reverse()
    n_events = 10000
    # n_events = 1000
    runs = [
        # {'n_protons_frac': 0.4, 'convergence_momentum': False, 'energy': 2, 'p_clust': 0.1},
        # {'n_protons_frac': 0.4, 'convergence_momentum': 0.001, 'energy': 2, 'p_clust': 0.1},

        # {'n_protons_frac': 0.4, 'convergence_momentum': 0.001, 'energy': 2, 'p_clust': False},
        # {'n_protons_frac': 0.4, 'convergence_momentum': 0.001, 'energy': 2, 'p_clust': 0.04},
        # {'n_protons_frac': 0.4, 'convergence_momentum': 0.001, 'energy': 2, 'p_clust': False},
        # {'n_protons_frac': 0.4, 'convergence_momentum': False, 'energy': 2, 'p_clust': 0.04},
        {'n_protons_frac': 0.4, 'convergence_momentum': False, 'energy': 2, 'p_clust': 0.06},
        {'n_protons_frac': 0.4, 'convergence_momentum': False, 'energy': 2, 'p_clust': 0.08},
        {'n_protons_frac': 0.4, 'convergence_momentum': False, 'energy': 2, 'p_clust': 0.1},
    ]
    y_max, pt_max, p_max = 0.5, 2, 2

    bin_widths = np.deg2rad([60, 72, 90, 120, 180, 288])
    resamples = 72

    threads = 12

    variations = len(runs)

    variation_num = 1
    for run in runs:
        n_protons_frac = run['n_protons_frac']
        convergence_momentum = run['convergence_momentum']
        energy = run['energy']
        p_clust = run['p_clust']
        var_dir = (f'{out_path}pfrac{int(n_protons_frac * 100 + 0.1)}_convp{int(convergence_momentum * 1000 + 0.1)}_'
                   f'energy{energy}_pclust{int(p_clust*100)}_nevent{int(n_events / 1000)}k/')
        if os.path.exists(var_dir):
            os.rmdir(var_dir)
        os.mkdir(var_dir)
        print(f'\n\nStarting variation {variation_num}/{variations} ({var_dir.split("/")[-2]}):\n')
        for n_i, n_track in enumerate(n_tracks):
            print(f'Starting {n_track} track events {n_i + 1}/{len(n_tracks)}:')
            n_protons = int(n_track * n_protons_frac + 0.5)
            rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_events)])
            jobs = [(n_track, next(rngs), energy, n_protons, y_max, pt_max, p_max, convergence_momentum, p_clust)
                    for i in range(n_events)]
            experiment_tracks = {}
            with Pool(threads) as pool:
                for tracks in tqdm.tqdm(pool.istarmap(sim_event, jobs), total=len(jobs)):
                    if len(tracks) not in experiment_tracks:
                        experiment_tracks.update({len(tracks): [tracks]})
                    else:
                        experiment_tracks[len(tracks)].append(tracks)

            n_protons_list = sorted([x for x in experiment_tracks.keys() if x > 1])
            n_rndms = len(n_protons_list) * len(bin_widths)
            rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_rndms)])
            for bin_width in bin_widths:
                dsig_2_list, dsig_2_err_list, v1_list, v2_list, v3_list = [], [], [], [], []
                for n_proton in n_protons_list:
                    p = bin_width / (2 * np.pi)
                    binom_var = n_proton * p * (1 - p)
                    data = bin_experiment_no_bs(experiment_tracks[n_proton], n_proton, bin_width, resamples,
                                                next(rngs), alg=3, sort=True)
                    stats = DistStats(data, unbinned=False)
                    delta_sig_2 = (stats.get_k_stat(2) - binom_var) / (n_proton * (n_proton - 1))
                    dsig_2_list.append(delta_sig_2.val)
                    dsig_2_err_list.append(delta_sig_2.err * np.sqrt(resamples))
                    track_combos = ak.combinations(ak.Array(experiment_tracks[n_proton]), 2)
                    combo_a, combo_b = ak.unzip(track_combos)
                    combos_dphi = combo_a - combo_b
                    v1_list.append(ak.mean(np.cos(combos_dphi)))
                    v2_list.append(ak.mean(np.cos(2 * combos_dphi)))
                    v3_list.append(ak.mean(np.cos(3 * combos_dphi)))

                with open(f'{var_dir}bin_width{int(np.rad2deg(bin_width) + 0.1)}.txt', 'a') as file:
                    file.write(f'Total Particles {n_track}\n')
                    file.write(', '.join([str(x) for x in n_protons_list]) + '\n')
                    file.write(', '.join([str(x) for x in dsig_2_list]) + '\n')
                    file.write(', '.join([str(x) for x in dsig_2_err_list]) + '\n')
                    file.write(', '.join([str(x) for x in v1_list]) + '\n')
                    file.write(', '.join([str(x) for x in v2_list]) + '\n')
                    file.write(', '.join([str(x) for x in v3_list]) + '\n')
                    file.write('\n')
        variation_num += 1


def sim_event(n_track, rng, energy, n_protons, y_max, pt_max, p_max, convergence_momentum, p_clust):
    vector.register_awkward()
    event = PConsSim(energy, n_track, rng=rng)
    if convergence_momentum is not False:
        event.convergence_momentum = convergence_momentum
        event.adaptive_alpha = True
        event.rotate_tracks()

        while event.init_net_momentum <= event.net_momentum:
            alpha = event.alpha
            event = PConsSim(energy, n_track, rng=rng)
            event.alpha = alpha / 2
            event.rotate_tracks()

    tracks = event.get_m_tracks(n_protons)
    if p_clust is not False:
        cluster_sim = ClusterSim(tracks, p_clust, rng=rng)
        cluster_sim.cluster_tracks()
        tracks = cluster_sim.tracks

    tracks = ak.Array({'px': tracks[:, 0], 'py': tracks[:, 1], 'pz': tracks[:, 2], 'M': tracks[:, 0] * 0},
                      with_name='Momentum4D')
    tracks = tracks[(abs(tracks.rapidity) < y_max) & (tracks.pt < pt_max) & (tracks.p < p_max)]
    tracks = np.array(tracks.phi)
    tracks[tracks < 0] += 2 * np.pi

    return tracks


def cluster_test():
    energy = 10
    n_particles = 20
    n_protons = 8
    dim = 3
    p_clust = 0.9
    frac_rotate = 0.95
    convergence_momentum = 0.01
    rng = np.random.default_rng()

    event = PConsSim(energy, n_particles, dim, rng=rng)
    event.convergence_momentum = convergence_momentum
    event.adaptive_alpha = True
    event.rotate_tracks()

    while event.init_net_momentum <= event.net_momentum:
        alpha = event.alpha
        event = PConsSim(energy, n_particles, rng=rng)
        event.alpha = alpha / 2
        event.rotate_tracks()

    proton_tracks = event.get_m_tracks(n_protons)
    cluster_sim = ClusterSim(proton_tracks, p_clust, frac_rotate, rng)
    cluster_sim.cluster_tracks()
    clustered_proton_tracks = cluster_sim.tracks

    print(f'Initial Proton Tracks: {proton_tracks}')
    print(f'Clustered Proton Tracks: {clustered_proton_tracks}')

    plot_vectors(event.momenta, 'Initial Tracks')
    plot_vectors(proton_tracks, 'Initial Proton Tracks')
    plot_vectors(clustered_proton_tracks, 'Clustered Proton Tracks')
    plt.show()


if __name__ == '__main__':
    main()
