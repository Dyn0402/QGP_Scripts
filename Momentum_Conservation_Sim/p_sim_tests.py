#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 29 12:02 PM 2023
Created in PyCharm
Created as QGP_Scripts/p_sim_tests

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector

from Analysis_POCs.poc_functions import bin_experiment, bin_experiment_no_bs
from PConsSim import PConsSim
from DistStats import DistStats


def main():
    vector.register_awkward()
    rng = np.random.default_rng(42)

    n_tracks = 100
    n_protons = 40
    n_events = 1000
    energy = 2  # Currently just the range of momenta
    y_max, pt_max, p_max = 0.5, 2, 2

    bin_width = np.deg2rad(120)
    resamples = 72

    experiment_tracks = {}
    for event_i in range(n_events):
        event = PConsSim(energy, n_tracks)
        # event.rotate_tracks()
        tracks = event.get_m_tracks(n_protons)
        tracks = ak.Array({'px': tracks[:, 0], 'py': tracks[:, 1], 'pz': tracks[:, 2], 'M': tracks[:, 0] * 0},
                          with_name='Momentum4D')
        # print(tracks.p)
        tracks = tracks[(abs(tracks.rapidity) < y_max) & (tracks.pt < pt_max) & (tracks.p < p_max)]
        # print(len(tracks))
        # tracks = np.array(tracks.phi)
        tracks = np.array(np.random.uniform(0, 2 * np.pi, size=len(tracks)))
        if len(tracks) not in experiment_tracks:
            experiment_tracks.update({len(tracks): [tracks]})
        else:
            experiment_tracks[len(tracks)].append(tracks)

    n_protons_list = sorted([x for x in experiment_tracks.keys() if x > 1])
    dsig_2_list, dsig_2_err_list = [], []
    for n_proton in n_protons_list:
        p = bin_width / (2 * np.pi)
        binom_var = n_proton * p * (1 - p)
        data = bin_experiment_no_bs(experiment_tracks[n_proton], n_proton, bin_width, resamples, rng, alg=3)
        stats = DistStats(data, unbinned=False)
        delta_sig_2 = (stats.get_k_stat(2) - binom_var) / (n_proton * (n_proton - 1))
        dsig_2_list.append(delta_sig_2.val)
        dsig_2_err_list.append(delta_sig_2.err)
        print(f'{n_proton} protons: mean: {stats.get_mean()}, variance: {stats.get_variance()}, delta_sig_2: {delta_sig_2}')

    plt.axhline(0, color='black')
    plt.grid()
    plt.errorbar(n_protons_list, dsig_2_list, yerr=dsig_2_err_list, marker='o', ls='none')
    plt.show()

    print('donzo')


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
