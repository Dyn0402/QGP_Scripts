#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 20 4:30 PM 2023
Created in PyCharm
Created as QGP_Scripts/main

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector


def main():
    corr_test()
    print('donzo')


def corr_test():
    # base_path = 'F:/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    base_path = '/media/ucla/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    energy = 62
    max_eta = 0.5
    ref3_min = 285  # Top 5%
    pid = 2212  # Proton?
    bin_width = np.deg2rad(120)
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    vector.register_awkward()

    phi_bins = np.linspace(-np.pi, np.pi, 100)
    phi_hist = np.histogram([], bins=phi_bins)[0]

    na, nb = [], []

    file_dir = f'{base_path}{energy}GeV/'
    files_paths = os.listdir(file_dir)
    for file_path in files_paths[:100]:
        print(file_path)
        # return
        with uproot.open(f'{file_dir}{file_path}') as file:
            tracks = file['tree'].arrays(read_branches)
            tracks = tracks[(tracks.refmult3 > ref3_min)]
            # print(tracks.pid)
            # tracks = tracks[(tracks.px >= pid)]  # Don't know why this doesn't work
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                            with_name='Momentum3D')
            tracks = tracks[(tracks['pid'] == pid) & (abs(tracks.eta) < max_eta)]
            # print(tracks)
            # print(len(tracks))
            # binned_protons = ak.count(tracks[(tracks.phi > 3) & (tracks.phi < 4)].phi, axis=1)
            # print(binned_protons)
            # print(len(binned_protons))
            # print(tracks.phi)
            # print(ak.flatten(tracks.phi))
            phi_hist += np.histogram(ak.flatten(tracks.phi), bins=phi_bins)[0]
            # print(tracks.phi + 2 * np.pi)
            # binned_protons2 = ak.count(tracks[((tracks.phi > 6) & (tracks.phi < 7)) |
            #                                   (tracks.phi + 2 * np.pi > 6) & (tracks.phi + 2 * np.pi < 7)].phi, axis=1)
            # print(binned_protons2)
            # print(len(binned_protons2))
            num_events = len(tracks)
            bins_a = np.random.rand(num_events) * 2 * np.pi - np.pi
            # bins_b = bins_a + np.deg2rad(65)
            bins_b = bins_a + np.deg2rad(65) * np.random.rand(num_events) + np.deg2rad(1)
            n_protons_a = ak.count(tracks[((tracks.phi > bins_a) & (tracks.phi <= bins_a + bin_width)) |
                                          (tracks.phi + 2 * np.pi > bins_a) &
                                          (tracks.phi + 2 * np.pi <= bins_a + bin_width)].phi, axis=1)
            n_protons_b = ak.count(tracks[((tracks.phi > bins_b) & (tracks.phi <= bins_b + bin_width)) |
                                          (tracks.phi + 2 * np.pi > bins_b) &
                                          (tracks.phi + 2 * np.pi <= bins_b + bin_width)].phi, axis=1)
            na.extend(n_protons_a)
            nb.extend(n_protons_b)

    print(phi_hist)
    phi_bin_centers = (phi_bins[1:] + phi_bins[:-1]) / 2
    phi_bin_widths = (phi_bins[1:] - phi_bins[:-1])
    plt.bar(phi_bin_centers, width=phi_bin_widths, height=phi_hist)

    plt.figure()
    plt.hist2d(na, nb, bins=(np.arange(-0.5, max(na) + 1), np.arange(-0.5, max(nb) + 1)))

    na, nb = np.array(na), np.array(nb)
    daa = np.mean(na ** 2) - np.mean(na) ** 2
    dbb = np.mean(nb ** 2) - np.mean(nb) ** 2
    dab = np.mean(na * nb) - np.mean(na) * np.mean(nb)
    sigma_c_2 = (daa + dbb - 2 * dab) / (np.mean(na + nb))

    print(list(zip(na, nb)))

    print(f'daa: {daa}')
    print(f'dbb: {dbb}')
    print(f'dab: {dab}')
    print(f'sigma_c_2: {sigma_c_2}')

    plt.show()


if __name__ == '__main__':
    main()
