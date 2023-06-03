#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 03 11:29 AM 2023
Created in PyCharm
Created as QGP_Scripts/write_test_samples.py

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector

from Azimuthal_Correlations.az_corr import get_ampt_ref3_edges


def main():
    base_in_path = 'F:/Research/'
    in_path = f'{base_in_path}AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    out_path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Jared_Tejes_Sample_Data/'
    ampt_cent_path = f'{base_in_path}Ampt_Centralities_New_Coalescence/string_melting/'
    energies = [7, 11, 19, 27, 39, 62]
    cent = 8  # Centrality bin
    n_events = 1000

    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    max_rapid = 0.5
    min_pt = 0.4  # GeV
    max_pt = 2.0  # GeV
    max_p = 3.0  # GeV
    proton_pid = 2212
    proton_mass = 0.09382721  # GeV

    out_data = {energy: [] for energy in energies}

    vector.register_awkward()
    for energy in energies:
        ref3_edges = get_ampt_ref3_edges(ampt_cent_path, energy)
        energy_dir = f'{in_path}{energy}GeV/'
        for file_name in os.listdir(energy_dir):
            if '.root' in file_name:
                with uproot.open(f'{energy_dir}{file_name}') as file:
                    tracks_all = file['tree'].arrays(read_branches)
                    tracks = tracks_all[(tracks_all.refmult3 <= ref3_edges[cent][0]) &
                                        (tracks_all.refmult3 > ref3_edges[cent][1])]
                    tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                                    with_name='Momentum3D')

                    protons = tracks[(tracks['pid'] == proton_pid)]
                    protons = protons[(protons.pt > min_pt) & (protons.pt < max_pt) & (protons.p < max_p)]
                    proton_rapid = rapidity(protons.px, protons.py, protons.pz, protons.px * 0 + proton_mass)
                    protons = protons[abs(proton_rapid) < max_rapid]
                    proton_rapid = rapidity(protons.px, protons.py, protons.pz, protons.px * 0 + proton_mass)
                    proton_phi = protons.phi
                    for event_index in range(len(proton_rapid)):
                        event_rapid, event_phi = np.array(proton_rapid[event_index]), np.array(proton_phi[event_index])
                        event_phi[event_phi < 0] += 2 * np.pi
                        event_rapid_phi = np.array(list(zip(event_rapid, event_phi)))
                        # print(event_rapid_phi)
                        out_data[energy].append(event_rapid_phi)
            if len(out_data[energy]) >= n_events:
                write_data(out_path, energy, out_data[energy][:n_events])
                break

    print('donzo')


def write_data(out_path, energy, data):
    with open(f'{out_path}{energy}GeV.txt', 'w') as file:
        file.write('Rapidity, Phi\n')
        for event_index, event in enumerate(data):
            file.write(f'\nEvent {event_index}  {len(event)} protons\n')
            for track in event:
                file.write(f'{track[0]}, {track[1]}\n')


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
