#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 29 3:20 PM 2021
Created in PyCharm
Created as QGP_Scripts/music_fist_net_protons

@author: Dylan Neff, dylan
"""
import matplotlib.ticker
import numpy as np
import matplotlib.pyplot as plt
import os

import uproot
import awkward as ak
import vector
import hist

from DistStats import DistStats


def main():
    vector.register_awkward()
    path = '/home/dylan/Research/CF_NetProton_Trees/'
    energies = [7, 19, 27, 39, 62]
    stats = {}
    for energy in energies:
        proton_stats, anti_proton_stats = get_stats(f'{path}{energy}GeV/')
        stats.update({energy: {'protons': proton_stats, 'anti-protons': anti_proton_stats}})

    cumulants_plot(stats)

    print('donzo')


def get_stats(path):
    track_atts = ['pid', 'px', 'py', 'pz']
    proton_pid = 2212
    anti_proton_pid = -2212
    rapid_max = 0.5
    pt_min = 0.4
    pt_max = 2.0
    m = 0.93827  # GeV, proton mass
    binning = np.arange(-0.5, 100.5)
    tree_name = 'tree'

    protons = np.zeros(len(binning) - 1, dtype=int)
    anti_protons = np.zeros(len(binning) - 1, dtype=int)

    files = 25
    file_num = 0
    for file_name in os.listdir(path):
        print(file_name)
        if '.root' not in file_name:
            continue
        root_path = os.path.join(path, file_name)
        with uproot.open(root_path) as file:
            tracks = file[tree_name].arrays(track_atts)
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz'],
                             'M': tracks['px'] * 0 + m}, with_name='Momentum4D')
            tracks = tracks[(abs(tracks.rapidity) < rapid_max) & (tracks.pt >= pt_min) & (tracks.pt <= pt_max)]  # cuts
            protons += np.histogram(np.sum(tracks.pid == proton_pid, axis=1), binning)[0]
            anti_protons += np.histogram(np.sum(tracks.pid == anti_proton_pid, axis=1), binning)[0]

        file_num += 1
        if file_num >= files:
            break

    # plt.step(binning[:-1], protons, where='post', color='blue')
    # plt.step(binning[:-1], anti_protons, where='post', color='red')
    # plt.show()

    bin_centers = (binning[1:] + binning[:-1]) / 2
    proton_stats = DistStats(dict(zip(bin_centers, protons)))
    anti_proton_stats = DistStats(dict(zip(bin_centers, anti_protons)))

    return proton_stats, anti_proton_stats


def cumulants_plot(stats):
    ranges = {'c2/c1 - 1': [-0.165, 0.03], 'c3/c2 - 1': [-0.45, 0.15]}
    y_ticks = {'c2/c1 - 1': ([-0.16, -0.12, -0.08, -0.04, 0.00], [-0.14, -0.10, -0.06, -0.02, 0.02]),
               'c3/c2 - 1': ([-0.4, -0.3, -0.2, -0.1, 0.0, 0.1], [-0.35, -0.25, -0.15, -0.05, 0.05])}
    x_range = (4, 400)
    ticks = [5, 10, 20, 50, 100, 200]

    energies = []
    c2_c1 = {p: {'low': [], 'high': []} for p in ['protons', 'anti-protons']}
    c3_c1 = {p: {'low': [], 'high': []} for p in ['protons', 'anti-protons']}

    for energy, stat in stats.items():
        energies.append(energy)
        for particle in c2_c1:
            meas31 = stat[particle].get_cumulant(3) / stat[particle].get_cumulant(1) - 1
            meas21 = stat[particle].get_cumulant(2) / stat[particle].get_cumulant(1) - 1
            for key, sign in {'low': -1, 'high': 1}.items():
                c3_c1[particle][key].append(meas31.val + sign * meas31.err)
                c2_c1[particle][key].append(meas21.val + sign * meas21.err)

    for title, y in {'c2/c1 - 1': c2_c1, 'c3/c2 - 1': c3_c1}.items():
        fig, ax = plt.subplots()
        ax.fill_between(energies, y['protons']['high'], y['protons']['low'], color='lightcoral')
        ax.fill_between(energies, y['anti-protons']['high'], y['anti-protons']['low'], color='lightgray')
        ax.set_ylim(ranges[title][0], ranges[title][1])
        ax.set_xlim(x_range)
        ax.axhline(0, color='b', ls='--')
        ax.set_xscale('log')
        ax.set_xticks(ticks)
        ax.set_yticks(y_ticks[title][0])
        ax.set_yticks(y_ticks[title][1], minor=True)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.tick_params(which='both', direction='in')
        ax.get_xaxis().set_ticks_position('both')
        ax.get_yaxis().set_ticks_position('both')
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel(title)
        ax.set_title(title)
        fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
