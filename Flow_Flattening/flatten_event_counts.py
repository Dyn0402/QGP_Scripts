#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 04 8:11 PM 2023
Created in PyCharm
Created as QGP_Scripts/flatten_event_counts

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import uproot
import awkward as ak
import vector


def main():
    finest_binning()
    # coarse_binning()
    print('donzo')


def finest_binning():
    energy = 27
    path = f'C:/Users/Dylan/phi_coefs_{energy}GeV.root'
    pp = PdfPages(f'C:/Users/Dylan/Desktop/flattener_counts_{energy}GeV.pdf')
    cents = np.arange(-1, 9, 1)
    eta_bin_num = 4
    eta_bins = np.arange(0, eta_bin_num, 1)
    hist_small_max = 1000
    tracks = {'protons': {cent: {eta: {} for eta in eta_bins} for cent in cents},
              'non-protons': {cent: {eta: {} for eta in eta_bins} for cent in cents}}
    num_tracks_list = []
    with uproot.open(path) as file:
        # print(file.keys())
        num_keys = len(file.keys())
        print(f'Number of plots {num_keys}')
        for i, key in enumerate(file.keys()):
            if not i % int(num_keys / 100):
                print(f'{float(i) / num_keys * 100:.1f}%')
            prof = file[key]
            name = prof.name
            name_split = name.split('_')
            if name_split[0] == 'cosine':
                continue
            run = int(name_split[-1])
            particle_type = name_split[2]
            cent = int(name_split[4])
            eta = int(name_split[7])
            num_tracks = prof.counts()[1]
            tracks[particle_type][cent][eta].update({run: num_tracks})
            num_tracks_list.append(num_tracks)
            # print(prof.name, name, name_split, run, prof.counts()[1])
            # input()
        # prof = file['sine_terms_non-protons_cent_-1_eta_bin_13_runkey_1113906;1']
        # print(prof.counts())
    max_num_tracks = max(num_tracks_list)
    fig, ax = plt.subplots()
    ax.hist(num_tracks_list, bins=100)
    ax.set_yscale('log')
    ax.set_xlabel('Number of tracks in phi flattening distribution')
    ax.set_title(f'{energy}GeV Tracks per Flattening Distribution Distribution')
    pp.savefig(fig)
    ax.set_yscale('linear')
    pp.savefig(fig)
    n_few_track = np.array(num_tracks_list)
    n_few_track = n_few_track[n_few_track < hist_small_max]
    fig, ax = plt.subplots()
    ax.hist(n_few_track, bins=100)
    ax.set_yscale('log')
    ax.set_xlabel('Number of tracks in phi flattening distribution')
    ax.set_title(f'{energy}GeV Truncated Tracks per Flattening Distribution Distribution')
    pp.savefig(fig)
    ax.set_yscale('linear')
    pp.savefig(fig)
    for particle_type in tracks:
        for cent in tracks[particle_type]:
            print(f'Plotting {particle_type} centrality {cent}')
            for eta in tracks[particle_type][cent]:
                runs = list(tracks[particle_type][cent][eta].keys())
                num_tracks = list(tracks[particle_type][cent][eta].values())
                fig, ax = plt.subplots()
                ax.set_yscale('log')
                ax.grid()
                ax.scatter(runs, num_tracks)
                ax.set_ylim(bottom=0.8, top=max_num_tracks)
                ax.set_xlabel('Run Number % 1000')
                ax.set_ylabel('Number of Tracks in Flattening Distribution')
                ax.set_title(f'{particle_type} {cent_convert(cent)} centrality {eta_bin_convert(eta, eta_bin_num)} eta')
                fig.tight_layout()
                pp.savefig(fig)
    pp.close()


def coarse_binning():
    path = 'C:/Users/Dylan/phi_coefs_7GeV.root'
    pp = PdfPages('C:/Users/Dylan/Desktop/flattener_counts_rebin.pdf')
    cents = np.arange(-1, 9, 1)
    eta_bin_num = 4
    eta_bins = np.arange(0, eta_bin_num, 1)
    tracks = {'protons': {cent: {eta: {} for eta in eta_bins} for cent in cents},
              'non-protons': {cent: {eta: {} for eta in eta_bins} for cent in cents}}
    with uproot.open(path) as file:
        # print(file.keys())
        num_keys = len(file.keys())
        print(f'Number of plots {num_keys}')
        for i, key in enumerate(file.keys()):
            if not i % int(num_keys / 100):
                print(f'{float(i) / num_keys * 100:.1f}%')
            prof = file[key]
            name = prof.name
            name_split = name.split('_')
            if name_split[0] == 'cosine':
                continue
            run = int(int(name_split[-1]) / 100)
            particle_type = name_split[2]
            cent = int(name_split[4])
            eta = eta_rebin(int(name_split[7]), eta_bin_num)
            num_tracks = prof.counts()[1]
            if run in tracks[particle_type][cent][eta]:
                tracks[particle_type][cent][eta][run] += num_tracks
            else:
                tracks[particle_type][cent][eta].update({run: num_tracks})
            # print(prof.name, name, name_split, run, eta, prof.counts()[1])
            # input()
        # prof = file['sine_terms_non-protons_cent_-1_eta_bin_13_runkey_1113906;1']
        # print(prof.counts())
    num_tracks_list = []
    for particle_type in tracks:
        for cent in tracks[particle_type]:
            for eta in tracks[particle_type][cent]:
                num_tracks = list(tracks[particle_type][cent][eta].values())
                print(particle_type, cent, eta, tracks[particle_type][cent][eta])
                num_tracks_list.extend(num_tracks)
    max_num_tracks = max(num_tracks_list)
    fig, ax = plt.subplots()
    ax.hist(num_tracks_list, bins=100)
    ax.set_yscale('log')
    ax.set_xlabel('Number of tracks in phi flattening distribution')
    pp.savefig(fig)
    ax.set_yscale('linear')
    pp.savefig(fig)
    for particle_type in tracks:
        for cent in tracks[particle_type]:
            print(f'Plotting {particle_type} centrality {cent}')
            for eta in tracks[particle_type][cent]:
                runs = list(tracks[particle_type][cent][eta].keys())
                num_tracks = list(tracks[particle_type][cent][eta].values())
                fig, ax = plt.subplots()
                ax.set_yscale('log')
                ax.grid()
                ax.scatter(runs, num_tracks)
                ax.set_ylim(bottom=0.8, top=max_num_tracks)
                ax.set_xlabel('Run Number / 1000')
                ax.set_ylabel('Number of Tracks in Flattening Distribution')
                ax.set_title(f'{particle_type} {cent_convert(cent)} centrality {eta_bin_convert(eta, eta_bin_num)} eta')
                fig.tight_layout()
                pp.savefig(fig)
    pp.close()


def eta_bin_convert(eta_bin, bins=20):
    eta_step = 1.0 / bins * 2
    eta_low = float(eta_bin) / bins * 2 - 1
    eta_high = eta_low + eta_step
    return f'({eta_low:.1f}, {eta_high:.1f})'


def eta_rebin(eta_bin, bins):
    return int(eta_bin / (20 / bins))


def cent_convert(cent_bin):
    if cent_bin == 8:
        return '0-5%'
    if cent_bin == 7:
        return '5-10%'
    cent_low = (6 - cent_bin) * 10 + 10
    return f'{cent_low}-{cent_low + 10}%'


if __name__ == '__main__':
    main()
