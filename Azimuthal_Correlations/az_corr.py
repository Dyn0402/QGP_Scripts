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
import pandas as pd
import matplotlib.pyplot as plt

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

import uproot
import awkward as ak
import vector

from Measure import Measure


def main():
    # corr_test()
    sigma_vs_mult()
    # random_test()
    print('donzo')


def corr_test():
    base_path = 'F:/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    # base_path = '/media/ucla/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
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
            # bins_b = bins_a + np.deg2rad(65) * np.random.rand(num_events) + np.deg2rad(1)
            bins_b = bins_a + bin_width
            # test_event = tracks[0]
            # print(f'Number of tracks: {len(test_event)}')
            # for track in test_event:
            #     print(track)
            # print(test_event)
            # print(f'bins_a: {bins_a}')
            # print(f'test_event phi: {test_event.phi}')
            # print(f'bins_a[0]: {bins_a[0]}')
            # print(test_event.phi > bins_a[0])
            # print(ak.sum(test_event.phi > bins_a[0]))
            # test_events = tracks[0:3]
            # print(f'bins_a[0:2]: {bins_a[0:3]}')
            # print(f'bins_a[0:2] + 1: {bins_a[0:3] + 1}')
            # print(f'test_events.phi: {test_events.phi}')
            # print(f'ak.count(test_events.phi): {ak.count(test_events.phi, axis=1)}')
            # print(test_events.phi > bins_a[0:3])
            # print(ak.sum(test_events.phi > bins_a[0:3]))
            # print(ak.sum(test_events.phi > bins_a[0:3], axis=1))
            # print(ak.sum((test_events.phi > bins_a[0:3]) & (test_events.phi <= bins_a[0:3] + 1), axis=1))
            # return
            n_protons_a = ak.sum(((tracks.phi > bins_a) & (tracks.phi <= bins_a + bin_width)) |
                                 (tracks.phi + 2 * np.pi > bins_a) & (tracks.phi + 2 * np.pi <= bins_a + bin_width),
                                 axis=1)
            n_protons_b = ak.sum(((tracks.phi > bins_b) & (tracks.phi <= bins_b + bin_width)) |
                                 (tracks.phi + 2 * np.pi > bins_b) & (tracks.phi + 2 * np.pi <= bins_b + bin_width),
                                 axis=1)
            # n_protons_a = ak.count(tracks[((tracks.phi > bins_a) & (tracks.phi <= bins_a + bin_width)) |
            #                               (tracks.phi + 2 * np.pi > bins_a) &
            #                               (tracks.phi + 2 * np.pi <= bins_a + bin_width)].phi, axis=1)
            # n_protons_b = ak.count(tracks[((tracks.phi > bins_b) & (tracks.phi <= bins_b + bin_width)) |
            #                               (tracks.phi + 2 * np.pi > bins_b) &
            #                               (tracks.phi + 2 * np.pi <= bins_b + bin_width)].phi, axis=1)
            na.extend(n_protons_a)
            nb.extend(n_protons_b)

    # print(phi_hist)
    phi_bin_centers = (phi_bins[1:] + phi_bins[:-1]) / 2
    phi_bin_widths = (phi_bins[1:] - phi_bins[:-1])
    plt.bar(phi_bin_centers, width=phi_bin_widths, height=phi_hist)

    plt.figure()
    plt.hist2d(na, nb, bins=(np.arange(-0.5, max(na) + 1), np.arange(-0.5, max(nb) + 1)), cmin=1, cmap='jet')

    na, nb = np.array(na), np.array(nb)
    daa = np.mean(na ** 2) - np.mean(na) ** 2
    dbb = np.mean(nb ** 2) - np.mean(nb) ** 2
    dab = np.mean(na * nb) - np.mean(na) * np.mean(nb)
    sigma_c_2 = (daa + dbb - 2 * dab) / (np.mean(na + nb))

    # print(list(zip(na, nb)))

    print(f'daa: {daa}')
    print(f'dbb: {dbb}')
    print(f'dab: {dab}')
    print(f'sigma_c_2: {sigma_c_2}')
    print(f'1 - sigma_c_2: {1 - sigma_c_2}')

    plt.show()


def get_ampt_ref3_edges(ampt_cent_path, energy):
    with open(f'{ampt_cent_path}{energy}GeV_Ampt_ref_bin_edge.txt', 'r') as file:
        edges = file.readlines()[2].strip().split()
        ref3_edges = ampt_str_edges_to_9(edges)
    return ref3_edges


def sigma_vs_mult():
    base_path = 'F:/Research/'
    # base_path = '/media/ucla/Research/'
    # base_path = '/media/ucla/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    data_path = f'{base_path}AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    ampt_cent_path = f'{base_path}Ampt_Centralities_New_Coalescence/string_melting/'
    fig_out_path = f'{base_path}Results/Two_Partition_Corr/'
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    energies = [7, 11, 19, 27, 39, 62]
    max_rapid = 0.5
    min_pt = 0.2  # GeV
    max_pt = 2.0  # GeV
    pid = 2212  # Proton?
    n_bootstraps = 100
    min_events = 3
    seed = 45
    bin_width = np.deg2rad(120)
    bin_sep = np.deg2rad(180)
    plot = False
    threads = 12
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']
    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}
    vector.register_awkward()

    df = []
    phi_bins = np.linspace(-np.pi, np.pi, 100)
    phi_bin_centers = (phi_bins[1:] + phi_bins[:-1]) / 2
    phi_bin_widths = (phi_bins[1:] - phi_bins[:-1])

    for energy in energies:
        print(f'\nStarting {energy}GeV:')
        ref3_edges = get_ampt_ref3_edges(ampt_cent_path, energy)

        nas, nbs = {cent: {} for cent in cents}, {cent: {} for cent in cents}
        phi_hists = {cent: {} for cent in cents}

        file_dir = f'{data_path}{energy}GeV/'
        files_paths = os.listdir(file_dir)

        seeds = iter(np.random.SeedSequence(seed).spawn(len(files_paths)))
        jobs = [(f'{file_dir}{path}', read_branches, cents, ref3_edges, bin_width, bin_sep, pid, max_rapid, min_pt,
                 max_pt, phi_bins, next(seeds)) for path in files_paths]

        with Pool(threads) as pool:
            for file_nas, file_nbs, file_phi_hists in tqdm.tqdm(pool.istarmap(read_file, jobs), total=len(jobs)):
                for cent in cents:
                    for n_proton in file_nas[cent]:
                        if n_proton in nas[cent]:
                            nas[cent][n_proton].extend(file_nas[cent][n_proton])
                            nbs[cent][n_proton].extend(file_nbs[cent][n_proton])
                            phi_hists[cent][n_proton] += file_phi_hists[cent][n_proton]
                        else:
                            nas[cent].update({n_proton: file_nas[cent][n_proton]})
                            nbs[cent].update({n_proton: file_nbs[cent][n_proton]})
                            phi_hists[cent].update({n_proton: file_phi_hists[cent][n_proton]})

        if plot:
            for cent in cents:
                for n_proton in nas:
                    na, nb = np.array(nas[cent][n_proton]), np.array(nbs[cent][n_proton])

                    plt.figure()
                    plt.bar(phi_bin_centers, width=phi_bin_widths, height=phi_hists[cent][n_proton])

                    plt.figure()
                    plt.hist2d(na, nb, bins=(np.arange(-0.5, max(na) + 1), np.arange(-0.5, max(nb) + 1)), cmin=1,
                               cmap='jet')

        jobs = [(np.array(nas[cent][n_proton]), np.array(nbs[cent][n_proton]), n_bootstraps, cent, n_proton)
                for cent in cents for n_proton in nas[cent] if len(nas[cent][n_proton]) >= min_events]
        with Pool(threads) as pool:
            for sig, sig_err, cent, n_proton in tqdm.tqdm(pool.istarmap(calc_bootstraps, jobs), total=len(jobs)):
                df.append({'energy': energy, 'cent': cent, 'n': n_proton, 'sig': sig, 'sig_err': sig_err})

    df = pd.DataFrame(df)
    for energy in energies:
        df_energy = df[df['energy'] == energy]
        fig_energy, ax_energy = plt.subplots(dpi=144)
        ax_energy.grid()
        ax_energy.axhline(1, color='black', ls='--')
        ax_energy.axhline(0, color='black', ls='-')
        ax_energy.set_xlabel('Total Protons in Event')
        ax_energy.set_ylabel(r'$\sigma(c)^2$')
        ax_energy.set_title(f'Ampt New Coalescence String Melting {energy}GeV')
        for cent in cents:
            df_cent = df_energy[df_energy['cent'] == cent]
            df_cent = df_cent.sort_values(by=['n'])

            if cent > 0:
                ax_energy.errorbar(df_cent['n'], df_cent['sig'], yerr=df_cent['sig_err'], ls='none', alpha=0.7,
                                   marker='o', label=cent_map[cent])

            fig, ax = plt.subplots(dpi=144)
            ax.errorbar(df_cent['n'], df_cent['sig'], yerr=df_cent['sig_err'], ls='none', marker='o')
            ax.axhline(1, color='black', ls='--')
            ax.axhline(0, color='black', ls='-')
            ax.grid()
            ax.set_xlabel('Total Protons in Event')
            ax.set_ylabel(r'$\sigma(c)^2$')
            ax.set_title(f'Ampt New Coalescence String Melting {energy}GeV {cent_map[cent]} Centrality')
            fig.tight_layout()
            fig.savefig(f'{fig_out_path}/Ampt_New_Coal_{energy}GeV_{cent}_cent.png')
        ax_energy.legend()
        fig_energy.tight_layout()
        fig_energy.savefig(f'{fig_out_path}/Ampt_New_Coal_{energy}GeV_allcents.png')
        ax_energy.set_ylim(0.75, 1.1)
        fig_energy.savefig(f'{fig_out_path}/Ampt_New_Coal_{energy}GeV_allcents_zoom.png')
    # fig, ax = plt.subplots(dpi=144)
    # ax.axhline(0, color='black', ls='--')
    # ax.set_xlabel('Centrality')
    # ax.set_ylabel(r'1 - $\sigma(c)^2$')
    # ax.set_title('Ampt New Coalescence String Melting')
    # for energy in pd.unique(df['energy']):
    #     df_energy = df[df['energy'] == energy]
    #     ax.errorbar(df_energy['cent'], df_energy['corr'], yerr=df_energy['corr_err'], ls='-', marker='o', alpha=0.8,
    #                 label=f'{energy}GeV')
    # ax.legend()
    # fig.tight_layout()

    # plt.show()


def random_test():
    n_tracks = 10
    n_events = 100000
    n_bootstraps = 100
    bin_upper_a = 1.0 / 3
    bin_upper_b = 2.0 / 3

    tracks = np.random.random((n_events, n_tracks))
    nas = np.sum(tracks < bin_upper_a, axis=1)
    nbs = np.sum((bin_upper_a <= tracks) & (tracks < bin_upper_b), axis=1)
    sigma = Measure(*calc_bootstraps(nas, nbs, n_bootstraps)[:-1])
    r = Measure(*calc_bootstraps(nas, nbs, n_bootstraps, func=calc_pearson_r)[:-1])
    print(f'sigma = {sigma}')
    print(f'r = {r}')


def read_file(file_path, read_branches, cents, ref3_edges, bin_width, bin_sep, pid, max_rapid, min_pt, max_pt, phi_bins,
              seed):
    rng = np.random.default_rng(seed)
    vector.register_awkward()
    proton_mass = 0.09382721  # GeV
    nas, nbs = {cent: {} for cent in cents}, {cent: {} for cent in cents}
    phi_hists = {cent: {} for cent in cents}
    with uproot.open(file_path) as file:
        tracks_all = file['tree'].arrays(read_branches)
        for cent in cents:
            tracks = tracks_all[(tracks_all.refmult3 <= ref3_edges[cent][0]) &
                                (tracks_all.refmult3 > ref3_edges[cent][1])]
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                            with_name='Momentum3D')
            protons = tracks[(tracks['pid'] == pid) & (tracks.pt > min_pt) & (tracks.pt < max_pt)]
            proton_rapid = rapidity(protons.px, protons.py, protons.pz, protons.px * 0 + proton_mass)
            protons = protons[abs(proton_rapid) < max_rapid]
            n_protons = ak.count(protons.pid, axis=1)
            for n_proton in np.unique(n_protons):
                protons_np = protons[n_protons == n_proton]
                n_protons_a, n_protons_b, phi_hist = bin_events(protons_np, bin_width, bin_sep, phi_bins, rng)
                nas[cent].update({n_proton: list(n_protons_a)})
                nbs[cent].update({n_proton: list(n_protons_b)})
                phi_hists[cent].update({n_proton: phi_hist})

    return nas, nbs, phi_hists


def bin_events(tracks, bin_width, bin_sep, phi_bins, rng):
    phi_hist = np.histogram(ak.flatten(tracks.phi), bins=phi_bins)[0]
    num_events = len(tracks)
    bins_a = rng.random(num_events) * 2 * np.pi - np.pi
    bins_b = bins_a + bin_sep
    n_protons_a = ak.sum(((tracks.phi > bins_a) & (tracks.phi <= bins_a + bin_width)) |
                         (tracks.phi + 2 * np.pi > bins_a) & (tracks.phi + 2 * np.pi <= bins_a + bin_width),
                         axis=1)
    n_protons_b = ak.sum(((tracks.phi > bins_b) & (tracks.phi <= bins_b + bin_width)) |
                         (tracks.phi + 2 * np.pi > bins_b) & (tracks.phi + 2 * np.pi <= bins_b + bin_width),
                         axis=1)

    return n_protons_a, n_protons_b, phi_hist


def bin_events_old(tracks, bin_width, pid, max_eta, min_pt, max_pt, phi_bins, rng):
    tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                    with_name='Momentum3D')
    tracks = tracks[(tracks['pid'] == pid) & (abs(tracks.eta) < max_eta) & (tracks.pt > min_pt) & (tracks.pt < max_pt)]
    phi_hist = np.histogram(ak.flatten(tracks.phi), bins=phi_bins)[0]
    num_events = len(tracks)
    bins_a = rng.random(num_events) * 2 * np.pi - np.pi
    # bins_b = bins_a + np.deg2rad(65)
    # bins_b = bins_a + np.deg2rad(65) * np.random.rand(num_events) + np.deg2rad(1)
    bins_b = bins_a + bin_width  # Second bin right next to first
    n_protons_a = ak.sum(((tracks.phi > bins_a) & (tracks.phi <= bins_a + bin_width)) |
                         (tracks.phi + 2 * np.pi > bins_a) & (tracks.phi + 2 * np.pi <= bins_a + bin_width),
                         axis=1)
    n_protons_b = ak.sum(((tracks.phi > bins_b) & (tracks.phi <= bins_b + bin_width)) |
                         (tracks.phi + 2 * np.pi > bins_b) & (tracks.phi + 2 * np.pi <= bins_b + bin_width),
                         axis=1)

    return n_protons_a, n_protons_b, phi_hist


def calc_corr(na, nb):
    return 1 - calc_sigma(na, nb)


def calc_sigma(na, nb):
    daa = np.mean(na ** 2) - np.mean(na) ** 2
    dbb = np.mean(nb ** 2) - np.mean(nb) ** 2
    dab = np.mean(na * nb) - np.mean(na) * np.mean(nb)
    sigma_c_2 = (daa + dbb - 2 * dab) / (np.mean(na + nb))

    return sigma_c_2


def calc_pearson_r(na, nb):
    na_cent = na - np.mean(na)
    nb_cent = nb - np.mean(nb)
    r = np.sum(na_cent * nb_cent) / np.sqrt(np.sum(na_cent ** 2) * np.sum(nb_cent ** 2))

    return r


def calc_bootstraps(na, nb, nbs=50, cent=0, n_proton=0, func=calc_sigma):
    corrs = []
    data_size = len(na)
    for bs_i in range(nbs):
        nabs, nbbs = [], []
        for j in range(data_size):
            index = np.random.randint(0, data_size)
            nabs.append(na[index])
            nbbs.append(nb[index])
        corrs.append(func(np.array(nabs), np.array(nbbs)))

    return func(na, nb), np.std(corrs), cent, n_proton


# def calc_bootstrap(na, nb):
#     data_size = len(na)
#     nabs, nbbs = [], []
#     for j in range(data_size):
#         index = np.random.randint(0, data_size)
#         nabs.append(na[index])
#         nbbs.append(nb[index])
#     return calc_corr(np.array(nabs), np.array(nbbs))


def ampt_str_edges_to_9(ref_str_edges):
    edges = [int(edge) for edge in ref_str_edges]
    edges = sorted(edges, reverse=True)
    edges_9bin = {8: [1000, edges[0]],  # Set upper edge of most central to very high value
                  7: [edges[0], edges[1]]}  # Do 0-5% and 5-10% manually

    cent_bin, edge_index = 6, 1
    while edge_index + 2 < len(edges):
        edges_9bin.update({cent_bin: [edges[edge_index], edges[edge_index + 2]]})
        cent_bin -= 1
        edge_index += 2

    return edges_9bin


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
