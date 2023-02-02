#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 31 6:22 PM 2023
Created in PyCharm
Created as QGP_Scripts/ampt_v2

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import uproot
import awkward as ak
import vector

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from Azimuthal_Correlations.az_corr import get_ampt_ref3_edges
from Sub_Resample.sub_event_resample_algorithm import get_resamples4


def main():
    test_file()
    print('donzo')


def test_file():
    base_path = 'F:/Research/'
    data_path = f'{base_path}AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    ampt_cent_path = f'{base_path}Ampt_Centralities_New_Coalescence/string_melting/'
    energy = 62
    bin_width = np.deg2rad(120)
    samples = 72
    max_eta = 0.5
    min_pt = 0.2  # GeV
    max_pt = 2.0  # GeV
    pid = 2212  # Proton?
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    seed = 44
    threads = 11

    file_dir = f'{data_path}{energy}GeV/'
    file_paths = os.listdir(file_dir)[:100]
    read_branches = ['pid', 'px', 'py', 'pz', 'refmult3']

    ref3_edges = get_ampt_ref3_edges(ampt_cent_path, energy)

    seeds = iter(np.random.SeedSequence(seed).spawn(len(file_paths)))
    jobs = [(f'{file_dir}{path}', read_branches, cents, ref3_edges, bin_width, pid, max_eta, min_pt, max_pt,
             samples, next(seeds)) for path in file_paths]

    data = {cent: {x: np.zeros(x + 1) for x in range(100)} for cent in cents}
    with Pool(threads) as pool:
        for file_data in tqdm.tqdm(pool.istarmap(read_file, jobs), total=len(jobs)):
            for cent in data:
                for nproton in data[cent]:
                    data[cent][nproton] += file_data[cent][nproton]

    plt.bar(np.arange(0, 11, 1), data[8][10], width=1)
    plt.show()


def read_file(file_path, read_branches, cents, ref3_edges, bin_width, pid, max_eta, min_pt, max_pt, samples, seed):
    rng = np.random.default_rng(seed)
    vector.register_awkward()
    data = {cent: {x: np.zeros(x + 1) for x in range(100)} for cent in cents}
    with uproot.open(file_path) as file:
        tracks_all = file['tree'].arrays(read_branches)
        for cent in cents:
            tracks = tracks_all[(tracks_all.refmult3 <= ref3_edges[cent][0]) &
                                (tracks_all.refmult3 > ref3_edges[cent][1])]
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                            with_name='Momentum3D')
            tracks = tracks[
                (tracks['pid'] == pid) & (abs(tracks.eta) < max_eta) & (tracks.pt > min_pt) & (tracks.pt < max_pt)]
            # print(f'\nCentrality {cent}')
            # print(len(tracks))
            # print(tracks[0])
            proton_mults = ak.count(tracks.phi, axis=1)
            # print(proton_mults)
            # print(len(proton_mults))
            all_phis = ak.to_numpy(ak.flatten(tracks.phi)) + np.pi
            # print(all_phis[:10])
            np.random.shuffle(all_phis)
            # print(all_phis[:10])
            mixed_phis = ak.unflatten(all_phis, proton_mults)
            # print(mixed_phis)
            # print(len(mixed_phis))
            for event in mixed_phis:
                # if cent == 8:
                #     print(ak.to_numpy(event))
                data[cent][len(event)] += get_resamples4(ak.to_numpy(event), bin_width, samples, rng)

    return data


if __name__ == '__main__':
    main()
