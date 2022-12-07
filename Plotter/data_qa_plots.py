#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 05 4:42 PM 2022
Created in PyCharm
Created as QGP_Scripts/data_qa_plots

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import uproot
import matplotlib.pyplot as plt


def main():
    read_trees()
    # read_qas()


def read_trees():
    base_path = 'F:/Research/BES1_Trees/'
    energies = [39]
    hist_name = ''
    track_attributes = ['refmult', 'refmult2', 'refmult3', 'btof_multi', 'btof_match', 'vx', 'vy', 'vz',
                        'proton.dca', 'proton.pt', 'proton.nsigma', 'proton.rapid', 'proton.charge',
                        'proton.nhits_fit', 'proton.m2']
    max_rapid = 0.5
    charge = 1
    max_dca = 1.0
    min_beta = 0.0
    min_nhits_fit = 20.0
    min_m2, max_m2 = 0.6, 1.2
    min_pt_notof, max_pt_notof = 0.4, 0.8
    min_pt_tof, max_pt_tof = 0.8, 2.0
    max_p_notof, max_p_tof = 1.0, 3.0
    fig_ref3_nproton, ax = plt.subplots(dpi=144)
    ref3_nproton = np.histogram2d([], [], bins=(np.arange(-0.5, 800.5, 1), np.arange(-0.5, 100.5, 1)))
    for energy in energies:
        max_nsigma = 2 if energy != 27 else 1
        trees_path = f'{base_path}{energy}GeV/'
        for file_name in os.listdir(trees_path):
            with uproot.open(f'{trees_path}{file_name}') as file:
                tree_index = 1
                while True:
                    tree_name = f'tree;{tree_index}'
                    if tree_name in file:
                        print(file.classnames())
                        print(file[tree_name].arrays(track_attributes))
                        input()
                        # tracks = file[tree_name].arrays(track_attributes)
                        # tracks = tracks[(abs(tracks['proton.rapid']) < max_rapid) &
                        #                 (tracks['proton.charge'] == charge) &
                        #                 (abs(tracks['proton.nsigma']) < max_nsigma) &
                        #                 (tracks['proton.dca'] < max_dca) & (tracks['proton.nhits_fit'] > min_nhits_fit)]
                        # tof_tracks = tracks[(tracks['proton.pt'] >= min_pt_tof) &
                        #                     (tracks['proton.pt'] <= max_pt_tof) &
                        #                     (tracks['proton.pt'] <= max_p_tof)]
                        # tof_tracks = tof_tracks[(tof_tracks['proton.beta'] > min_beta) &
                        #                         (tof_tracks['proton.m2'] > min_m2) & (tof_tracks['proton.m2'] < max_m2)]
                        # tpc_tracks = tracks[(tracks['proton.pt'] >= min_pt_notof) &
                        #                     (tracks['proton.pt'] <= max_pt_notof) &
                        #                     (tracks['proton.p'] <= max_p_notof)]
                        #
                    else:
                        break
    print('donzo')


def read_qas():
    base_path = 'F:/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
    energies = [7]
    hist_name = ''
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
    for energy in energies:
        qa_path = f'{base_path}{energy}GeV/QA_{energy}GeV.root'
        with uproot.open(qa_path) as file:
            print(file.classnames())
            print(file['btof_multi_pile_can;1'].to_hist())
    print('donzo')


if __name__ == '__main__':
    main()
