#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 2:17 PM 2021
Created in PyCharm
Created as QGP_Scripts/find_bad_ampt_event

@author: Dylan Neff, dylan
"""

import ROOT
# from ROOT import TFile, TTree
import awkward as ak
import uproot
import vector
import os
from datetime import datetime
from multiprocessing import Pool


def main():
    uproot_finder()
    print('donzo')


def root_finder():
    # Found: AMPT_Trees/min_bias/string_melting/7GeV/data_741821621.root EVENT:1617 - 1618
    path = '/home/dylan/Research/AMPT_Trees/min_bias/string_melting/7GeV/'
    tree_num = 0
    for root_path in os.listdir(path):
        file = ROOT.TFile(path + root_path, 'READ')
        tree_num += 1
        print(f'Tree #{tree_num}: {root_path}')
        for event in file.tree:
            num_protons = 0
            for track_index in range(len(event.pid)):
                if event.pid[track_index] == 2212:  # Is proton
                    vec = ROOT.TVector3(event.px[track_index], event.py[track_index], event.pz[track_index])
                    if -0.5 < vec.PseudoRapidity() < 0.5:
                        num_protons += 1
            if num_protons > 80:
                print(f'protons: {num_protons}, ref3: {event.refmult3}, event_num: {event.event}')
                event.Show()


def uproot_finder():
    """ Looks like this runs just as fast as compiled ROOT code, just slow 2-p comparisons """
    print(f'Start {datetime.now()}\n')
    vector.register_awkward()

    # out_file_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    # path = '/home/dylan/Research/AMPT_Trees/'
    path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/7GeV/'
    # path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/string_melting/7GeV/'
    threads = 16
    tree_name = 'tree'
    max_eta = 1
    count_pid = 2212
    track_attributes = ['event', 'refmult3', 'pid', 'px', 'py', 'pz']

    root_paths = []
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            root_paths.append(root_path)

    n_files = len(root_paths)
    print(f'{n_files} files found to be checked.')

    with Pool(threads) as pool:
        bad_trees = pool.starmap(check_file,
                                 [(root_path, tree_name, track_attributes, max_eta, count_pid, num, n_files) for
                                  num, root_path in enumerate(root_paths)])

    for tree in bad_trees:
        if len(tree) > 0:
            print(tree)

    print(f'\nEnd {datetime.now()}')


def check_file(root_path, tree_name, track_attributes, max_eta, count_pid, root_num=-1, root_num_total=-1):
    print(f'{root_num}/{root_num_total} {root_path} :\t{datetime.now()}')
    bad_events = []
    with uproot.open(root_path) as file:
        tracks_og = file[tree_name].arrays(track_attributes)
        tracks = tracks_og[(tracks_og.refmult3 >= 156) & (tracks_og.refmult3 <= 500)]
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        tracks = tracks[abs(tracks.eta) < max_eta]
        tracks = tracks[tracks['pid'] == count_pid]
        tracks = tracks['pid']
        count = ak.count(tracks, axis=1)
        if ak.sum(count > 80) > 0:
            print(count)
            bad_events.append({'tracks': tracks_og[count > 80], 'file': root_path})

        # print(len(count), count)
        return bad_events


if __name__ == '__main__':
    main()
