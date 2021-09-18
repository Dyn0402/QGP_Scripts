#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 2:17 PM 2021
Created in PyCharm
Created as QGP_Scripts/find_bad_ampt_event

@author: Dylan Neff, dylan
"""

import numpy as np
import uproot
import awkward as ak
import vector
import os


def main():
    macro_name = 'ampt_identical_track_finder_single.cpp'
    macro_path = f'C:/Users/Dylan/Source/Repos/Dyn0402/Root_Macros/src/{macro_name}'




def uproot_finder():
    vector.register_awkward()
    # Found: AMPT_Trees/min_bias/string_melting/7GeV/data_741821621.root EVENT:1617 - 1618
    # path = '/home/dylan/Research/AMPT_Trees/min_bias/string_melting/7GeV/'
    # path = '/home/dylan/Research/AMPT_Bad_Event/'
    # path = 'D:/Research/AMPT_Trees/min_bias/string_melting/7GeV/'

    bad_events = []
    path = 'D:/Research/Ampt_Bad_Event/'
    tree_num = 0
    tree_name = 'tree'
    max_eta = 1
    bad_pids = [313, 111]
    track_attributes = ['pid', 'px', 'py', 'pz']
    for root_path in os.listdir(path)[0:100]:
        with uproot.open(path + root_path) as file:
            tracks = file[tree_name].arrays(track_attributes)
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                                with_name='Momentum3D')
            tracks = tracks[abs(tracks.eta) < max_eta]
            for bad_pid in bad_pids:
                tracks = tracks[abs(tracks['pid']) != bad_pid]
            tree_num += 1
            print(f'Tree #{tree_num}: {root_path}')
            track_a, track_b = ak.unzip(ak.combinations(tracks, 2))
            ident_p = track_a == track_b  # Check if momenta same
            ident_pid = ak.where(track_a['pid'] == track_b['pid'], True, False)
            ident = ident_p * ident_pid
            num_ident = [{'event_num': x, 'num_identical': y} for x, y in enumerate(ak.sum(ident, axis=1)) if y > 0]
            if len(num_ident) > 0:
                bad_events.append({})
                id_tracks = track_a[ident]
                print(num_ident)
                print(np.unique(id_tracks[ak.num(id_tracks) > 0].pid, return_counts=True))

    print('donzo')


if __name__ == '__main__':
    main()
