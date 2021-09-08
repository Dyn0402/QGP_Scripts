#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 2:17 PM 2021
Created in PyCharm
Created as QGP_Scripts/find_bad_ampt_event

@author: Dylan Neff, dylan
"""

import numpy as np
import ROOT
# from ROOT import TFile, TTree
import os


def main():
    # Found: AMPT_Trees/min_bias/string_melting/7GeV/data_741821621.root EVENT:1617 - 1618
    # path = '/home/dylan/Research/AMPT_Trees/min_bias/string_melting/7GeV/'
    path = '/home/dylan/Research/AMPT_Bad_Event/'
    tree_num = 0
    track_attributes = ['pid', 'px', 'py', 'pz']
    for root_path in os.listdir(path):
        file = ROOT.TFile(path + root_path, 'READ')
        tree_num += 1
        print(f'Tree #{tree_num}: {root_path}')
        for event in file.tree:
            for track_index_i in range(len(event.pid) - 1):
                for track_index_j in range(track_index_i+1, len(event.pid)):
                    track_eqs = [getattr(event, x)[track_index_i] == getattr(event, x)[track_index_j]
                                 for x in track_attributes]
                    if np.prod(track_eqs):
                        print(f'Identical particles in event_num: {event.event}')
                        print(f'\ttrack indices: {track_index_i}, {track_index_j}')
                        event.Show()

    print('donzo')


if __name__ == '__main__':
    main()
