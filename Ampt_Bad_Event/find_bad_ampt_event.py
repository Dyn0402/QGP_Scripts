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
import os


def main():
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
            if num_protons > 100:
                print(f'protons: {num_protons}, ref3: {event.refmult3}, event_num: {event.event}')
                event.Show()

    print('donzo')


if __name__ == '__main__':
    main()
