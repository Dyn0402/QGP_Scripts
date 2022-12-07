#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 06 8:41 PM 2022
Created in PyCharm
Created as QGP_Scripts/TreeReader

@author: Dylan Neff, Dylan
"""

import uproot
import awkward as ak
import vector


class TreeReader:
    def __init__(self, path):
        vector.register_awkward()
        self.path = path
        self.file = uproot.open(path)
        self.tree_name = 'tree'
        self.tracks = []
        self.event_attributes = ['refmult', 'refmult2', 'refmult3', 'btof_multi', 'btof_match', 'vx', 'vy', 'vz',
                                 'proton.dca', 'proton.pt', 'proton.nsigma', 'proton.rapid', 'proton.charge',
                                 'proton.nhits_fit', 'proton.m2']
        self.proton_mass = 0.9  # Fix this

    def read(self):
        tree_index = 1
        while True:
            tree_name = f'{self.tree_name};{tree_index}'
            if tree_name in self.file.classnames():
                tracks = self.file[tree_name].arrays(self.event_attributes)
                tracks =
                # Calculate p and rapidity
                self.tracks.append(tracks)
            else:
                break
            tree_index += 1

    def event_cuts(self):
        pass
