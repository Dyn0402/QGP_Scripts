#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 03 1:38 PM 2021
Created in PyCharm
Created as QGP_Scripts/sub_event_analysis.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import uproot
import awkward as ak
import vector


def main():
    path = '/home/dylan/Research/BES1_Trees/7GeV/'
    tree_name = 'tree'
    event_branches = ['refmult3']
    particle = 'proton'
    track_branches = ['phi']
    branch_names = event_branches + [f'{particle}.{branch}' for branch in track_branches]
    files = 1

    file_num = 0
    for file_name in os.listdir(path):
        print(file_name)
        if '.root' not in file_name:
            continue
        root_path = os.path.join(path, file_name)
        with uproot.open(root_path) as file:
            branches = file[tree_name].arrays(branch_names)
            # Cuts
            # Phi selection
            # phi = branches['proton.phi']
            # Add to mixer


            # Azimuthal Binning
            bin_count = np.sum((0.2 < branches['proton.phi']) & (branches['proton.phi'] < 1.2), axis=1)

            print((0.2 < branches['proton.phi']) & (branches['proton.phi'] < 1.2))
            print(bin_count)

        file_num += 1
        if file_num >= files:
            break


if __name__ == '__main__':
    main()
