#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 20 4:30 PM 2023
Created in PyCharm
Created as QGP_Scripts/main

@author: Dylan Neff, Dylan
"""
import os

import uproot
import awkward as ak
import vector


def main():
    corr_test()
    print('donzo')


def corr_test():
    base_path = 'F:/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    energy = 7
    max_eta = 0.5
    vector.register_awkward()

    file_dir = f'{base_path}{energy}GeV/'
    files_paths = os.listdir(file_dir)
    for file_path in files_paths:
        print(file_path)
        # return
        with uproot.open(f'{file_dir}{file_path}') as file:
            tracks = file['tree'].arrays(['pid', 'px', 'py', 'pz'])
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                            with_name='Momentum3D')
            tracks = tracks[abs(tracks.eta) < max_eta]
            print(tracks)
        return


if __name__ == '__main__':
    main()
