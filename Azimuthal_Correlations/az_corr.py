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
    base_path = 'F:/Research/AMPT_Trees/min_bias/string_melting/'
    energy = 7

    energy_dir = f'{base_path}{energy}GeV/'
    files_paths = os.listdir(energy_dir)
    for file_path in files_paths:
        print(file_path)
        return
        # with uproot.open(f'{base_path}{energy}GeV')


if __name__ == '__main__':
    main()
