#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 05 4:42 PM 2022
Created in PyCharm
Created as QGP_Scripts/data_qa_plots

@author: Dylan Neff, Dylan
"""
import os

import uproot
import matplotlib.pyplot as plt


def main():
    # read_trees()
    read_qas()


def read_trees():
    base_path = 'F:/Research/BES1_Trees/'
    energies = [7]
    hist_name = ''
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
    for energy in energies:
        trees_path = f'{base_path}{energy}GeV/'
        for file_name in os.listdir(trees_path):
            with uproot.open(f'{trees_path}{file_name}') as file:
                print(file.classnames())
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
