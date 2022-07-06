#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 28 8:01 PM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_ref3_compare.py

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import uproot

from multiprocessing import Pool
import tqdm
import istarmap


def main():
    base_path = 'D:/Transfer/Research/'
    ampt_version_paths = {
        'Baryon_First': 'AMPT_Trees_Baryon_First',
        'New_Coalescence': 'AMPT_Trees',
        'Meson_First': 'AMPT_Trees_Meson_First'
    }
    mid_path = '/min_bias/string_melting/'
    energies = [7]
    tree_name = 'tree'
    tree_attribute = 'refmult3'
    bins = np.arange(-0.5, 500.5, 1)
    threads = 6

    ampt_version_dists = {}

    for energy in energies:
        fig, ax = plt.subplots()
        ax.set_title(f'{energy}GeV refmult3')
        fig.canvas.manager.set_window_title(f'{energy}GeV refmult3')
        for v_name, v_path in ampt_version_paths.items():
            print(f'Starting {v_name}')
            ref3_dist = np.array([])
            path = f'{base_path}{v_path}{mid_path}{energy}GeV/'
            jobs = [(f'{path}{root_name}', tree_name, tree_attribute) for root_name in os.listdir(path)]
            with Pool(threads) as pool:
                for file_ref3 in tqdm.tqdm(pool.istarmap(read_file_ref3, jobs), total=len(jobs)):
                    ref3_dist = np.append(ref3_dist, file_ref3)
            # for root_name in os.listdir(path):
            #     print(root_name)
            #     with uproot.open(path + root_name) as file:
            #         file_ref3 = file[tree_name].arrays(tree_attribute, library='np')[tree_attribute]
            #         ref3_dist = np.append(ref3_dist, file_ref3)
            ampt_version_dists.update({v_name: np.histogram(ref3_dist, bins=bins)})
            ax.hist(ref3_dist, bins=bins, density=True, histtype='step', alpha=0.8, label=v_name)
        ax.legend()
        fig.tight_layout()

    # fig, ax = plt.subplots()
    # ampt_version_dists['Baryon_First'] / ampt_version_dists['Meson_First']
    # ampt_version_dists['New_Coalescence'] / ampt_version_dists['Meson_First']

    plt.show()
    print('donzo')


def read_file_ref3(root_path, tree_name, tree_attribute):
    with uproot.open(root_path) as file:
        file_ref3 = file[tree_name].arrays(tree_attribute, library='np')[tree_attribute]
    return file_ref3


if __name__ == '__main__':
    main()
