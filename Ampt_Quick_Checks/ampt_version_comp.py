#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 28 8:01 PM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_version_comp.py

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import uproot

from multiprocessing import Pool
import tqdm
try:
    import istarmap
except ModuleNotFoundError:
    try:
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Analyzer'))
        import istarmap
    except ModuleNotFoundError:
        print('Can\'t find istarmap!')


def main():
    # base_path = 'D:/Transfer/Research/'
    # out_path = ''
    # ampt_version_paths = {
    #     'Baryon_First': 'AMPT_Trees_Baryon_First',
    #     'New_Coalescence': 'AMPT_Trees',
    #     'Meson_First': 'AMPT_Trees_Meson_First'
    # }
    # mid_path = '/min_bias/string_melting/'
    base_path = '/gpfs01/star/pwg/dneff/data/AMPT/'
    out_dir = '/star/u/dneff/ampt_version_comp/'
    ampt_version_paths = {
        'Baryon_First': 'min_bias_baryon_first',
        'New_Coalescence': 'min_bias',
        'Meson_First': 'min_bias_meson_first'
    }
    mid_path = '/string_melting/'
    energies = [7]
    tree_name = 'tree'
    tree_attribute = 'imp'
    # bins = np.arange(-0.5, 430.5, 1)
    bins = np.linspace(-0.1, 25, 200)
    density = True
    y_log = False
    num_files = None  # None for all
    threads = 8

    ampt_version_dists = {}

    for energy in energies:
        fig, ax = plt.subplots()
        ax.set_title(f'{energy}GeV {tree_attribute}')
        fig.canvas.manager.set_window_title(f'{energy}GeV {tree_attribute}')
        for v_name, v_path in ampt_version_paths.items():
            print(f'Starting {v_name}')
            ref3_dist = np.array([])
            path = f'{base_path}{v_path}{mid_path}{energy}GeV/'
            jobs = [(f'{path}{root_name}', tree_name, tree_attribute) for root_name in os.listdir(path)[:num_files]]
            with Pool(threads) as pool:
                for file_ref3 in tqdm.tqdm(pool.istarmap(read_file_ref3, jobs), total=len(jobs)):
                    ref3_dist = np.append(ref3_dist, file_ref3)
            ampt_version_dists.update({v_name: np.histogram(ref3_dist, bins=bins)})
            ax.hist(ref3_dist, bins=bins, density=density, histtype='step', alpha=0.8, label=v_name)
        if y_log:
            ax.set_yscale('log')
        ax.legend()
        fig.tight_layout()

    fig.savefig(f'{out_dir}ampt_{tree_attribute}_version_comp.png')
    plt.show()
    print('donzo')


def read_file_ref3(root_path, tree_name, tree_attribute):
    with uproot.open(root_path) as file:
        file_ref3 = file[tree_name].arrays(tree_attribute, library='np')[tree_attribute]
    return file_ref3


if __name__ == '__main__':
    main()
