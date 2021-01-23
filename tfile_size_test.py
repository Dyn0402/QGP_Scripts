#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 22 7:41 PM 2021
Created in PyCharm
Created as QGP_Scripts/tfile_size_test.py

@author: Dylan Neff, dylan
"""


import os
import matplotlib.pyplot as plt
import numpy as np


def main():
    path = '/home/dylan/Research/TFile_Size_Test3/'
    file_types_colors = {'array': 'red', 'vector': 'blue', 'junk': 'green'}
    sizes = {ftype: {'nvecs': [], 'nvals': [], 'size': [], 'color': color}
             for ftype, color in file_types_colors.items()}

    sizes = get_file_sizes(path, sizes)
    plot_file_sizes(sizes)
    print('donzo')


def get_file_sizes(path, sizes):
    files = os.listdir(path)
    for file in files:
        for ftype in sizes:
            if ftype in file:
                nvecs = int(file.strip('.root').split('_')[-2])
                nvals = int(file.strip('.root').split('_')[-1])
                sizes[ftype]['nvecs'].append(nvecs)
                sizes[ftype]['nvals'].append(nvals)
                sizes[ftype]['size'].append(os.path.getsize(path+file))

    print(sizes)
    return sizes


def plot_file_sizes(sizes):
    nvals = sizes['array']['nvals'][0]
    float_size = 4  # Bytes

    fix, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for ftype, vals in sizes.items():
        x, y = zip(*sorted(zip(vals['nvecs'], vals['size'])))
        if ftype == 'junk':
            ax2.plot(x, y, color=vals['color'], label=ftype)
        else:
            ax1.plot(x, y, color=vals['color'], label=ftype)
    ax1.plot(x, float_size * nvals * np.asarray(x), color='purple', label='minimum')
    plt.title(f'File size with {nvals} value vectors/arrays stored')
    ax1.set_xlabel(f'Number of {nvals} value vectors/arrays stored')
    ax1.set_ylabel(f'Size of file (bytes)')
    # plt.yscale('log')
    ax1.grid()
    ax2.set_ylabel('Size of junk file (bytes)')
    ax1.legend(loc='lower right')
    ax2.legend(loc='upper left')
    plt.show()


if __name__ == '__main__':
    main()
