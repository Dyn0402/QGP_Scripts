#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 26 6:46 PM 2023
Created in PyCharm
Created as QGP_Scripts/poisson_dist_plots

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    root_file_path = 'F:/Research/Results/Rand_Tests/root_dist.txt'
    sample_poisson_file_path = 'F:/Research/Results/Rand_Tests/sample_poisson_dist.txt'

    with open(root_file_path) as file:
        root_dist = [float(x) for x in file.read().split()]

    with open(sample_poisson_file_path) as file:
        sample_poisson_dist = [float(x) for x in file.read().split()]

    fig, ax = plt.subplots(dpi=144)
    ax.hist(root_dist, label='root_dist', bins=np.arange(-0.5, 20.5, 1), alpha=0.4)
    ax.hist(sample_poisson_dist, label='sample_poisson_dist', bins=np.arange(-0.5, 20.5, 1), alpha=0.4)
    plt.legend()
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
