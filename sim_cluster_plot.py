#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 20 1:08 PM 2020
Created in PyCharm
Created as QGP_Scripts/sim_cluster_plot.py

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt
import numpy as np


def main():
    mu = 3.0
    sigma = [0.02, 0.05, 0.1]
    a = 1.0
    x = np.linspace(0, 2*np.pi, 10000)
    for sig in sigma:
        y = gaus(x, mu, sig, a)
        plt.plot(x, y, label=f'sigma {sig}')
    plt.legend()
    plt.show()
    print('donzo')


def gaus(x, mu, sigma, a):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2))


if __name__ == '__main__':
    main()
