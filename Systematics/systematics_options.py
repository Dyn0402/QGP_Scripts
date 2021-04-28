#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 12 2:00 PM 2021
Created in PyCharm
Created as QGP_Scripts/systematics_options.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    dca_sd = 0.2

    dca_points_rand = np.random.uniform(0.5, 1.5, 40)
    plot_dca(dca_points_rand, dca_sd)

    dca_points_rep = np.linspace(0.5, 1.5, 4)
    plot_dca(np.array(list(dca_points_rep)*10), dca_sd)

    m_points_rand = np.random.uniform(0.8, 1.2, 40)

    plot_dca_m(dca_points_rand, m_points_rand)

    m_points_rep = np.linspace(0.8, 1.2, 4)

    dca_m_points = [(x, y) for x in dca_points_rep for y in m_points_rep]*2

    plot_dca_m(*np.array(dca_m_points).T)


def plot_dca(dca_points, dca_sd):
    dca_stat_true = dca_points * 1.1 + 0.3
    dca_stat_rand = [np.random.normal(x, dca_sd) for x in dca_stat_true]
    plt.scatter(dca_points, dca_stat_rand, label='Sampled Moment Mean')
    plt.plot(dca_points, dca_stat_true, 'r', label='True Moment Mean')
    plt.xlabel('dca cut (cm)')
    plt.ylabel('arbitrary moment value')
    plt.title(f'{len(dca_points)} Runs')
    plt.legend()
    plt.show()


def plot_dca_m(dca_points, m_points):
    plt.scatter(dca_points, m_points, alpha=0.5)
    plt.xlabel('dca cut (cm)')
    plt.ylabel('m^2 cut (GeV)')
    plt.title(f'{len(dca_points)} Runs')
    plt.show()


if __name__ == '__main__':
    main()
