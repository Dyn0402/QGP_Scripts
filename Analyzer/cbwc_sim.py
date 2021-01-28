#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 28 12:27 PM 2021
Created in PyCharm
Created as QGP_Scripts/cbwc_sim.py

@author: Dylan Neff, dylan
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy


def main():
    ref3_vals = range(200, 300)
    ref3_events = np.asarray([int(np.exp(-5 * ref3)) for ref3 in ref3_vals])
    binom_n = np.asarray([20 for ref3 in ref3_vals])

    plot_input_dists(ref3_vals, ref3_events, binom_n)
    simulate(ref3_vals, ref3_events, binom_n)
    calc_moments()
    plot_moments()

    print('donzo')


def plot_input_dists(ref3_vals, ref3_events, binom_n):
    fig1, ax1 = plt.subplots()
    ax1.plot(ref3_vals, ref3_events)

    fig2, ax2 = plt.subplots()



if __name__ == '__main__':
    main()
