#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 30 9:59 PM 2020
Created in PyCharm
Created as QGP_Scripts/normed_binomial.py

@author: Dylan Neff, dylan
"""

from scipy.stats import binom
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def main():
    divs = 3
    tracks = [[50, 'r'], [25, 'b'], [3, 'g']]
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for track in tracks:
        x = range(0, track[0]+2)
        x2 = np.linspace(0, 1, len(x))
        ax1.plot(x, binom.pmf(x, track[0], 1.0/divs), color=track[1])
        ax2.plot(x2, binom.pmf(x, track[0], 1.0 / divs), color=track[1])
    plt.show()


if __name__ == '__main__':
    main()
