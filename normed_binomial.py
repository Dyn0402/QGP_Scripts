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
    tracks = [[50, 'r'], [20, 'b'], [10, 'm'], [3, 'g']]
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for track in tracks:
        x = np.asarray(range(0, track[0]+2))
        x2 = x / x[-1]  # np.linspace(0, 1, len(x))
        x3 = x - x[-1] / float(divs)
        ax1.plot(x, binom.pmf(x, track[0], 1.0/divs), color=track[1], marker='.', label=f'{track[0]} trials')
        ax2.plot(x2, binom.pmf(x, track[0], 1.0 / divs), color=track[1], marker='.', label=f'{track[0]} trials')
        ax3.plot(x3, binom.pmf(x, track[0], 1.0 / divs), color=track[1], marker='.', label=f'{track[0]} trials')

    ax1.title.set_text('Binomial Distributions')
    ax1.set_xlabel('Number of Successes')
    ax1.set_ylabel('Probability')
    ax1.legend(loc='upper right')
    ax2.title.set_text('Compressed Binomial Distributions')
    ax2.set_xlabel('Number of Successes / Number of Trials')
    ax2.set_ylabel('Probability')
    ax2.legend(loc='upper right')
    ax3.title.set_text('Shifted Binomial Distributions')
    ax3.set_xlabel('Number of Successes - Number of Trials / Divisions')
    ax3.set_ylabel('Probability')
    ax3.legend(loc='upper right')
    plt.show()


if __name__ == '__main__':
    main()
