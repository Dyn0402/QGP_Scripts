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
    binoms_tracks()
    binoms_divs()
    binoms_divs_stitch()


def binoms_tracks():
    divs = 3
    tracks = [[50, 'r'], [20, 'b'], [10, 'm'], [3, 'g']]
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for track in tracks:
        x = np.asarray(range(0, track[0]+2))
        x2 = x / x[-1]  # np.linspace(0, 1, len(x))
        x3 = x - x[-1] / float(divs)
        ax1.plot(x, binom.pmf(x, track[0], 1.0 / divs), color=track[1], marker='.', label=f'{track[0]} trials')
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


def binoms_divs():
    tracks = 20
    divs = [[60, 'b'], [180, 'g'], [300, 'r']]
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for div in divs:
        x = np.asarray(range(0, tracks + 1))
        x2 = x / x[-1]  # np.linspace(0, 1, len(x))
        x3 = x - x[-1] * float(div[0]) / 360.0
        ax1.plot(x, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')
        ax2.plot(x2, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')
        ax3.plot(x3, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')

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


def binoms_divs_stitch():
    tracks = 20
    divs = [[60, 'b'], [180, 'g'], [300, 'r']]
    fig1, ax1 = plt.subplots()
    # fig2, ax2 = plt.subplots()
    # fig3, ax3 = plt.subplots()
    x = np.asarray(range(0, tracks + 1))
    print(x)
    x_low = x - tracks * 300.0 / 360.0 + 10
    y_low = binom.pmf(x, tracks, 300.0 / 360.0)
    x_high = x - tracks * 60.0 / 360.0 + 10
    y_high = binom.pmf(x, tracks, 60.0 / 360.0)
    norm = (1 - binom.cdf(60.0 / 360.0 * tracks, tracks, 60.0 / 360.0)) + binom.cdf(300.0 / 360.0 * tracks, tracks, 300.0 / 360.0)
    ax1.plot(x_low / tracks, y_low / norm, marker='.', label='300')
    ax1.plot(x_high / tracks, y_high / norm, marker='.', label='60')
    ax1.plot(x / tracks, binom.pmf(x, tracks, 180.0 / 360.0), marker='.', label='180')
    ax1.legend()
    # y = [binom.pmf(xi + tracks * 300.0 / 360.0 - tracks / 2, tracks, 300.0 / 360.0) if xi <= 10 else
    #      binom.pmf(xi - tracks * 60.0 / 360.0 + tracks / 2, tracks, 60.0 / 360.0) for xi in x]
    # print(y)
    # y1 = [binom.pmf(xi + tracks * 300.0 / 360.0 - tracks / 2, tracks, 300.0 / 360.0) for xi in x]
    # y2 = [binom.pmf(xi - tracks * 60.0 / 360.0 + tracks / 2, tracks, 60.0 / 360.0) for xi in x]
    # xa = [xi + tracks * 300.0 / 360.0 - tracks / 2 for xi in x]
    # print(y1)
    # print(y2)
    # print(xa)
    # x2 = x / x[-1]
    # ax1.plot(x2, y, marker='.')
    # ax1.plot(x2, y1, marker='.')
    # ax1.plot(x2, y2, marker='.')
    # for div in divs:
    #
    #     x2 = x / x[-1]  # np.linspace(0, 1, len(x))
    #     x3 = x - x[-1] * float(div[0]) / 360.0
    #     ax1.plot(x, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')
    #     ax2.plot(x2, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')
    #     ax3.plot(x3, binom.pmf(x, tracks, div[0] / 360.0), color=div[1], marker='.', label=f'{div[0]}$^\circ$')
    #
    # ax1.title.set_text('Binomial Distributions')
    # ax1.set_xlabel('Number of Successes')
    # ax1.set_ylabel('Probability')
    # ax1.legend(loc='upper right')
    # ax2.title.set_text('Compressed Binomial Distributions')
    # ax2.set_xlabel('Number of Successes / Number of Trials')
    # ax2.set_ylabel('Probability')
    # ax2.legend(loc='upper right')
    # ax3.title.set_text('Shifted Binomial Distributions')
    # ax3.set_xlabel('Number of Successes - Number of Trials / Divisions')
    # ax3.set_ylabel('Probability')
    # ax3.legend(loc='upper right')
    plt.show()


if __name__ == '__main__':
    main()
