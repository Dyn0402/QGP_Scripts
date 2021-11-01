
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 29 2:33 PM 2021
Created in PyCharm
Created as QGP_Scripts/anti_cluster

@author: Dylan Neff, Dyn04
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def main():
    mean = 1.5
    sd = 0.002
    dist = norm(mean, sd)
    dist_high = norm(mean + 2 * np.pi, sd)
    dist_low = norm(mean - 2 * np.pi, sd)
    m = dist.pdf(mean)
    print(dist.pdf(mean))
    r = np.linspace(-5, 10, 10000)
    p = r * 0 + 1.0 / r[-1]

    wrap1_dist = dist.pdf(r) + dist_low.pdf(r) + dist_high.pdf(r)
    m1 = dist.pdf(mean) + dist_low.pdf(mean) + dist_high.pdf(mean)
    # p1 = 1 - (dist.pdf(r) + dist_low.pdf(r) + dist_high.pdf(r)) / m
    p1 = 1 - wrap1_dist / m1
    # p2 = 1
    p = p1 * p

    fig1, ax1 = plt.subplots()
    ax1.plot(r, p)
    ax1.axhline(0, color='r', ls='--')
    ax1.axvline(0, ls='--')
    ax1.axvline(2 * np.pi, ls='--')

    fig2, ax2 = plt.subplots()
    ax2.plot(r, p)
    ax2.axhline(0, color='r', ls='--')
    ax2.set_xlim([0, 2 * np.pi])

    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
