#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 20 10:20 PM 2023
Created in PyCharm
Created as QGP_Scripts/poisson_moments.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    means = np.linspace(1, 30, 1000)
    vars = means
    skews = 1 / np.sqrt(means)
    kurts = 1 / means

    c3 = skews * np.power(vars, 1.5)
    c4 = kurts * vars * vars

    plt.plot(means, c4, label='C4')
    plt.plot(means, c3, label='C3')
    plt.xlabel('C1')
    plt.legend()
    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
