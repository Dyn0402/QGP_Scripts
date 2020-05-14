#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 14 10:47 AM 2020
Created in PyCharm
Created as QGP_Scripts/poisson_vs_binomial

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom, poisson


def main():
    n = 60
    p = 0.5
    x = np.arange(60)
    plt.plot(x, binom.pmf(x, n, p), label='binomial')
    plt.plot(x, poisson.pmf(x, n * p), label='poisson')
    plt.legend()
    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()