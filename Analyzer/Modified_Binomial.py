#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 08 2:33 PM 2020
Created in PyCharm
Created as QGP_Scripts/Modified_Binomial

@author: Dylan Neff, Dyn04
"""

import matplotlib.pyplot as plt
from scipy import stats
import math
import numpy as np


def main():
    n = 50
    p = 0.2
    c = -0.01
    k = np.arange(0, n, 1)
    mod_b = []
    for ki in k:
        mod_b.append(mod_binom(ki, n, p, c))
    mod_b = np.asarray(mod_b) / sum(mod_b)
    plt.plot(k, p + abs(k - n*p) * c / n)
    plt.show()
    plt.scatter(k, stats.binom.pmf(k, n, p), color='red', label='Binomial')
    plt.scatter(k, mod_b, color='blue', label='Mod_Binomial')
    plt.legend()
    plt.show()
    print('donzo')


def mod_binom(k, n, p, c):
    norm = math.factorial(n) / (math.factorial(k) * math.factorial(n - k))
    p_mod = p + abs(k - n*p) * c / n
    q_mod = 1 - p

    return norm * p_mod**k * q_mod**(n-k)


if __name__ == '__main__':
    main()
