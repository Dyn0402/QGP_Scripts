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
    c = 1.0
    k = np.arange(0, n, 1)
    mod_b = []
    mod_p_plot = []
    for ki in k:
        mod_b.append(mod_binom(ki, n, p, c))
        mod_p_plot.append(mod_p(ki, n, p, c))
        print(f'k: {ki}  mod: {mod_b[-1]}')
    mod_b = np.asarray(mod_b) / sum(mod_b)
    plt.plot(k, mod_p_plot)
    plt.show()
    plt.scatter(k, stats.binom.pmf(k, n, p), color='red', label='Binomial')
    plt.scatter(k, mod_b, color='blue', label='Mod_Binomial')
    plt.legend()
    plt.show()
    plt.scatter(k, mod_b - stats.binom.pmf(k, n, p))
    plt.axhline(0, ls='--', color='red')
    plt.show()
    print('donzo')


def mod_binom(k, n, p, c):
    comb_factor = math.factorial(n) / (math.factorial(k) * math.factorial(n - k))
    p = mod_p(k, n, p, c)
    q = 1 - p

    return comb_factor * p**k * q**(n-k)


def mod_p(k, n, p, c):
    k = np.asanyarray(k)
    if k >= n*p:
        p_mod = p + (k - n*p) / (n - n*p) * c / n
    else:
        p_mod = p + (n*p - k) / (n*p) * c / n

    return p_mod


if __name__ == '__main__':
    main()
