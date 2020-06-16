#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 16 10:40 AM 2020
Created in PyCharm
Created as QGP_Scripts/gen_function_Nihal.py

@author: Dylan Neff, dylan
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


def main():
    x = range(0, 25)
    dist = 400*(2.3*norm.pdf(x, 9.3, 2.7) + norm.pdf(x, 16.4, 1.22)) + 3
    print(dist)
    out = ''
    for x, y in zip(x, dist):
        out += f'{{{x},{int(y)}}}, '
    print(out)
    plt.plot(x, dist)
    plt.show()
    print('donzo')


if __name__ == '__main__':
    main()
