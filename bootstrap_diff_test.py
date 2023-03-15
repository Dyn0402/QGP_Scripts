#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 13 11:18 PM 2023
Created in PyCharm
Created as QGP_Scripts/bootstrap_diff_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from Measure import Measure


def main():
    size = 25
    raw = norm(5, 0.7).rvs(size)
    mix = norm(3, 0.2).rvs(size)
    raw_err = np.std(raw)
    mix_err = np.std(mix)
    diff_err = np.std([r - m for r in raw for m in mix])
    print('Raw - Mix Err: ', np.sqrt(raw_err**2 + mix_err**2))
    print('Diff Err: ', diff_err)

    div_err = np.std([r / m for r in raw for m in mix])
    print('Raw / Mix Err: ', (Measure(5, raw_err) / Measure(3, mix_err)).err)
    print('Div Err: ', div_err)

    # print(raw)
    print('donzo')


if __name__ == '__main__':
    main()
