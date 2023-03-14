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


def main():
    size = 25
    raw = np.random.rand(size)
    mix = np.random.rand(size)
    raw_err = np.std(raw)
    mix_err = np.std(mix)

    print(raw)
    print('donzo')


if __name__ == '__main__':
    main()
