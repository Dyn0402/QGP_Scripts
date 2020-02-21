#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 14 10:20 PM 2020
Created in PyCharm
Created as QGP_Scripts/se_standard_deviation

@author: Dylan Neff, dylan
"""

import numpy as np


def main():
    data = [10, 11, 13, 13, 14, 15]
    mean = np.mean(data)
    mu2 = sum((data - mean)**2) / len(data)
    mu4 = sum((data - mean)**4) / len(data)
    m4 = mu4 / mu2**2
    sd_err_root = np.sqrt(mu2/(2*len(data)))
    sd_err_delta = np.sqrt((mu4 - mu2**2) / len(data))
    sd_err_delta2 = np.sqrt((m4 - 1) * mu2 / (4 * len(data)))
    sd_err_delta3 = np.sqrt((mu4 / mu2**2 - 1) / (4 * len(data)))

    print(f'Root err: {sd_err_root}  |  Delta err: {sd_err_delta}  |  Delta err2: {sd_err_delta2}  |  Delta err3: {sd_err_delta3}')


if __name__ == '__main__':
    main()
