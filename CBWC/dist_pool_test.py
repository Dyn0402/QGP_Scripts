#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 22 5:28 PM 2021
Created in PyCharm
Created as QGP_Scripts/dist_pool_test

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
from scipy.stats import binom, poisson, skellam


def main():
    dist = skellam(20, 1.2)
    events = [{'events': 5, 'rndm': np.random.randint(0, 2147483647)}, {'events': 8, 'rndm': np.random.randint(0, 2147483647)}]
    func = partial(test_func, dist)
    with Pool(2) as p:
        res = p.map(func, events)

    print(res)
    print('donzo')


def test_func(dist, events):
    np.random.seed(seed=events['rndm'])
    return np.random.random(events['events'])
    # return dist.rvs(size=events)


if __name__ == '__main__':
    main()
