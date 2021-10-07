#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 07 12:38 PM 2021
Created in PyCharm
Created as QGP_Scripts/random_timing

@author: Dylan Neff, dylan
"""


import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
import timeit


def main():

    setup = 'from __main__ import rand_indiv\n' \
            'from __main__ import rand_set\n' \
            'from __main__ import set_events\n'\
            'import numpy as np\n'\
            'import random\n'
    individual = 'rand_indiv(set_events())'
    sets = 'rand_set(set_events())'
    print(sum(timeit.repeat(setup=setup, stmt=individual, repeat=5, number=1)))
    print(sum(timeit.repeat(setup=setup, stmt=sets, repeat=5, number=1)))
    print('donzo')


def set_events():
    seed = 20
    event_num = 100000
    mean = 10
    rng = np.random.default_rng(seed)
    event_dist = poisson(mean)

    return event_dist.rvs(size=event_num, random_state=rng)


def rand_indiv(events):
    z = []
    for x in events:
        z += [random.random() for y in range(x)]


def rand_set(events):
    [np.random.rand(x) for x in events]


if __name__ == '__main__':
    main()
