#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 9:48 PM 2021
Created in PyCharm
Created as QGP_Scripts/cluster_simulation.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import poisson
import timeit

import awkward as ak


def main():
    n_events = 10
    for event in n_events:
        pass
        sds = []
        skews = []
        kurts = []
        for event in events:
            hist = get_hist(event, bin_width, samples)
            stats = DistStats(hist)
            kurts.append(stats.get_kurtosis().val)
            skews.append(stats.get_skewness().val)
            sds.append(stats.get_sd().val)

        kurt = np.mean(kurts)
        skew = np.mean(skews)
        sd = np.mean(sds)
    print('donzo')


def gen_events(num_events, num_tracks):
    events = []
    poisson.rvs(10, size=num_events)
    for event in range(num_events):

        event = gen_event(num_tracks)
        events.append(event)

    return events


def gen_event(num_tracks):
    phis = []
    for track in range(num_tracks):
        phis.append(random.random() * 360)

    return phis


if __name__ == '__main__':
    main()
