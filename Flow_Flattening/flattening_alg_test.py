#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 05 4:50 PM 2023
Created in PyCharm
Created as QGP_Scripts/flattening_alg_test

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


# Seems like it doesn't matter if I shift the distribution to -pi to +pi or not, flattening works either way?
def main():
    # bounds = (0, 2 * np.pi)
    n_max = 12
    bounds = (8 * np.pi, 8.1 * np.pi)
    # dist = norm(0.8, 0.1)
    dist = norm(8.05 * np.pi, 0.01)
    x = dist.rvs(100000)
    x = impose_periodic_boundary(x, *bounds)
    fig, ax = plt.subplots()
    ax.hist(x, bins=20)
    a, b = get_coefs(x, n_max)
    x_flat = impose_periodic_boundary(flatten(x, a, b), *bounds)
    fig, ax = plt.subplots()
    ax.hist(x_flat, bins=20)
    ax.set_title('Flat1')
    a2, b2 = get_coefs2(x, n_max, *bounds)
    x_flat2 = impose_periodic_boundary(flatten2(x, a2, b2, *bounds), *bounds)
    fig, ax = plt.subplots()
    ax.hist(x_flat2, bins=20)
    ax.set_title('Flat2')
    plt.show()
    print('donzo')


def get_coefs(dist_vals, n_max):
    a, b = [], []
    for n in range(n_max + 1):
        a.append(np.mean(np.cos(n * dist_vals)))
        b.append(np.mean(np.sin(n * dist_vals)))

    return a, b


def get_coefs2(dist_vals, n_max, x_min, x_max):
    l = (x_max - x_min) / 2
    x_avg = (x_max + x_min) / 2
    a, b = [], []
    for n in range(n_max + 1):
        a.append(np.mean(np.cos(n * np.pi * (dist_vals - x_avg) / l)))
        b.append(np.mean(np.sin(n * np.pi * (dist_vals - x_avg) / l)))

    return a, b


def flatten(dist_vals, a, b):
    assert len(a) == len(b)
    shift_vals = dist_vals.copy()
    for n in range(1, len(a)):
        shift_vals += 2. / n * (a[n] * np.sin(n * dist_vals) - b[n] * np.cos(n * dist_vals))

    return shift_vals


def flatten2(dist_vals, a, b, x_min, x_max):
    l = (x_max - x_min) / 2
    x_avg = (x_max + x_min) / 2
    assert len(a) == len(b)
    dist_vals = np.pi * (dist_vals - x_avg) / l
    shift_vals = dist_vals.copy()
    for n in range(1, len(a)):
        shift_vals += 2. / n * (a[n] * np.sin(n * dist_vals) - b[n] * np.cos(n * dist_vals))
    shift_vals = l * shift_vals / np.pi + x_avg

    return shift_vals


def impose_periodic_boundary(dist_vals, bound_low, bound_high):
    period = bound_high - bound_low
    for i in range(len(dist_vals)):
        while dist_vals[i] >= bound_high:
            dist_vals[i] -= period
        while dist_vals[i] < bound_low:
            dist_vals[i] += period

    return dist_vals

if __name__ == '__main__':
    main()
