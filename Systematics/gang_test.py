#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 14 8:23 PM 2021
Created in PyCharm
Created as QGP_Scripts/gang_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from Measure import Measure


def main():
    # sigy = 0.8
    dylan_test()
    gang_test()


    print('donzo')


def dylan_test():
    sigy_list = np.arange(0, 5.0, 0.5)
    reps = 100
    m = 0.0012
    b = 4.7
    d_min, d_max = 0, 5
    def_d = 4
    var_d_low = 3
    var_d_high = 5
    n_points = 1000
    plot = False
    n_dsamples = 50

    sigy_plot = []
    sys_plot = []
    sys_mean = []
    for sigy in sigy_list:
        sys_sig = []
        for rep in range(reps):
            d_points = np.random.uniform(d_min, d_max, n_points)
            y_points = [np.random.normal(ytrue_d(d, m, b), sigy) for d in d_points]

            if plot:
                plot_data(d_points, y_points, ytrue_d(d_points, m, b))

            y0 = yd_cut_mean(y_points, d_points, def_d)
            yis = []
            for i in range(n_dsamples):
                yis.append(yd_cut_mean(y_points, d_points, np.random.uniform(var_d_low, var_d_high)).val)
            sys = np.std(yis)
            sigy_plot.append(sigy)
            sys_plot.append(sys)
            sys_sig.append(sys)
        sys_mean.append(get_mean(sys_sig))

    print(f'Guess: {(def_d - var_d_low) * m / 2}')
    # print(f'Sys: {sys / np.sqrt(12)}')

    fig, ax = plt.subplots()
    ax.set_title('Dylan Test')
    ax.plot(sigy_list, [x.val for x in sys_mean], color='r', marker='.')
    ax.fill_between(sigy_list, [x.val - x.err for x in sys_mean], [x.val + x.err for x in sys_mean],
                    color='r', alpha=0.3)
    ax.scatter(sigy_plot, sys_plot, alpha=0.3)
    plt.grid()
    plt.show()


def gang_test():
    sigy_list = np.arange(0, 5.0, 0.5)
    reps = 100
    m = 0.0012
    b = 4.7
    d_min, d_max = 0, 5
    def_d = 4
    var_d = 3
    n_points = 1000
    plot = False

    sigy_plot = []
    sys_plot = []
    sys_mean = []
    for sigy in sigy_list:
        sys_sig = []
        for rep in range(reps):
            d_points = np.random.uniform(d_min, d_max, n_points)
            y_points = [np.random.normal(ytrue_d(d, m, b), sigy) for d in d_points]

            if plot:
                plot_data(d_points, y_points, ytrue_d(d_points, m, b))

            y0 = yd_cut_mean(y_points, d_points, def_d)
            yi = yd_cut_mean(y_points, d_points, var_d)
            sys = calc_sys(y0, yi)
            sigy_plot.append(sigy)
            sys_plot.append(sys)
            sys_sig.append(sys)
        sys_mean.append(get_mean(sys_sig))

    print(f'Guess: {(def_d - var_d) * m / 2}')
    # print(f'Sys: {sys / np.sqrt(12)}')

    fig, ax = plt.subplots()
    ax.set_title('Gang Test')
    ax.plot(sigy_list, [x.val for x in sys_mean], color='r', marker='.')
    ax.fill_between(sigy_list, [x.val - x.err for x in sys_mean], [x.val + x.err for x in sys_mean],
                    color='r', alpha=0.3)
    ax.scatter(sigy_plot, sys_plot, alpha=0.3)
    plt.grid()
    plt.show()


def plot_data(d_points, y_points, ytrue):
    fig1, ax1 = plt.subplots()
    ax1.hist(d_points)
    ax1.set_xlabel('d')

    fig2, ax2 = plt.subplots()
    ax2.scatter(d_points, y_points, alpha=0.5, label='y sample')
    ax2.plot(d_points, ytrue, color='red', ls='--', label='y true')
    ax2.set_xlabel('d')
    ax2.set_ylabel('y')
    ax2.grid()
    ax2.legend()
    plt.show()


def ytrue_d(d, m=0.12, b=4.7):
    return d * m + b


def yd_cut_mean(ys, ds, d_max):
    ys, ds = np.array([(y, d) for y, d in zip(ys, ds) if d < d_max]).T

    return get_mean(ys)


def get_mean(xs):
    return Measure(np.mean(xs), np.std(xs) / np.sqrt(len(xs)))


def calc_sys(y0, yi):
    a = abs(y0.val - yi.val)
    b = np.sqrt(abs(yi.err**2 - y0.err**2))
    if a > b:
        sys = np.sqrt(a**2 - b**2) / np.sqrt(12)
    else:
        sys = 0
    # print(f'y0: {y0}, yi: {yi}, a: {a}, b: {b}, sys: {sys}')

    return sys


if __name__ == '__main__':
    main()
