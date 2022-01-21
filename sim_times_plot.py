#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 20 3:45 PM 2022
Created in PyCharm
Created as QGP_Scripts/sim_times_plot.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from datetime import datetime
from scipy.optimize import curve_fit as cf


def main():
    path = 'C:/Users/Dylan/Desktop/sim_times.txt'
    expected_num = 493
    start_time = datetime(2022, 1, 18, 14, 0)
    count_time = datetime(2022, 1, 18, 18, 0)
    fit_time = datetime(2022, 1, 20, 2, 0)
    times = []
    nums = [expected_num]
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            day = int(line[6])
            hour, minute = [int(x) for x in line[7].split(':')]
            dt = datetime(2022, 1, day, hour, minute)
            times.append(dt) if dt > count_time else 0
    times = sorted(times)

    # nums.append(nums[-1] - 1)
    # times.append(dt)
    #
    num_left = expected_num - np.arange(len(times))
    x_fit, y_fit = zip(*[(x.timestamp(), y) for x, y in zip(times, num_left) if x > fit_time])

    popt, pcov = cf(lin, x_fit, y_fit)

    fig, ax = plt.subplots()
    formatter = DateFormatter('%a %H:%M')
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_tick_params(rotation=30)
    ax.grid()
    ax.scatter(times, num_left, marker='.')
    x_plt = np.array([fit_time, datetime(2022, 1, 21, 5)])
    print(x_plt)
    ax.plot(x_plt, lin(np.array([x.timestamp() for x in x_plt]), *popt), color='red', ls='--')
    ax.axhline(0, color='black')
    # ax.axvline(start_time, ls='--', color='green')
    ax.set_ylabel('Jobs Remaining')
    plt.show()


def lin(x, a, b):
    return a * x + b


if __name__ == '__main__':
    main()
