#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 27 10:15 PM 2021
Created in PyCharm
Created as QGP_Scripts/tree_read_speed_analysis.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def main():
    path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Other/Ampt_Tree_Speed_Tests/'
    platforms = {'ubuntu': {'file': path+'ubuntu_speed.txt'}, 'windows': {'file': path+'windows_speed.txt'}}

    for platform in platforms:
        platforms[platform].update({'set_times': read_times(platforms[platform]['file'])})
        platforms[platform].update({'90_10_times': get_times(platforms[platform]['set_times'])})

    plot_times(platforms)
    plot_time_ratios(platforms)


def read_times(path):
    set_times = {}
    set_num = 0
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'Reading 7GeV trees.' in line:
                set_num += 1
                set_times.update({set_num: {}})
            if '% complete | time: ' in line:
                time = line.strip().split()[-4].strip('s')
                percent = line[:line.find('%')].split()[-1]
                print(f'set_num: {set_num}, percent: {percent}, time: {time}')
                set_times[set_num].update({int(percent): float(time)})

    return set_times


def get_times(set_times, upper_percent=90, lower_percent=10):
    times = {'sets': [], 'time': []}
    for sets, set_data in set_times.items():
        times['sets'].append(sets)
        times['time'].append(set_data[upper_percent] - set_data[lower_percent])

    return times


def plot_times(platforms):
    fig, ax = plt.subplots()
    for platform, platform_data in platforms.items():
        x = platform_data['90_10_times']['sets']
        y = platform_data['90_10_times']['time']
        popt, pcov = curve_fit(lin, x, y)
        ax.scatter(x, y, label=platform)
        ax.plot(x, lin(np.array(x), *popt), label=f'y = {popt[0]:.2f} * x + {popt[1]:.2f}')
    ax.set_xlabel('Number of Sets Run Together')
    ax.set_ylabel('Time (s) from 10% to 90% Complete')
    ax.legend()
    ax.grid()
    fig.tight_layout()
    plt.show()


def plot_time_ratios(platforms):
    fig, ax = plt.subplots()
    x = platforms['windows']['90_10_times']['sets']
    y = []
    for u, w in zip(platforms['ubuntu']['90_10_times']['time'], platforms['windows']['90_10_times']['time']):
        y.append(w / u)
    ax.scatter(x, y)
    ax.set_xlabel('Number of Sets Run Together')
    ax.set_ylabel('Ubuntu/Windows Time (s) from 10% to 90% Complete')
    ax.grid()
    fig.tight_layout()
    plt.show()


def lin(x, a, b):
    return a * x + b


if __name__ == '__main__':
    main()
