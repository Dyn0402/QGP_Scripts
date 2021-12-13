#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 10 4:33 PM 2021
Created in PyCharm
Created as QGP_Scripts/resample_alg_timing

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    base_path = '/home/dylan/Research/Results/Resample_Timing_Tests/'
    tracks = ['10_tracks', '30_tracks', '10_tracks_vec', '30_tracks_vec', '10_tracks_alg3', '30_tracks_alg3']
    plot_time_vs_samples(base_path, tracks)
    print('donzo')


def plot_time_vs_samples(base_path, tracks):
    fig, ax = plt.subplots()
    for track in tracks:
        samples = []
        times = []
        time_errs = []
        with open(base_path + f'{track}.txt') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split()
                samples.append(int(line[0]))
                times.append(float(line[1]))
                time_errs.append(float(line[2]))

        times = np.asarray(times)
        time_errs = np.asarray(time_errs)
        ax.fill_between(samples, times + time_errs, times - time_errs, label=f'{track} tracks')
    ax.set_xlabel('Samples')
    ax.set_ylabel('Run Time (s)')
    ax.grid()
    ax.legend(loc='upper left')
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
