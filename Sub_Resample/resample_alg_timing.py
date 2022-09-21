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
    # base_path = '/home/dylan/Research/Results/Resample_Timing_Tests/'
    base_path = 'D:/Transfer/Research/Results/Resample_Timing_Tests/'
    tracks = [{'file_name': '10_tracks_alg2', 'color': 'green', 'lightness': 0.8},
              {'file_name': '30_tracks_alg2', 'color': 'green', 'lightness': 1},
              {'file_name': '10_tracks_alg3', 'color': 'blue', 'lightness': 0.8},
              {'file_name': '30_tracks_alg3', 'color': 'blue', 'lightness': 1},
              {'file_name': '10_tracks_alg4', 'color': 'red', 'lightness': 0.8},
              {'file_name': '30_tracks_alg4', 'color': 'red', 'lightness': 1},
              {'file_name': '10_tracks_alg4_norand', 'color': 'orange', 'lightness': 0.8},
              {'file_name': '30_tracks_alg4_norand', 'color': 'orange', 'lightness': 1},
              {'file_name': '10_tracks_alg4_norand_nosort', 'color': 'salmon', 'lightness': 0.8},
              {'file_name': '30_tracks_alg4_norand_nosort', 'color': 'salmon', 'lightness': 1},
              {'file_name': '10_tracks_alg5', 'color': 'purple', 'lightness': 0.8},
              {'file_name': '30_tracks_alg5', 'color': 'purple', 'lightness': 1},
              ]
    plot_time_vs_samples(base_path, tracks)
    print('donzo')


def plot_time_vs_samples(base_path, tracks):
    fig, ax = plt.subplots()
    for track in tracks:
        samples = []
        times = []
        time_errs = []
        with open(base_path + f'{track["file_name"]}.txt') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split()
                samples.append(int(line[0]))
                times.append(float(line[1]))
                time_errs.append(float(line[2]))

        times = np.asarray(times)
        time_errs = np.asarray(time_errs)
        c = lighten_color(track['color'], track['lightness'])
        ax.fill_between(samples, times + time_errs, times - time_errs, color=c, alpha=0.5,
                        label=f'{track["file_name"]}')
    ax.set_xlabel('Samples')
    ax.set_ylabel('Run Time (s)')
    ax.grid()
    ax.legend(loc='upper left')
    fig.tight_layout()
    plt.show()


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


if __name__ == '__main__':
    main()
