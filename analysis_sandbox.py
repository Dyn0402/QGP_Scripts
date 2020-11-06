#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 04 4:20 PM 2020
Created in PyCharm
Created as QGP_Scripts/analysis_sandbox.py

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import kurtosis


def main():
    mode = 'save'
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.4, 0.8], 'kurt_lim': [-1.5, -0.2], 'max_sample': 250, 'dpi': 70},
                {'bin_width': 120, 'sd_lim': [0.9, 1.4], 'kurt_lim': [-1.5, -0.2], 'max_sample': 250, 'dpi': 70}]
    for gif in gif_sets:
        make_gif(gif, mode)
    print('donzo')


def make_gif(gif, mode='show'):
    angles = [50.01, 58, 144.55, 222.8, 275, 309.99, 315, 351]

    fig, axs = plt.subplots(3, 1, figsize=(10, 7), dpi=gif['dpi'])
    fig.set_tight_layout(True)

    print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))

    plot_samples = []
    stds = []
    kurts = []

    anim = FuncAnimation(fig, update, frames=np.arange(1, gif['max_sample']), interval=50, repeat_delay=1000,
                         fargs=(axs, angles, gif['bin_width'], plot_samples, stds, kurts,
                                gif['sd_lim'], gif['kurt_lim'], gif['max_sample']))
    if mode == 'save':
        anim.save(f'/home/dylan/Research/Results/Presentations/11-6-20/{gif["bin_width"]}_degree.gif', dpi=gif['dpi'],
                  writer='imagemagick')
    else:
        plt.show()


def update(i, axs, angles, bin_width, plot_samples, stds, kurts, sd_lim, kurt_lim, max_sample):
    samples = i
    print(f'{samples} Samples')
    axs[0].clear()
    axs[0].set_title(f'{bin_width}° Bin  |  {samples} Samples  <-->  {360/samples:.2f}° Steps')
    axs[0].set_xlabel('Particles in Bin')
    axs[0].set_ylabel('Probability Density')
    axs[0].set_ylim(0, 1)
    hist = get_hist(angles, bin_width, samples)
    axs[0].hist(hist, range(-1, len(angles)+1), density=True)

    plot_samples.append(i)
    stds.append(np.std(hist))
    kurts.append(kurtosis(hist))

    axs[1].clear()
    axs[1].set_xlim(0, max_sample)
    axs[1].set_ylim(sd_lim[0], sd_lim[1])
    axs[1].set_ylabel('Standard Deviation')
    axs[1].plot(plot_samples, stds)

    axs[2].clear()
    axs[2].set_xlim(0, max_sample)
    axs[2].set_ylim(kurt_lim[0], kurt_lim[1])
    axs[2].set_ylabel('Kurtosis')
    axs[2].set_xlabel('Number of Samples')
    axs[2].plot(plot_samples, kurts, color='red')

    return axs


def get_hist(angles, bin_width, samples):
    bin_step = 360 / samples

    bin_content = []

    for start in np.arange(0, 360, bin_step):
        end = start + bin_width
        count = 0
        for angle in angles:
            if start <= angle < end:
                count += 1
        if end > 360:
            for angle in angles:
                if start <= angle + 360 < end:
                    count += 1
        bin_content.append(count)

    return bin_content


if __name__ == '__main__':
    main()
