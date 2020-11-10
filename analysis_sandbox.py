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
import random


def main():
    stats_with_events()
    # stats_with_samples()
    # gif_events_maker()
    # gif_maker()
    print('donzo')


def stats_with_samples():
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.4, 0.8], 'kurt_lim': [-1.5, -0.2], 'max_sample': 150, 'dpi': 70,
                 'events': 100, 'tracks': 5}]
    gif = gif_sets[0]
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), dpi=gif['dpi'])
    fig.set_tight_layout(True)

    print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))

    events = gen_events(gif['events'], gif['tracks'])
    plot_samples = []
    stds = []
    kurts = []
    for sample in np.arange(1, gif['max_sample'] + 1):
        plot_samples.append(sample)
        events_hist = hist_events(events, gif['tracks'], gif['bin_width'], sample)
        stds.append(events_hist[1])
        kurts.append(events_hist[2])
    print('Data generated. Animating...')

    n = gif['tracks']
    p = gif['bin_width'] / 360
    q = 1 - p

    axs[0].set_xlim(0, gif['max_sample'])
    axs[0].set_ylabel('Standard Deviation')
    axs[0].plot(plot_samples, stds)
    binom_sd = np.sqrt(n * p * q)
    axs[0].set_ylim(binom_sd * 0.8, binom_sd * 1.2)
    axs[0].axhline(y=binom_sd, color='black', linestyle='--')
    axs[0].plot([1/p], [stds[int(1/p + 0.5) - 1]], marker='o', markersize=8, color='green')

    axs[1].set_xlim(0, gif['max_sample'])
    axs[1].set_ylabel('Kurtosis')
    axs[1].set_xlabel('Number of Samples')
    axs[1].plot(plot_samples, kurts, color='red')
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    axs[1].set_ylim((binom_kurt + 3) * 0.8 - 3, (binom_kurt + 3) * 1.2 - 3)
    axs[1].axhline(y=binom_kurt, color='black', linestyle='--')
    axs[1].plot([1 / p], [kurts[int(1 / p + 0.5) - 1]], marker='o', markersize=8, color='green')
    plt.show()


def stats_with_events():
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 50, 'dpi': 70,
                 'events': 400, 'tracks': 5, 'color': 'red'},
                {'bin_width': 30, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 12, 'dpi': 70,
                 'events': 400, 'tracks': 5, 'color': 'green'}
                ]
    dpi = 70
    gif = gif_sets[0]
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), dpi=dpi)
    fig.set_tight_layout(True)

    events = gen_events(gif['events'], gif['tracks'])
    for gif in gif_sets:
        plot_events = []
        stds = []
        kurts = []
        for event_num in np.arange(1, gif['events'] + 1):
            plot_events.append(event_num)
            events_hist = hist_events(events[:event_num], gif['tracks'], gif['bin_width'], gif['sample'])
            stds.append(events_hist[1])
            kurts.append(events_hist[2])
        print('Data generated. Animating...')

        n = gif['tracks']
        p = gif['bin_width'] / 360
        q = 1 - p

        axs[0].set_xlim(0, gif['events'])
        axs[0].set_ylabel('Standard Deviation')
        axs[0].plot(plot_events, stds, color=gif['color'], label=f'{gif["sample"]} samples')
        binom_sd = np.sqrt(n * p * q)
        axs[0].set_ylim(binom_sd * gif['sd_lim'][0], binom_sd * gif['sd_lim'][1])
        axs[0].axhline(y=binom_sd, color='black', linestyle='--')

        axs[1].set_xlim(0, gif['events'])
        axs[1].set_ylabel('Kurtosis')
        axs[1].set_xlabel('Number of Events')
        axs[1].plot(plot_events, kurts, color=gif['color'], label=f'{gif["sample"]} samples')
        binom_kurt = (1 - 6 * p * q) / (n * p * q)
        axs[1].set_ylim((binom_kurt + 3) * gif['kurt_lim'][0] - 3, (binom_kurt + 3) * gif['kurt_lim'][1] - 3)
        axs[1].axhline(y=binom_kurt, color='black', linestyle='--')
    plt.legend()
    plt.show()


def gif_events_maker():
    mode = 'show'
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.4, 0.8], 'kurt_lim': [-1.5, -0.2], 'max_sample': 150, 'dpi': 70,
                 'events': 1000, 'tracks': 5}]
    for gif in gif_sets:
        make_events_gif(gif, mode)


def gif_maker():
    mode = 'show'
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.4, 0.8], 'kurt_lim': [-1.5, -0.2], 'max_sample': 250, 'dpi': 70},
                {'bin_width': 120, 'sd_lim': [0.9, 1.4], 'kurt_lim': [-1.5, -0.2], 'max_sample': 250, 'dpi': 70}]
    for gif in gif_sets:
        make_gif(gif, mode)


def make_gif(gif, mode='show'):
    angles = [50.01, 58, 144.55, 222.8, 275, 309.99, 315, 351]

    fig, axs = plt.subplots(3, 1, figsize=(10, 7), dpi=gif['dpi'])
    fig.set_tight_layout(True)

    print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))

    print('Generating data...')
    hists = []
    stds = []
    kurts = []
    plot_samples = []
    for sample in range(1, gif['max_sample']+1):
        plot_samples.append(sample)
        hist = get_hist(angles, gif['bin_width'], sample)
        stds.append(np.std(hist))
        kurts.append(kurtosis(hist))
        hists.append(np.histogram(hist, range(-1, len(angles)+1), density=True)[0])
    hist_bins = np.arange(-0.5, len(angles) + 0.5)
    print('Data generated. Animating...')

    anim = FuncAnimation(fig, update, frames=np.arange(1, gif['max_sample']), interval=50, repeat_delay=10000,
                         fargs=(axs, np.asarray(hists), hist_bins, plot_samples, np.asarray(stds), np.asarray(kurts),
                                gif['bin_width'], gif['sd_lim'], gif['kurt_lim'], gif['max_sample'], len(angles)))

    if mode == 'save':
        anim.save(f'/home/dylan/Research/Results/Presentations/11-6-20/{gif["bin_width"]}_degree.gif', dpi=gif['dpi'],
                  writer='imagemagick')
    else:
        plt.show()


def update(i, axs, hists, hist_bins, plot_samples, stds, kurts, bin_width, sd_lim, kurt_lim, max_sample, max_tracks):
    samples = i
    print(f'{samples} Samples')
    axs[0].clear()
    axs[0].set_title(f'{bin_width}° Bin  |  {samples} Samples  <-->  {360/samples:.2f}° Steps')
    axs[0].set_xlabel('Particles in Bin')
    axs[0].set_ylabel('Probability Density')
    axs[0].set_ylim(0, 1)
    # axs[0].hist(hists[samples], range(-1, max_tracks+1), density=True)
    axs[0].bar(hist_bins, hists[samples], width=1, align='center')

    axs[1].clear()
    axs[1].set_xlim(0, max_sample)
    axs[1].set_ylim(sd_lim[0], sd_lim[1])
    axs[1].set_ylabel('Standard Deviation')
    axs[1].plot(plot_samples[:samples], stds[:samples])

    axs[2].clear()
    axs[2].set_xlim(0, max_sample)
    axs[2].set_ylim(kurt_lim[0], kurt_lim[1])
    axs[2].set_ylabel('Kurtosis')
    axs[2].set_xlabel('Number of Samples')
    axs[2].plot(plot_samples[:samples], kurts[:samples], color='red')

    return axs


def make_events_gif(gif, mode):
    fig, axs = plt.subplots(3, 1, figsize=(10, 7), dpi=gif['dpi'])
    fig.set_tight_layout(True)

    print('fig size: {0} DPI, size in inches {1}'.format(fig.get_dpi(), fig.get_size_inches()))

    events = gen_events(gif['events'], gif['tracks'])
    hists = []
    plot_samples = []
    stds = []
    kurts = []
    for sample in np.arange(1, gif['max_sample'] + 1):
        plot_samples.append(sample)
        events_hist = hist_events(events, gif['tracks'], gif['bin_width'], sample)
        stds.append(events_hist[1])
        kurts.append(events_hist[2])
        hists.append(events_hist[0])
    hist_bins = np.arange(0, gif['tracks'] + 1)
    print('Data generated. Animating...')

    anim = FuncAnimation(fig, update_binom, frames=range(gif['max_sample']), interval=50, repeat_delay=10000,
                         fargs=(
                             axs, np.asarray(hists), hist_bins, plot_samples, np.asarray(stds), np.asarray(kurts),
                             gif['bin_width'], gif['sd_lim'], gif['kurt_lim'], gif['max_sample'], gif['tracks']))

    if mode == 'save':
        anim.save(f'/home/dylan/Research/Results/Presentations/11-6-20/{gif["bin_width"]}_degree.gif',
                  dpi=gif['dpi'],
                  writer='imagemagick')
    else:
        plt.show()


def update_binom(i, axs, hists, hist_bins, plot_samples, stds, kurts, bin_width, sd_lim, kurt_lim, max_sample, max_tracks):
    samples = i + 1
    print(f'{samples} Samples')
    axs[0].clear()
    axs[0].set_title(f'{bin_width}° Bin  |  {samples} Samples  <-->  {360/samples:.2f}° Steps')
    axs[0].set_xlabel('Particles in Bin')
    axs[0].set_ylabel('Probability Density')
    axs[0].set_ylim(0, 1)
    # axs[0].hist(hists[samples], range(-1, max_tracks+1), density=True)
    axs[0].bar(hist_bins, hists[i], width=1, align='center')

    n = max_tracks
    p = bin_width / 360
    q = 1 - p

    axs[1].clear()
    axs[1].set_xlim(0, max_sample)
    axs[1].set_ylabel('Standard Deviation')
    axs[1].plot(plot_samples[:i], stds[:i])
    binom_sd = np.sqrt(n * p * q)
    axs[1].set_ylim(binom_sd * 0.8, binom_sd * 1.2)
    axs[1].axhline(y=binom_sd, color='black', linestyle='-')

    axs[2].clear()
    axs[2].set_xlim(0, max_sample)
    axs[2].set_ylabel('Kurtosis')
    axs[2].set_xlabel('Number of Samples')
    axs[2].plot(plot_samples[:i], kurts[:i], color='red')
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    axs[2].set_ylim((binom_kurt + 3) * 0.8 - 3, (binom_kurt + 3) * 1.2 - 3)
    axs[2].axhline(y=binom_kurt, color='black', linestyle='-')

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


def hist_events(events, num_tracks, bin_width, samples):
    hist = []
    for event in events:
        hist.extend(get_hist(event, bin_width, samples))
    sd = np.std(hist)
    kurt = kurtosis(hist)
    hist = np.histogram(hist, np.arange(-0.5, num_tracks + 1.5), density=True)[0]
    print(f'{bin_width}° bins,  {samples} samples,  {len(events)} events,  {num_tracks} tracks,  sd = {sd},  kurt = {kurt}')

    return np.asarray(hist), sd, kurt


def gen_events(num_events, num_tracks):
    events = []
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
