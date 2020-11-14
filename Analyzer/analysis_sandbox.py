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
import datetime
from DistStats import DistStats
import math
import seaborn as sns
import warnings


def main():
    # angles = [50.01, 58, 144.55, 222.8, 275, 309.99, 315, 351]
    # print(get_hist(angles, 120, 3))
    # events = gen_events(10, 10)
    # hist_events_ebe(events, 10, 120, 50)

    stat_deviation_with_events()
    # stats_with_events()
    # stats_with_samples()
    # gif_events_maker()
    # gif_maker()
    print('donzo')


def stats_with_samples():
    gif_sets = [{'bin_width': 120, 'sd_lim': [0.8, 1.2], 'kurt_lim': [0.5, 1.5], 'max_sample': 100, 'dpi': 70,
                 'events': 1000, 'tracks': 150}]
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

    axs[0].set_title(f'{gif["bin_width"]}° Bin  {gif["tracks"]} Tracks  {gif["events"]} Events')
    axs[0].set_xlim(0, gif['max_sample'])
    axs[0].set_ylabel('Standard Deviation')
    axs[0].plot(plot_samples, stds)
    binom_sd = np.sqrt(n * p * q)
    axs[0].set_ylim(binom_sd * gif['sd_lim'][0], binom_sd * gif['sd_lim'][1])
    axs[0].axhline(y=binom_sd, color='black', linestyle='--')
    axs[0].plot([1/p], [stds[int(1/p + 0.5) - 1]], marker='o', markersize=8, color='green')

    axs[1].set_xlim(0, gif['max_sample'])
    axs[1].set_ylabel('Kurtosis')
    axs[1].set_xlabel('Number of Samples')
    axs[1].plot(plot_samples, kurts, color='red')
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    axs[1].set_ylim((binom_kurt + 3) * gif['kurt_lim'][0] - 3, (binom_kurt + 3) * gif['kurt_lim'][1] - 3)
    axs[1].axhline(y=binom_kurt, color='black', linestyle='--')
    axs[1].plot([1 / p], [kurts[int(1 / p + 0.5) - 1]], marker='o', markersize=8, color='green')
    plt.show()


def stats_with_events():
    gif_sets = [{'bin_width': 120, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 50, 'dpi': 70,
                 'events': 400, 'tracks': 250, 'color': 'red', 'mode': 'normal'},
                {'bin_width': 120, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 3, 'dpi': 70,
                 'events': 400, 'tracks': 250, 'color': 'green', 'mode': 'normal'},
                # {'bin_width': 120, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 50, 'dpi': 70,
                #  'events': 40, 'tracks': 250, 'color': 'blue', 'mode': 'ebe'},
                {'bin_width': 120, 'sd_lim': [0.7, 1.3], 'kurt_lim': [0.5, 1.5], 'sample': 50, 'dpi': 70,
                 'events': 400, 'tracks': 250, 'color': 'purple', 'mode': 'parallel'}
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
            if gif['mode'] == 'ebe':
                events_hist = hist_events_ebe(events[:event_num], gif['tracks'], gif['bin_width'], gif['sample'])
                stds.append(events_hist[0])
                kurts.append(events_hist[1])
            elif gif['mode'] == 'parallel':
                if event_num > 3:
                    events_hist = hist_events_parallel(events[:event_num], gif['tracks'], gif['bin_width'], gif['sample'])
                    stds.append(events_hist[1])
                    kurts.append(events_hist[2])
                else:
                    plot_events.pop()
            else:
                events_hist = hist_events(events[:event_num], gif['tracks'], gif['bin_width'], gif['sample'])
                stds.append(events_hist[1])
                kurts.append(events_hist[2])
        print('Data generated. Animating...')

        n = gif['tracks']
        p = gif['bin_width'] / 360
        q = 1 - p

        axs[0].set_title(f'{gif["bin_width"]}° Bin  {gif["tracks"]} Tracks  {gif["events"]} Events')
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


def stat_deviation_with_events():
    gif_sets = [{'sample': 3, 'color': 'green', 'mode': 'normal'},
                {'sample': 100,  'color': 'red', 'mode': 'normal'},
                {'sample': 100, 'color': 'purple', 'mode': 'parallel'},
                {'sample': 100, 'color': 'blue', 'mode': 'ebe'},
                {'sample': 50, 'color': 'black', 'mode': 'ebe'}
                ]
    dpi = 70
    alpha = 0.3

    tracks = 15
    bin_width = 120
    n_events = 200
    trials = 100
    parallel_min_event = 2

    fig, axs = plt.subplots(3, 1, figsize=(12, 10), dpi=dpi, sharex=True)
    fig.set_tight_layout(True)

    n = tracks
    p = bin_width / 360
    q = 1 - p
    binom_sd = np.sqrt(n * p * q)
    binom_skew = (q - p) / np.sqrt(n * p * q)
    binom_kurt = (1 - 6 * p * q) / (n * p * q)

    deviations = [{'sd': {}, 'skew': {}, 'kurt': {}} for gif in gif_sets]
    dev_errs = [{'sd': {}, 'skew': {}, 'kurt': {}} for gif in gif_sets]
    for gif_num, gif in enumerate(deviations):
        for stat in gif:
            for event_num in np.arange(1, n_events + 1):
                if gif_sets[gif_num]['mode'] == 'parallel' and event_num < parallel_min_event:
                    continue
                gif[stat].update({event_num: []})
                dev_errs[gif_num][stat].update({event_num: []})

    for trial in range(1, trials + 1):
        print(f'Trial {trial} of {trials}  {datetime.datetime.now().strftime("%H:%M:%S")}')
        events = gen_events(n_events, tracks)
        for gif_num, gif in enumerate(gif_sets):
            for event_num in np.arange(1, n_events + 1):
                if gif['mode'] == 'ebe':
                    sd, skew, kurt = hist_events_ebe(events[:event_num], tracks, bin_width, gif['sample'])
                    if not math.isnan(sd):
                        deviations[gif_num]['sd'][event_num].append((sd - binom_sd)**2)
                    if not math.isnan(skew):
                        deviations[gif_num]['skew'][event_num].append((skew - binom_skew) ** 2)
                    if not math.isnan(kurt):
                        deviations[gif_num]['kurt'][event_num].append((kurt - binom_kurt)**2)
                elif gif['mode'] == 'parallel':
                    if event_num >= parallel_min_event:
                        his, sd, kurt, skew = hist_events_parallel(events[:event_num], tracks, bin_width, gif['sample'])
                        if not math.isnan(sd):
                            deviations[gif_num]['sd'][event_num].append((sd - binom_sd) ** 2)
                        if not math.isnan(sd):
                            deviations[gif_num]['skew'][event_num].append((skew - binom_skew) ** 2)
                        if not math.isnan(sd):
                            deviations[gif_num]['kurt'][event_num].append((kurt - binom_kurt) ** 2)
                else:
                    hist, sd, kurt, skew = hist_events(events[:event_num], tracks, bin_width, gif['sample'])
                    if not math.isnan(sd):
                        deviations[gif_num]['sd'][event_num].append((sd - binom_sd) ** 2)
                    if not math.isnan(sd):
                        deviations[gif_num]['skew'][event_num].append((skew - binom_skew) ** 2)
                    if not math.isnan(sd):
                        deviations[gif_num]['kurt'][event_num].append((kurt - binom_kurt) ** 2)

    axs[0].set_title(f'Deviation from Binomial Moments\n'
                     f'{bin_width}° Bin  {tracks} Tracks  {n_events} Events  {trials} Trials')
    axs[0].set_xlim(0, n_events)
    axs[0].set_ylabel('Standard Deviation')
    # axs[0].set_ylim(0, 0.5)

    axs[1].set_xlim(0, n_events)
    axs[1].set_ylabel('Skewness')
    # axs[0].set_ylim(0, 0.5)

    axs[2].set_xlim(0, n_events)
    axs[2].set_ylabel('Kurtosis')
    axs[2].set_xlabel('Number of Events')
    # axs[0].set_ylim(0, 0.5)

    for gif_num, gif in enumerate(gif_sets):
        for stat in deviations[gif_num]:
            for event_num in deviations[gif_num][stat]:
                # Should double check these calculations to be sure they're right, especially the error
                dev_errs[gif_num][stat][event_num] = np.sqrt(np.std(deviations[gif_num][stat][event_num]) /
                                                             np.sqrt(len(deviations[gif_num][stat][event_num])))
                deviations[gif_num][stat][event_num] = np.sqrt(np.mean(deviations[gif_num][stat][event_num]))

        sd_x = np.asarray(list(deviations[gif_num]['sd'].keys()))
        sd_y = np.asarray(list(deviations[gif_num]['sd'].values()))
        sd_err = np.asarray(list(dev_errs[gif_num]['sd'].values()))
        axs[0].plot(sd_x, sd_y, color=gif['color'], label=f'{gif["mode"]}, {gif["sample"]} samples')
        axs[0].fill_between(sd_x, sd_y - sd_err, sd_y + sd_err, facecolor=gif['color'], alpha=alpha, interpolate=True)

        skew_x = np.asarray(list(deviations[gif_num]['skew'].keys()))
        skew_y = np.asarray(list(deviations[gif_num]['skew'].values()))
        skew_err = np.asarray(list(dev_errs[gif_num]['skew'].values()))
        axs[1].plot(skew_x, skew_y, color=gif['color'], label=f'{gif["mode"]}, {gif["sample"]} samples')
        axs[1].fill_between(skew_x, skew_y - skew_err, skew_y + skew_err, facecolor=gif['color'], alpha=alpha)

        kurt_x = np.asarray(list(deviations[gif_num]['kurt'].keys()))
        kurt_y = np.asarray(list(deviations[gif_num]['kurt'].values()))
        kurt_err = np.asarray(list(dev_errs[gif_num]['kurt'].values()))
        axs[2].plot(kurt_x, kurt_y, color=gif['color'], label=f'{gif["mode"]}, {gif["sample"]} samples')
        axs[2].fill_between(kurt_x, kurt_y - kurt_err, kurt_y + kurt_err, facecolor=gif['color'], alpha=alpha)

    axs[0].legend()
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    plt.show()


def gif_events_maker():
    mode = 'save'
    gif_sets = [{'bin_width': 120, 'sd_lim': [0.8, 1.2], 'kurt_lim': [0.5, 1.5], 'max_sample': 150, 'dpi': 70,
                 'events': 10, 'tracks': 15},
                {'bin_width': 120, 'sd_lim': [0.8, 1.2], 'kurt_lim': [0.5, 1.5], 'max_sample': 150, 'dpi': 70,
                 'events': 1, 'tracks': 15}
                ]
    for gif in gif_sets:
        make_events_gif(gif, mode)


def gif_maker():
    mode = 'show'
    gif_sets = [{'bin_width': 30, 'sd_lim': [0.8, 1.2], 'kurt_lim': [0.5, 1.5], 'max_sample': 250, 'dpi': 70},
                {'bin_width': 120, 'sd_lim': [0.8, 1.2], 'kurt_lim': [0.5, 1.5], 'max_sample': 250, 'dpi': 70}]
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
        stats = DistStats(hist)
        kurts.append(stats.get_kurtosis().val)
        stds.append(stats.get_sd().val)
        norm = sum(hist.values())
        for key in hist:
            hist[key] /= norm
        hists.append(hist)
    hist_bins = np.arange(-0.5, len(angles) + 0.5)
    print('Data generated. Animating...')

    anim = FuncAnimation(fig, update, frames=np.arange(1, gif['max_sample']), interval=50, repeat_delay=10000,
                         fargs=(axs, hists, hist_bins, plot_samples, np.asarray(stds), np.asarray(kurts),
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
    hist = []
    for hist_bin in hist_bins:
        if hist_bin in hists[samples]:
            hist.append(hists[samples][hist_bin])
        else:
            hist.append(0)
    axs[0].bar(hist_bins, hist, width=1, align='center')

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
        norm = sum(events_hist[0].values())
        for key in events_hist[0]:
            events_hist[0][key] /= norm
        hists.append(events_hist[0])
    hist_bins = np.arange(0, gif['tracks'] + 1)
    print('Data generated. Animating...')

    anim = FuncAnimation(fig, update_binom, frames=range(gif['max_sample']), interval=50, repeat_delay=10000,
                         fargs=(
                             axs, hists, hist_bins, plot_samples, np.asarray(stds), np.asarray(kurts),
                             gif['bin_width'], gif['sd_lim'], gif['kurt_lim'], gif['max_sample'], gif['tracks'],
                        gif['events']))

    if mode == 'save':
        anim.save(f'/home/dylan/Research/Results/Presentations/11-10-20/'
                  f'{gif["bin_width"]}_degree_{gif["events"]}_events.gif', dpi=gif['dpi'], writer='imagemagick')
    else:
        plt.show()


def update_binom(i, axs, hists, hist_bins, plot_samples, stds, kurts, bin_width, sd_lim, kurt_lim, max_sample,
                 max_tracks, events):
    samples = i + 1
    print(f'{samples} Samples')
    axs[0].clear()
    axs[0].set_title(f'{bin_width}° Bin  {max_tracks} Tracks  {events} Events |'
                     f'  {samples} Samples  <-->  {360/samples:.2f}° Steps')
    axs[0].set_xlabel('Particles in Bin')
    axs[0].set_ylabel('Probability Density')
    axs[0].set_ylim(0, 1)
    hist = []
    for hist_bin in hist_bins:
        if hist_bin in hists[i]:
            hist.append(hists[i][hist_bin])
        else:
            hist.append(0)
    axs[0].bar(hist_bins, hist, width=1, align='center')

    n = max_tracks
    p = bin_width / 360
    q = 1 - p

    axs[1].clear()
    axs[1].set_xlim(0, max_sample)
    axs[1].set_ylabel('Standard Deviation')
    axs[1].plot(plot_samples[:i], stds[:i])
    binom_sd = np.sqrt(n * p * q)
    axs[1].set_ylim(binom_sd * sd_lim[0], binom_sd * sd_lim[1])
    axs[1].axhline(y=binom_sd, color='black', linestyle='-')

    axs[2].clear()
    axs[2].set_xlim(0, max_sample)
    axs[2].set_ylabel('Kurtosis')
    axs[2].set_xlabel('Number of Samples')
    axs[2].plot(plot_samples[:i], kurts[:i], color='red')
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    axs[2].set_ylim((binom_kurt + 3) * kurt_lim[0] - 3, (binom_kurt + 3) * kurt_lim[1] - 3)
    axs[2].axhline(y=binom_kurt, color='black', linestyle='-')

    return axs


def get_hist(angles, bin_width, samples):
    bin_step = 360 / samples

    bin_content = {}

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
        if count in bin_content:
            bin_content[count] += 1
        else:
            bin_content.update({count: 1})

    return bin_content


def get_hist_parallel(angles, bin_width, samples):
    bin_step = 360 / samples

    bin_content = []

    for start in np.arange(0, 360, bin_step):
        bin_content.append({})
        end = start + bin_width
        count = 0
        for angle in angles:
            if start <= angle < end:
                count += 1
        if end > 360:
            for angle in angles:
                if start <= angle + 360 < end:
                    count += 1
        if count in bin_content[-1]:
            bin_content[-1][count] += 1
        else:
            bin_content[-1].update({count: 1})

    return bin_content


def hist_events(events, num_tracks, bin_width, samples):
    hist = {}
    for event in events:
        for key, val in get_hist(event, bin_width, samples).items():
            if key in hist:
                hist[key] += val
            else:
                hist.update({key: val})
    stats = DistStats(hist)
    kurt = stats.get_kurtosis().val
    skew = stats.get_skewness().val
    sd = stats.get_sd().val
    # print(f'{bin_width}° bins,  {samples} samples,  {len(events)} events,  {num_tracks} tracks,  sd = {sd},  kurt = {kurt}')

    return hist, sd, kurt, skew


def hist_events_ebe(events, num_tracks, bin_width, samples):
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

    return sd, kurt, skew


def hist_events_parallel(events, num_tracks, bin_width, samples):
    hists = [{} for i in range(samples)]
    for event in events:
        for i, hist in enumerate(get_hist_parallel(event, bin_width, samples)):
            for key, val in hist.items():
                if key in hists[i]:
                    hists[i][key] += val
                else:
                    hists[i].update({key: val})
    kurts = []
    skews = []
    sds = []
    for hist in hists:
        stats = DistStats(hist)
        kurts.append(stats.get_kurtosis().val)
        skews.append(stats.get_skewness().val)
        sds.append(stats.get_sd().val)

    kurt = np.mean(kurts)
    skew = np.mean(skews)
    sd = np.mean(sds)

    return hist, sd, kurt, skew


def gen_events(num_events, num_tracks):
    events = []
    for event in range(num_events):
        event = gen_event(num_tracks)
        events.append(event)

    return events


def gen_events_stats(num_events, num_tracks, bin_width, samples):
    event_stats = {'mean': [], 'sd': [], 'skew': [], 'kurt': []}
    for event in range(num_events):
        event = gen_event(num_tracks)
        stats = DistStats(get_hist(event, bin_width, samples))
        event_stats['kurt'].append(stats.get_kurtosis().val)
        event_stats['skew'].append(stats.get_skewness().val)
        event_stats['sd'].append(stats.get_sd().val)
        event_stats['mean'].append(stats.get_mean().val)

    return event_stats


def gen_event(num_tracks):
    phis = []
    for track in range(num_tracks):
        phis.append(random.random() * 360)

    return phis


if __name__ == '__main__':
    main()
