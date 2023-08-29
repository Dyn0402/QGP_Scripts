#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 05 9:04 AM 2022
Created in PyCharm
Created as QGP_Scripts/pdf_examples.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from ClustDist import ClustDist, ClustDistPlusMinus
from sub_event_resample_algorithm import get_resamples4, plot_event_nobin


def main():
    # run_dynamic()
    # run_dynamic_plus_minus()
    cluster_voids_example_plots()
    plot_example_event_tracks()
    print('donzo')


def plot_example_event_tracks():
    n_tracks = 25
    bin_width = np.radians(120)
    resamples = 1  # Not used
    seed = 94

    # Clustering
    amp = 0.5
    spread = 1.0
    wraps = 8
    clust_pars = (amp, spread, wraps)
    tracks = sim_event_dynamic(clust_pars, n_tracks, bin_width, resamples, seed, return_tracks=True)[-1]
    print('Clustering')
    plot_event_nobin(tracks)

    # Anti-Clustering
    amp = -0.5
    aclust_pars = (amp, spread, wraps)
    tracks = sim_event_dynamic(aclust_pars, n_tracks, bin_width, resamples, seed, return_tracks=True)[-1]
    print('Anti-Clustering')
    plot_event_nobin(tracks)


def run_static():
    phis = np.linspace(0, 2 * np.pi, 25)[:-1]
    # phis = [1, 5]
    sd = 0.01
    cl_amp = 90
    bin_width = np.radians(120)
    tracks = 25
    events = 100
    samples = 72
    threads = 12
    seed = 3

    info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions, {tracks} Tracks/Event, {events} Events, {samples} Samples/Event'
    pdf = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=5)
    plot_pdf(pdf, f'Track PDF\n{info_string}')
    plt.show()

    hist = sim_pdf_static(pdf, events, tracks, bin_width, samples, threads, seed)
    plot_pdf(pdf, f'PDF {info_string}')
    plot_sim(hist, bin_width, f'Azimuthal Partition Multiplicity Distribution\n{info_string}')
    plt.show()


def run_dynamic():
    phis = np.linspace(0, 2 * np.pi, 25)[:-1]
    # phis = [1, 5]
    sd = 1.0
    cl_amp = -0.5
    bin_width = np.radians(120)
    tracks = 25
    events = 100
    samples = 72
    threads = 11
    seed = 61

    sim_pars = (cl_amp, sd, 5)

    # sim_event_dynamic_plot(sim_pars, tracks, bin_width, samples, seed)
    # return

    info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions, {tracks} Tracks/Event, {events} Events, {samples} Samples/Event'
    # pdf = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=5)
    # plot_pdf(pdf, f'Track PDF\n{info_string}')
    # plt.show()

    hist = sim_pdf_dynamic(sim_pars, events, tracks, bin_width, samples, threads, seed)
    # plot_pdf(pdf, f'PDF {info_string}')
    plot_sim(hist, bin_width, f'Azimuthal Partition Multiplicity Distribution\n{info_string}')
    plt.show()


def cluster_voids_example_plots():
    bin_width = np.radians(120)
    tracks = 25
    events = 100
    samples = 72
    threads = 11
    seed = 61

    plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.4

    fig, [ax_minus, ax_plus] = plt.subplots(1, 2, figsize=(13.33, 5), dpi=144)

    sd = 1.0
    cl_amp = -0.5
    sim_pars = (cl_amp, sd, 5)
    # info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions, {tracks} Tracks/Event, {events} Events, {samples} Samples/Event'
    info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions'
    hist = sim_pdf_dynamic(sim_pars, events, tracks, bin_width, samples, threads, seed)
    # title = f'Azimuthal Partition Multiplicity Distribution: Negative Correlation\n{info_string}'
    title = f'Negative Correlation'
    plot_sim(hist, bin_width, title, ax_sim=ax_minus)

    sd = 1.0
    cl_amp = 0.5
    sim_pars = (cl_amp, sd, 5)
    # info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions, {tracks} Tracks/Event, {events} Events, {samples} Samples/Event'
    info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions'
    hist = sim_pdf_dynamic(sim_pars, events, tracks, bin_width, samples, threads, seed)
    # title = f'Azimuthal Partition Multiplicity Distribution: Positive Correlation\n{info_string}'
    title = f'Positive Correlation'
    plot_sim(hist, bin_width, title, ax_sim=ax_plus)

    prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", shrinkA=0, shrinkB=0, color='black', linewidth=2)

    ax_minus.annotate("", xy=(3.9, 228), xytext=(2.15, 595), arrowprops=prop)
    ax_minus.text(2.15, 595 + 40, 'No voids', ha='center', size='large')
    ax_minus.annotate("", xy=(13.87, 182), xytext=(15.75, 595), arrowprops=prop)
    ax_minus.text(15.75, 595 + 40, 'No clusters', ha='center', size='large')
    # ax_minus.text(20.5, 1850, 'A = -0.5\nσ =  1.0', size='large', va='top')
    # ax_minus.text(19.75, 1590, 'A = -0.5\nσ =  1.0', size='large', va='top')

    ax_plus.annotate("", xy=(2.5, 665), xytext=(2.5, 950), arrowprops=prop)
    ax_plus.text(2.5, 950 + 20, 'Excess voids', ha='center', size='large')
    ax_plus.annotate("", xy=(19.5, 230), xytext=(22, 410), arrowprops=prop)
    ax_plus.text(22, 410 + 20, 'Excess clusters', ha='center', size='large')
    # ax_plus.text(20.5, 1000, 'A = +0.5\nσ =  1.0', size='large', va='top')
    # ax_plus.text(19.75, 875, 'A = +0.5\nσ =  1.0', size='large', va='top')

    fig.tight_layout()
    fig.subplots_adjust(top=0.94, bottom=0.105, left=0.065, right=0.995, wspace=0.175)

    plt.show()


def run_dynamic_plus_minus():
    # phis = np.linspace(0, 2 * np.pi, 2)[:-1]
    phis = [3]
    sd_minus, sd_plus = 2.0, 0.2
    amp_minus, amp_plus = -0.2, 0.2
    bin_width = np.radians(120)
    tracks = 2
    events = 100
    samples = 72
    threads = 15
    seed = 61

    sim_pars = (amp_minus, amp_plus, sd_minus, sd_plus, 5)

    # sim_event_dynamic_plot(sim_pars, tracks, bin_width, samples, seed)
    # return

    info_string = f'{int(np.rad2deg(bin_width) + 0.5)}° Partitions, {tracks} Tracks/Event, {events} Events, ' \
                  f'{samples} Samples/Event'
    pdf = ClustDistPlusMinus(phis, sd_minus, sd_plus, amp_minus, amp_plus, a=0, b=2 * np.pi, wrap_num=5)
    plot_pdf(pdf, f'Track PDF\n{info_string}')
    plt.show()

    hist = sim_pdf_dynamic(sim_pars, events, tracks, bin_width, samples, threads, seed)
    # plot_pdf(pdf, f'PDF {info_string}')
    plot_sim(hist, bin_width, f'Azimuthal Partition Multiplicity Distribution\n{info_string}')
    plt.show()


def plot_pdf(pdf, title='Pdf'):
    fig_pdf, ax_pdf = plt.subplots()
    ax_pdf.axhline(0, ls='--', color='black')
    ax_pdf.plot(pdf.x, pdf.pdf(pdf.x))
    ax_pdf.set_xlabel('Phi')
    ax_pdf.set_title(title)
    fig_pdf.tight_layout()
    fig_pdf.canvas.manager.set_window_title(title.replace('\n', '_'))


def sim_pdf_static(pdf, n_events, n_tracks, bin_width, samples, threads=1, seed=42):
    hist = np.zeros(n_tracks + 1, dtype=int)
    seeds = iter(np.random.SeedSequence(seed).spawn(n_events))
    jobs = [(pdf, n_tracks, bin_width, samples, next(seeds)) for event_i in range(n_events)]

    with Pool(threads) as pool:
        for hist_i in tqdm.tqdm(pool.istarmap(sim_event_static, jobs), total=len(jobs)):
            hist += hist_i

    return hist


def sim_event_static(pdf, n_tracks, bin_width, samples, seed):
    rng = np.random.default_rng(seed)
    tracks = np.sort(pdf.rvs(size=n_tracks, random_state=rng))
    hist = get_resamples4(tracks, bin_width, samples)
    # hist = np.histogram(hist, bins=np.arange(-0.5, n_tracks + 1.5, 1))[0]

    return hist


def sim_pdf_dynamic(sim_pars, n_events, n_tracks, bin_width, samples, threads=1, seed=42):
    hist = np.zeros(n_tracks + 1, dtype=int)
    seeds = iter(np.random.SeedSequence(seed).spawn(n_events))
    jobs = [(sim_pars, n_tracks, bin_width, samples, next(seeds)) for event_i in range(n_events)]

    with Pool(threads) as pool:
        for hist_i in tqdm.tqdm(pool.istarmap(sim_event_dynamic, jobs), total=len(jobs)):
            hist += hist_i

    return hist


def sim_event_dynamic(sim_pars, n_tracks, bin_width, samples, seed, return_tracks=False):
    rng = np.random.default_rng(seed)
    cl_amp, sd, wrap_num = sim_pars
    phis = []
    prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=wrap_num)
    while len(phis) < n_tracks:
        prob_dist.random_state = rng
        phis.append(prob_dist.rvs())
        prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=wrap_num)
    tracks = np.sort(phis)
    hist = get_resamples4(tracks, bin_width, samples, rng=rng)
    # print(hist)
    # hist = np.histogram(hist, bins=np.arange(-0.5, n_tracks + 1.5, 1))[0]
    # print(hist)

    if return_tracks:
        return hist, tracks

    return hist


def sim_pdf_dynamic_plus_minus(sim_pars, n_events, n_tracks, bin_width, samples, threads=1, seed=42):
    hist = np.zeros(n_tracks + 1, dtype=int)
    seeds = iter(np.random.SeedSequence(seed).spawn(n_events))
    jobs = [(sim_pars, n_tracks, bin_width, samples, next(seeds)) for event_i in range(n_events)]

    with Pool(threads) as pool:
        for hist_i in tqdm.tqdm(pool.istarmap(sim_event_dynamic_plus_minus, jobs), total=len(jobs)):
            hist += hist_i

    return hist


def sim_event_dynamic_plus_minus(sim_pars, n_tracks, bin_width, samples, seed):
    rng = np.random.default_rng(seed)
    amp_minus, amp_plus, sd_minus, sd_plus, wrap_num = sim_pars
    phis = []
    prob_dist = ClustDistPlusMinus(phis, sd_minus, sd_plus, amp_minus, amp_plus, a=0, b=2 * np.pi, wrap_num=wrap_num)
    while len(phis) < n_tracks:
        prob_dist.random_state = rng
        phis.append(prob_dist.rvs())
        prob_dist = ClustDistPlusMinus(phis, sd_minus, sd_plus, amp_minus, amp_plus, a=0, b=2 * np.pi,
                                       wrap_num=wrap_num)
    tracks = np.sort(phis)
    hist = get_resamples4(tracks, bin_width, samples)
    # print(hist)
    # hist = np.histogram(hist, bins=np.arange(-0.5, n_tracks + 1.5, 1))[0]
    # print(hist)

    return hist


def plot_sim(hist, bin_width, title='Protons in Bin Distribution', ax_sim=None):
    if ax_sim is None:
        fig_sim, ax_sim = plt.subplots(figsize=(6.67, 5), dpi=144)
        ax_passed = False
    else:
        ax_passed = True
    bin_centers = np.arange(0, len(hist), 1)
    ax_sim.bar(bin_centers, hist, width=1, align='center', label='Simulation')
    ax_sim.plot(bin_centers, binom.pmf(bin_centers, len(hist), bin_width / (2 * np.pi)) * np.sum(hist), color='red',
                label='Binomial')
    ax_sim.set_xlabel('Protons in Partition')
    ax_sim.set_ylabel('Number of Partitions')
    ax_sim.set_title(title)
    ax_sim.legend()
    if not ax_passed:
        fig_sim.tight_layout()
        fig_sim.canvas.manager.set_window_title(title.replace('\n', '_'))

    return ax_sim


if __name__ == '__main__':
    main()
