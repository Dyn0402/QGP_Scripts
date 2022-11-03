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
from sub_event_resample_algorithm import get_resamples4


def main():
    # run_dynamic()
    run_dynamic_plus_minus()
    print('donzo')


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
    threads = 15
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


def run_dynamic_plus_minus():
    # phis = np.linspace(0, 2 * np.pi, 2)[:-1]
    phis = [3]
    sd_minus, sd_plus = 1.0, 0.2
    amp_minus, amp_plus = -0.5, 0.2
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


def sim_event_dynamic(sim_pars, n_tracks, bin_width, samples, seed):
    rng = np.random.default_rng(seed)
    cl_amp, sd, wrap_num = sim_pars
    phis = []
    prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=wrap_num)
    while len(phis) < n_tracks:
        prob_dist.random_state = rng
        phis.append(prob_dist.rvs())
        prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=wrap_num)
    tracks = np.sort(phis)
    hist = get_resamples4(tracks, bin_width, samples)
    # print(hist)
    # hist = np.histogram(hist, bins=np.arange(-0.5, n_tracks + 1.5, 1))[0]
    # print(hist)

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


def plot_sim(hist, bin_width, title='Protons in Bin Distribution'):
    fig_sim, ax_sim = plt.subplots(figsize=(6.67, 5), dpi=144)
    bin_centers = np.arange(0, len(hist), 1)
    ax_sim.bar(bin_centers, hist, width=1, align='center', label='Simulation')
    ax_sim.plot(bin_centers, binom.pmf(bin_centers, len(hist), bin_width / (2 * np.pi)) * np.sum(hist), color='red',
                label='Binomial')
    ax_sim.set_xlabel('Protons in Partition')
    ax_sim.set_ylabel('Number of Partitions')
    ax_sim.set_title(title)
    ax_sim.legend()
    fig_sim.tight_layout()
    fig_sim.canvas.manager.set_window_title(title.replace('\n', '_'))


if __name__ == '__main__':
    main()
