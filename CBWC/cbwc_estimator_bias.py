#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 03 4:07 PM 2021
Created in PyCharm
Created as QGP_Scripts/cbwc_estimator_bias

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import binom, poisson, skellam
from Analyzer.DistStats import DistStats
import time
from multiprocessing import Pool
from functools import partial


def main():
    plt.rcParams['figure.autolayout'] = True
    # num_events = np.asarray(np.arange(2, 101, 1))
    # percentiles = [5, 30, 50, 70, 95]

    # n, p = 20, 0.4
    # q = 1 - p
    # dist = binom(n, q)
    # binning = np.arange(-0.5, 20 + 1.5, 1)
    # trials = 100
    # moment_pars = {'c2': {'method': lambda x: x.get_cumulant(2), 'true': n*p*q},
    #                'c3': {'method': lambda x: x.get_cumulant(3), 'true': n*p*q*(1-2*p)},
    #                'c4': {'method': lambda x: x.get_cumulant(4),
    #                       'true': n * p * q * (1 + (3 * n - 6) * p * q) - 3 * (n*p*q)**2},
    #                'k2': {'method': lambda x: x.get_k_stat(2), 'true': n*p*q},
    #                'k3': {'method': lambda x: x.get_k_stat(3), 'true': n*p*q*(1-2*p)},
    #                'k4': {'method': lambda x: x.get_k_stat(4),
    #                       'true': n * p * q * (1 + (3 * n - 6) * p * q) - 3 * (n*p*q)**2}}

    # num_events = np.asarray(np.arange(10, 101, 1))
    # percentiles = []
    #
    # mu = 5
    # dist = poisson(mu)
    # binning = np.arange(-0.5, mu * 10 + 1.5, 1)
    # trials = 1000
    # moment_pars = {'c2': {'method': lambda x: x.get_cumulant(2).val, 'true': mu},
    #                'c3': {'method': lambda x: x.get_cumulant(3).val, 'true': mu},
    #                'c4': {'method': lambda x: x.get_cumulant(4).val,
    #                       'true': mu},
    #                'k2': {'method': lambda x: x.get_k_stat(2).val, 'true': mu},
    #                'k3': {'method': lambda x: x.get_k_stat(3).val, 'true': mu},
    #                'k4': {'method': lambda x: x.get_k_stat(4).val,
    #                       'true': mu},
    #                'c4/c2': {'method': lambda x: x.get_cumulant(4).val / x.get_cumulant(2).val,
    #                          'true': 1},
    #                'k4/k2': {'method': lambda x: x.get_k_stat(4).val / x.get_k_stat(2).val,
    #                          'true': 1}
    #                }

    num_events = np.asarray(np.arange(10, 250, 5))
    threads = 13
    # num_events = np.asarray(np.arange(50, 5000, 10))
    percentiles = []

    mu1, mu2 = 20.9, 0.2
    mu_bar = (mu1 + mu2) / 2
    mu_delta = mu1 - mu2
    dist = skellam(mu1, mu2)
    binning = np.arange(-0.5, mu_bar * 2 * 10 + 1.5, 1)
    trials = 10000
    moment_pars = {'c2': {'method': get_c2, 'method_single': get_c2_meas, 'true': 2 * mu_bar},
                   'c3': {'method': get_c3, 'method_single': get_c3_meas, 'true': mu_delta},
                   'c4': {'method': get_c4, 'method_single': get_c4_meas, 'true': 2 * mu_bar},
                   'c5': {'method': get_c5, 'method_single': get_c5_meas, 'true': mu_delta},
                   'c6': {'method': get_c6, 'method_single': get_c6_meas, 'true': 2 * mu_bar},
                   'k2': {'method': get_k2, 'method_single': get_k2_meas, 'true': 2 * mu_bar},
                   'k3': {'method': get_k3, 'method_single': get_k3_meas, 'true': mu_delta},
                   'k4': {'method': get_k4, 'method_single': get_k4_meas, 'true': 2 * mu_bar},
                   'k5': {'method': get_k5, 'method_single': get_k5_meas, 'true': mu_delta},
                   'k6': {'method': get_k6, 'method_single': get_k6_meas, 'true': 2 * mu_bar},
                   'c4/c2': {'method': get_c4_div_c2, 'method_single': get_c4_div_c2_meas, 'true': 1},
                   'k4/k2': {'method': get_k4_div_k2, 'method_single': get_k4_div_k2_meas, 'true': 1},
                   'c6/c2': {'method': get_c6_div_c2, 'method_single': get_c6_div_c2_meas, 'true': 1},
                   'k6/k2': {'method': get_k6_div_k2, 'method_single': get_k6_div_k2_meas, 'true': 1},
                   'c4/c2 - k4/k2': {'method': get_c4_div_c2_sub_k4_div_k2,
                                     'method_single': get_c4_div_c2_sub_k4_div_k2_meas, 'true': 1},
                   'c6/c2 - k6/k2': {'method': get_c6_div_c2_sub_k6_div_k2,
                                     'method_single': get_c6_div_c2_sub_k6_div_k2_meas, 'true': 1},
                   }

    save_path = '/home/dylan/Desktop/'

    # demo_plots(dist)
    sim_single_trial(dist, num_events, binning, moment_pars, threads)
    # sim_trials(dist, num_events, trials, binning, percentiles, moment_pars, threads)
    # start = time.time()
    # event_means, event_errs, event_percs = simulate(dist, num_events, trials, binning, moment_pars, percentiles)
    # print(f'Simulation time: {time.time() - start}s')
    # # plot_moments(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    # plot_moments_together(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    # # plot_cumulants(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    # plot_cumulant_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    # plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)

    print('donzo')


def sim_trials(dist, num_events, trials, binning, percentiles, moment_pars, threads):
    start = time.time()
    event_means, event_errs, event_percs = simulate(dist, num_events, trials, binning, moment_pars, percentiles, threads)
    print(f'Simulation time: {time.time() - start}s')
    # plot_moments(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    plot_moments_together(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    # plot_cumulants(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    plot_cumulant_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)


def sim_single_trial(dist, num_events, binning, moment_pars, threads):
    start = time.time()
    event_means, event_errs = simulate_single(dist, num_events, binning, moment_pars, threads)
    print(f'Simulation time: {time.time() - start}s')
    plot_moments_together_single(num_events, event_means, event_errs, moment_pars)
    plot_ratios_single(num_events, event_means, event_errs, moment_pars)


def demo_plots(dist):
    plot_pmf(dist)
    plot_pmf_sample(dist, 100)
    plot_pmf_sample(dist, 1000)
    plot_c2_vs_trials(dist, 10000, 200)
    plot_k2_vs_trials(dist, 10000, 200)
    # plot_k2_est_vs_n()


def simulate(dist, num_events, trials, binning, moment_pars, percentiles, threads):
    event_means = {x: [] for x in moment_pars.keys()}
    event_errs = {x: [] for x in moment_pars.keys()}
    event_percs = {x: {y: [] for y in percentiles} for x in moment_pars.keys()}

    # 0 element events, 1 random state
    nevents_state = list(zip(num_events, [np.random.RandomState() for i in range(len(num_events))]))
    func = partial(sim_events, dist, trials, moment_pars, binning, percentiles)
    with Pool(threads) as p:
        trial_stats = p.map(func, nevents_state)

    # trial_stats = []
    # for n in num_events:
    #     trial_stats.append(sim_events(dist, trials, moment_pars, binning, percentiles, n))

    print("Combining events and plotting...")

    for trial_mean, trial_err, trial_perc in trial_stats:  # Iterate over num_events
        for moment, mean in trial_mean.items():
            event_means[moment].append(mean)
            event_errs[moment].append(trial_err[moment])
            for perc, perc_val in trial_perc[moment].items():
                event_percs[moment][perc].append(perc_val)

    for moment, means in event_means.items():
        event_means[moment] = np.asarray(means)
        event_errs[moment] = np.asarray(event_errs[moment])
        for perc, perc_val in event_percs[moment].items():
            event_percs[moment][perc] = np.asarray(perc_val)

    return event_means, event_errs, event_percs


def sim_events(dist, trials, moment_pars, binning, percentiles, nevents_state):
    events = nevents_state[0]
    print(f'Events: {events}')
    trial_moments = {x: [] for x in moment_pars.keys()}
    for trial in range(trials):
        hist, bin_edges = np.histogram(dist.rvs(size=events, random_state=nevents_state[1]), binning)
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        for moment, moment_val in calc_moments(hist, moment_pars).items():
            trial_moments[moment].append(moment_val)
    return comb_trials(trial_moments, percentiles)


def simulate_single(dist, num_events, binning, moment_pars, threads):
    event_means = {x: [] for x in moment_pars.keys()}
    event_errs = {x: [] for x in moment_pars.keys()}

    # 0 element events, 1 random state
    nevents_state = list(zip(num_events, [np.random.RandomState() for i in range(len(num_events))]))
    func = partial(sim_events_single, dist, moment_pars, binning)
    with Pool(threads) as p:
        stats = p.map(func, nevents_state)

    # stats = []
    # for n in num_events:
    #     stats.append(sim_events_single(dist, moment_pars, binning, n))

    print("Combining events and plotting...")

    for means, errs, events in stats:  # Iterate over num_events, getting means and errs from each
        print(events)
        for moment, mean in means.items():  # Iterate over different moments for num_events
            event_means[moment].append(mean)
            event_errs[moment].append(errs[moment])

    for moment, means in event_means.items():
        event_means[moment] = np.asarray(means)
        event_errs[moment] = np.asarray(event_errs[moment])

    return event_means, event_errs


def sim_events_single(dist, moment_pars, binning, nevents_state):
    events = nevents_state[0]
    print(f'Events: {events}')
    means, errs = {}, {}
    hist, bin_edges = np.histogram(dist.rvs(size=events, random_state=nevents_state[1]), binning)
    hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
    for moment, moment_meas in calc_moments_single(hist, moment_pars).items():
        means.update({moment: moment_meas.val})
        errs.update({moment: moment_meas.err})

    return means, errs, events


def calc_moments(hist, moment_pars):
    moments = {}
    stats = DistStats(hist)
    for moment, pars in moment_pars.items():
        moments.update({moment: pars['method'](stats)})

    return moments


def calc_moments_single(hist, moment_pars):
    moments = {}
    stats = DistStats(hist)
    for moment, pars in moment_pars.items():
        if 'method_single' in pars:
            moments.update({moment: pars['method_single'](stats)})

    return moments


def comb_trials(trial_moments, percentiles):
    means, errs, percs = {}, {}, {}
    for moment, trial_vals in trial_moments.items():
        means.update({moment: np.mean(trial_vals)})
        errs.update({moment: np.std(trial_vals) / np.sqrt(len(trial_vals))})
        perc = np.percentile(trial_vals, percentiles)
        percs.update({moment: dict(zip(percentiles, perc))})

    return means, errs, percs


def plot_moments(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    for i in range(2, 7):
        c, k = 'c' + str(i), 'k' + str(i)
        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title('Order ' + str(i))

        ax.axhline(moment_pars[k]['true'], color='r', ls='--', label='Analytical')

        for j in range(int(len(percs) / 2)):
            ax.fill_between(num_events, event_percs[c][percs[j]], event_percs[c][percs[len(percs) - 1 - j]],
                            color='b', alpha=0.3)
            ax.fill_between(num_events, event_percs[k][percs[j]], event_percs[k][percs[len(percs) - 1 - j]],
                            color='g', alpha=0.3)

        ax.fill_between(num_events, event_means[c] + event_errs[c], event_means[c] - event_errs[c],
                        color='b', alpha=0.3)
        ax.fill_between(num_events, event_means[k] + event_errs[k], event_means[k] - event_errs[k],
                        color='g', alpha=0.3)

        ax.plot(num_events, event_means[c], color='b', label=f'{c}')
        ax.plot(num_events, event_means[k], color='g', label=f'{k}')
        ax.legend()
        ax.set_xlabel('Sample Size n')
    plt.show()


def plot_moments_together(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    fig = plt.figure()
    gs1 = gridspec.GridSpec(5, 1)
    gs1.update(hspace=0.0)
    # axs = gs1.subplots(sharex=True)
    # # fig, axs = plt.subplots(5, 1, sharex=True)
    # plt.subplots_adjust(hspace=0.0, wspace=0.0)
    fig.canvas.manager.set_window_title('Cumulant KStat Compare Threads')
    for i in range(2, 7):
        ax = plt.subplot(gs1[i - 2])
        c, k = 'c' + str(i), 'k' + str(i)

        if i == 2:
            ax.axhline(moment_pars[k]['true'], color='r', ls='--', label='Analytical')
        else:
            ax.axhline(moment_pars[k]['true'], color='r', ls='--')

        for j in range(int(len(percs) / 2)):
            ax.fill_between(num_events, event_percs[c][percs[j]], event_percs[c][percs[len(percs) - 1 - j]],
                            color='b', alpha=0.3)
            ax.fill_between(num_events, event_percs[k][percs[j]], event_percs[k][percs[len(percs) - 1 - j]],
                            color='g', alpha=0.3)

        ax.fill_between(num_events, event_means[c] + event_errs[c], event_means[c] - event_errs[c],
                        color='b', alpha=0.3)
        ax.fill_between(num_events, event_means[k] + event_errs[k], event_means[k] - event_errs[k],
                        color='g', alpha=0.3)

        ax.plot(num_events, event_means[c], color='b', label=f'{c}')
        ax.plot(num_events, event_means[k], color='g', label=f'{k}')
        ax.legend()
        if i == 6:
            ax.set_xlabel('Sample Size n')
    plt.show()


def plot_moments_together_single(num_events, event_means, event_errs, moment_pars):
    fig = plt.figure()
    gs1 = gridspec.GridSpec(5, 1)
    gs1.update(hspace=0.0)
    # axs = gs1.subplots(sharex=True)
    # # fig, axs = plt.subplots(5, 1, sharex=True)
    # plt.subplots_adjust(hspace=0.0, wspace=0.0)
    fig.canvas.manager.set_window_title('Cumulant KStat Compare')
    for i in range(2, 7):
        ax = plt.subplot(gs1[i - 2])
        c, k = 'c' + str(i), 'k' + str(i)

        if i == 2:
            ax.axhline(moment_pars[k]['true'], color='r', ls='--', label='Analytical')
        else:
            ax.axhline(moment_pars[k]['true'], color='r', ls='--')

        ax.errorbar(num_events, event_means[c], yerr=event_errs[c], color='b', marker='_', ls='none', label=f'{c}')
        ax.errorbar(num_events, event_means[k], yerr=event_errs[k], color='g', marker='_', ls='none', label=f'{k}')
        ax.legend()
        if i == 6:
            ax.set_xlabel('Sample Size n')
    plt.show()


def plot_cumulants(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    fig, axs = plt.subplots(5, 1, sharex=True, hspace=0.0)
    fig.canvas.manager.set_window_title('Cumulants')
    for i in range(2, 7):
        c = 'c' + str(i)

        axs[i - 2].axhline(moment_pars[c]['true'], color='r', ls='--', label='Analytical')

        for j in range(int(len(percs) / 2)):
            axs[i - 2].fill_between(num_events, event_percs[c][percs[j]], event_percs[c][percs[len(percs) - 1 - j]],
                                    color='b', alpha=0.3)

        axs[i - 2].fill_between(num_events, event_means[c] + event_errs[c], event_means[c] - event_errs[c],
                                color='b', alpha=0.3)

        axs[i - 2].plot(num_events, event_means[c], color='b', label='Simulation')
        axs[i - 2].set_ylabel(c)

    axs[0].legend()
    axs[2].set_xlabel('Sample Size n')
    plt.show()


def plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    fig1, axs1 = plt.subplots(2, 1, sharex=True)
    fig1.canvas.manager.set_window_title('C4 / C2 with K4 / K2 Threads')

    axs1[0].axhline(moment_pars['k4/k2']['true'], color='r', ls='--', label='Analytical')
    axs1[0].fill_between(num_events, event_means['c4/c2'] - event_errs['c4/c2'], event_means['c4/c2'] + event_errs['c4/c2'],
                    color='b', alpha=0.3)
    axs1[0].fill_between(num_events, event_means['k4/k2'] - event_errs['k4/k2'], event_means['k4/k2'] + event_errs['k4/k2'],
                    color='g', alpha=0.3)

    axs1[0].plot(num_events, event_means['c4/c2'], color='b', label=f'C4/C2')
    axs1[0].plot(num_events, event_means['k4/k2'], color='g', label=f'K4/K2')
    axs1[0].legend()

    axs1[1].axhline(moment_pars['k6/k2']['true'], color='r', ls='--', label='Analytical')
    axs1[1].fill_between(num_events, event_means['c6/c2'] - event_errs['c6/c2'],
                         event_means['c6/c2'] + event_errs['c6/c2'],
                         color='b', alpha=0.3)
    axs1[1].fill_between(num_events, event_means['k6/k2'] - event_errs['k6/k2'],
                         event_means['k6/k2'] + event_errs['k6/k2'],
                         color='g', alpha=0.3)

    axs1[1].plot(num_events, event_means['c6/c2'], color='b', label=f'C6/C2')
    axs1[1].plot(num_events, event_means['k6/k2'], color='g', label=f'K6/K2')
    axs1[1].legend()
    axs1[1].set_xlabel('Sample Size n')

    fig2, axs2 = plt.subplots(2, 1, sharex=True)
    fig2.canvas.manager.set_window_title('(C4 / C2) - (K4 / K2) Threads')
    axs2[0].axhline(0, color='black', ls='--')
    axs2[0].fill_between(num_events, event_means['c4/c2 - k4/k2'] - event_errs['c4/c2 - k4/k2'],
                     event_means['c4/c2 - k4/k2'] + event_errs['c4/c2 - k4/k2'],
                     color='red', alpha=0.3)
    axs2[0].plot(num_events, event_means['c4/c2 - k4/k2'], color='r', label='C4/C2 - K4/K2')
    axs2[0].legend()
    axs2[0].grid()

    axs2[1].axhline(0, color='black', ls='--')
    axs2[1].fill_between(num_events, event_means['c6/c2 - k6/k2'] - event_errs['c6/c2 - k6/k2'],
                         event_means['c6/c2 - k6/k2'] + event_errs['c6/c2 - k6/k2'],
                         color='red', alpha=0.3)
    axs2[1].plot(num_events, event_means['c6/c2 - k6/k2'], color='r', label='C6/C2 - K6/K2')
    axs2[1].legend()
    axs2[1].grid()
    axs2[1].set_xlabel('Sample Size n')

    plt.show()


def plot_ratios_single(num_events, event_means, event_errs, moment_pars):
    fig1, axs1 = plt.subplots(2, 1, sharex=True)
    fig1.canvas.manager.set_window_title('C4 / C2 with K4 / K2')

    axs1[0].axhline(moment_pars['k4/k2']['true'], color='r', ls='--', label='Analytical')

    axs1[0].errorbar(num_events, event_means['c4/c2'], yerr=event_errs['c4/c2'], marker='o', ls='none', color='b',
                     label=f'C4/C2')
    axs1[0].errorbar(num_events, event_means['k4/k2'], yerr=event_errs['k4/k2'], marker='o', ls='none', color='g',
                     label=f'K4/K2')
    axs1[0].legend()

    axs1[1].axhline(moment_pars['k6/k2']['true'], color='r', ls='--', label='Analytical')

    axs1[1].errorbar(num_events, event_means['c6/c2'], yerr=event_errs['c6/c2'], marker='o', ls='none', color='b',
                     label=f'C6/C2')
    axs1[1].errorbar(num_events, event_means['k6/k2'], yerr=event_errs['k6/k2'], marker='o', ls='none', color='g',
                     label=f'K6/K2')
    axs1[1].legend()
    axs1[1].set_xlabel('Sample Size n')

    fig2, axs2 = plt.subplots(2, 1, sharex=True)
    fig2.canvas.manager.set_window_title('(C4 / C2) - (K4 / K2)')
    axs2[0].axhline(0, color='black', ls='--')
    axs2[0].errorbar(num_events, event_means['c4/c2 - k4/k2'], yerr=event_errs['c4/c2 - k4/k2'], marker='o', ls='none',
                     color='r', label='C4/C2 - K4/K2')
    axs2[0].legend()
    axs2[0].grid()

    axs2[1].axhline(0, color='black', ls='--')
    axs2[1].errorbar(num_events, event_means['c6/c2 - k6/k2'], yerr=event_errs['c6/c2 - k6/k2'], marker='o', ls='none',
                     color='r', label='C6/C2 - K6/K2')
    axs2[1].legend()
    axs2[1].grid()
    axs2[1].set_xlabel('Sample Size n')

    plt.show()


def plot_cumulant_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title('C4 / C2')

    ax.axhline(moment_pars['c4/c2']['true'], color='r', ls='--', label='Analytical')
    ax.fill_between(num_events, event_means['c4/c2'] - event_errs['c4/c2'],
                    event_means['c4/c2'] + event_errs['c4/c2'],
                    color='b', alpha=0.3)

    ax.plot(num_events, event_means['c4/c2'], color='b', label=f'C4/C2')
    ax.legend()
    ax.set_xlabel('Sample Size n')


def plot_pmf(dist):
    fig, ax = plt.subplots()
    x_plot = np.arange(dist.ppf(0.001), dist.ppf(0.999))
    ax.plot(x_plot, dist.pmf(x_plot), color='red', lw=4,
            label=fr'Skellam PDF $\mu_1$={dist.args[0]}, $\mu_2$={dist.args[1]}')
    ax.set_xlabel('Distribution Value')
    ax.set_ylabel('Probability Mass')
    ax.legend()
    ax.grid()
    plt.show()


def plot_pmf_sample(dist, num):
    fig, ax = plt.subplots()
    x_plot = np.arange(dist.ppf(0.001), dist.ppf(0.999))
    ax.plot(x_plot, dist.pmf(x_plot) * num, color='red', lw=4,
            label=fr'Skellam PDF $\mu_1$={dist.args[0]}, $\mu_2$={dist.args[1]}')
    x_bin = list(x_plot - 0.5)
    x_bin.append(x_bin[-1] + 1)
    ax.hist(dist.rvs(size=num), np.asarray(x_bin), color='blue', label=f'Sample with n={num}')
    ax.set_xlabel('Distribution Value')
    ax.set_ylabel('Counts')
    ax.legend()
    plt.show()


def plot_c2_vs_trials(dist, num, trials):
    x_plot = np.arange(dist.ppf(0.000000001), dist.ppf(0.999999999))
    x_bin = list(x_plot - 0.5)
    x_bin.append(x_bin[-1] + 1)
    c2s = []
    for trial in range(trials):
        hist, bin_edges = np.histogram(dist.rvs(size=num), x_bin)
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        stats = DistStats(hist)
        c2s.append(stats.get_cumulant(2).val)
    fig, ax = plt.subplots()
    ax.axhline(sum(dist.args), color='red', ls='--', label=f'True Value: {sum(dist.args):.1f}')
    mean = np.mean(c2s)
    sem = np.std(c2s) / np.sqrt(len(c2s))
    ax.axhspan(mean + sem, mean - sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
    ax.axhline(mean, color='green', ls='--')
    ax.scatter(range(trials), c2s)
    ax.set_xlabel('Trial Number (arb)')
    ax.set_ylabel('C2')
    ax.legend()

    fig2, ax2 = plt.subplots()
    ax2.hist(c2s)
    ax2.axvline(sum(dist.args), color='red', ls='--', label=f'True Value: {sum(dist.args):.1f}')
    ax2.axvspan(mean + sem, mean - sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
    ax2.axvline(mean, color='green', ls='--')
    ax2.set_xlabel('C2')
    ax2.legend()
    plt.show()


def plot_k2_vs_trials(dist, num, trials):
    x_plot = np.arange(dist.ppf(0.000000001), dist.ppf(0.999999999))
    x_bin = list(x_plot - 0.5)
    x_bin.append(x_bin[-1] + 1)
    k2s = []
    for trial in range(trials):
        hist, bin_edges = np.histogram(dist.rvs(size=num), x_bin)
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        stats = DistStats(hist)
        k2s.append(stats.get_k_stat(2).val)
    fig, ax = plt.subplots()
    ax.axhline(sum(dist.args), color='red', ls='--', label=f'True Value: {sum(dist.args):.1f}')
    mean = np.mean(k2s)
    sem = np.std(k2s) / np.sqrt(len(k2s))
    ax.axhspan(mean + sem, mean - sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
    ax.axhline(mean, color='green', ls='--')
    ax.scatter(range(trials), k2s)
    ax.set_xlabel('Trial Number (arb)')
    ax.set_ylabel('K2')
    ax.legend()

    fig2, ax2 = plt.subplots()
    ax2.hist(k2s)
    ax2.axvline(sum(dist.args), color='red', ls='--', label=f'True Value: {sum(dist.args):.1f}')
    ax2.axvspan(mean + sem, mean - sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
    ax2.axvline(mean, color='green', ls='--')
    ax2.set_xlabel('K2')
    ax2.legend()

    plt.show()


def plot_k2_est_vs_n():
    mu2 = 11
    mu4 = 372
    n = np.arange(500, 5000, 10)
    m4 = (n - 1) * (3 * (2 * n - 3) * mu2 ** 2 + (n ** 2 - 3 * n + 3) * mu4) / n ** 3
    m2 = (n - 1) * mu2 / n
    k4s = n ** 2 * ((n + 1) * m4 - 3 * (n - 1) * m2 ** 2) / ((n - 1) * (n - 2) * (n - 3))
    k2s = n / (n - 1) * m2

    fig, ax = plt.subplots()
    ax.axhline(1.0, color='red', ls='--')
    ax.plot(n, (k4s / k2s) / 0.8)
    plt.show()


def get_c2(x):
    return x.get_cumulant(2).val


def get_c3(x):
    return x.get_cumulant(3).val


def get_c4(x):
    return x.get_cumulant(4).val


def get_c5(x):
    return x.get_cumulant(5).val


def get_c6(x):
    return x.get_cumulant(6).val


def get_k2(x):
    return x.get_k_stat(2).val


def get_k3(x):
    return x.get_k_stat(3).val


def get_k4(x):
    return x.get_k_stat(4).val


def get_k5(x):
    return x.get_k_stat(5).val


def get_k6(x):
    return x.get_k_stat(6).val


def get_c4_div_c2(x):
    return x.get_cumulant(4).val / x.get_cumulant(2).val


def get_k4_div_k2(x):
    return x.get_k_stat(4).val / x.get_k_stat(2).val


def get_c6_div_c2(x):
    return x.get_cumulant(6).val / x.get_cumulant(2).val


def get_k6_div_k2(x):
    return x.get_k_stat(6).val / x.get_k_stat(2).val


def get_c4_div_c2_sub_k4_div_k2(x):
    return x.get_cumulant(4).val / x.get_cumulant(2).val - x.get_k_stat(4).val / x.get_k_stat(2).val


def get_c6_div_c2_sub_k6_div_k2(x):
    return x.get_cumulant(6).val / x.get_cumulant(2).val - x.get_k_stat(6).val / x.get_k_stat(2).val


def get_c2_meas(x):
    return x.get_cumulant(2)


def get_c3_meas(x):
    return x.get_cumulant(3)


def get_c4_meas(x):
    return x.get_cumulant(4)


def get_c5_meas(x):
    return x.get_cumulant(5)


def get_c6_meas(x):
    return x.get_cumulant(6)


def get_k2_meas(x):
    return x.get_k_stat(2)


def get_k3_meas(x):
    return x.get_k_stat(3)


def get_k4_meas(x):
    return x.get_k_stat(4)


def get_k5_meas(x):
    return x.get_k_stat(5)


def get_k6_meas(x):
    return x.get_k_stat(6)


def get_c4_div_c2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_cumulant(4) / x.get_cumulant(2)


def get_k4_div_k2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_k_stat(4) / x.get_k_stat(2)


def get_c6_div_c2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_cumulant(6) / x.get_cumulant(2)


def get_k6_div_k2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_k_stat(6) / x.get_k_stat(2)


def get_c4_div_c2_sub_k4_div_k2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_cumulant(4) / x.get_cumulant(2) - x.get_k_stat(4) / x.get_k_stat(2)


def get_c6_div_c2_sub_k6_div_k2_meas(x):  # Delta theorem errors will be wrong, need to apply directly to quantity
    return x.get_cumulant(6) / x.get_cumulant(2) - x.get_k_stat(6) / x.get_k_stat(2)


if __name__ == '__main__':
    main()
