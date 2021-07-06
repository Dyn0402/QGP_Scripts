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
import matplotlib
from scipy.stats import binom, poisson, skellam
from Analyzer.DistStats import DistStats
import time


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

    num_events = np.asarray(np.arange(10, 101, 1))
    percentiles = []

    mu1, mu2 = 20.9, 0.2
    mu_bar = (mu1 + mu2) / 2
    mu_delta = mu1 - mu2
    dist = skellam(mu1, mu2)
    binning = np.arange(-0.5, mu_bar * 2 * 10 + 1.5, 1)
    trials = 1000
    moment_pars = {'c2': {'method': lambda x: x.get_cumulant(2).val, 'true': 2 * mu_bar},
                   'c3': {'method': lambda x: x.get_cumulant(3).val, 'true': mu_delta},
                   'c4': {'method': lambda x: x.get_cumulant(4).val,
                          'true': 2 * mu_bar},
                   'k2': {'method': lambda x: x.get_k_stat(2).val, 'true': 2 * mu_bar},
                   'k3': {'method': lambda x: x.get_k_stat(3).val, 'true': mu_delta},
                   'k4': {'method': lambda x: x.get_k_stat(4).val,
                          'true': 2 * mu_bar},
                   'c4/c2': {'method': lambda x: x.get_cumulant(4).val / x.get_cumulant(2).val,
                             'true': 1},
                   'k4/k2': {'method': lambda x: x.get_k_stat(4).val / x.get_k_stat(2).val,
                             'true': 1}
                   }

    save_path = '/home/dylan/Desktop/'

    # demo_plots(dist)
    start = time.time()
    event_means, event_errs, event_percs = simulate(dist, num_events, trials, binning, moment_pars, percentiles)
    print(f'Simulation time: {time.time() - start}s')
    plot_moments(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)

    print('donzo')


def demo_plots(dist):
    # plot_pmf(dist)
    # plot_pmf_sample(dist, 100)
    # plot_pmf_sample(dist, 1000)
    plot_c2_vs_n(dist, 10, 1000)
    plot_k2_vs_n(dist, 10, 1000)


def simulate(dist, num_events, trials, binning, moment_pars, percentiles):
    event_means = {x: [] for x in moment_pars.keys()}
    event_errs = {x: [] for x in moment_pars.keys()}
    event_percs = {x: {y: [] for y in percentiles} for x in moment_pars.keys()}
    for events in num_events:
        print(f'Events: {events}')
        trial_moments = {x: [] for x in moment_pars.keys()}
        for trial in range(trials):
            hist, bin_edges = np.histogram(dist.rvs(size=events), binning)
            hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
            for moment, moment_val in calc_moments(hist, moment_pars).items():
                trial_moments[moment].append(moment_val)
        trial_mean, trial_err, trial_perc = comb_trials(trial_moments, percentiles)
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


def calc_moments(hist, moment_pars):
    moments = {}
    stats = DistStats(hist)
    for moment, pars in moment_pars.items():
        moments.update({moment: pars['method'](stats)})

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
    for i in range(2, 5):
        c, k = 'c' + str(i), 'k' + str(i)
        fig, ax = plt.subplots()
        fig.canvas.manager.set_window_title('Order ' + str(i))

        ax.axhline(moment_pars[k]['true'], color='r', ls='--')

        for j in range(int(len(percs) / 2)):
            ax.fill_between(num_events, event_percs[c][percs[j]], event_percs[c][percs[len(percs)-1-j]],
                            color='b', alpha=0.3)
            ax.fill_between(num_events, event_percs[k][percs[j]], event_percs[k][percs[len(percs)-1-j]],
                            color='g', alpha=0.3)

        ax.fill_between(num_events, event_means[c] + event_errs[c], event_means[c] - event_errs[c],
                        color='b', alpha=0.3)
        ax.fill_between(num_events, event_means[k] + event_errs[k], event_means[k] - event_errs[k],
                        color='g', alpha=0.3)

        ax.plot(num_events, event_means[c], color='b')
        ax.plot(num_events, event_means[k], color='g')
    plt.show()


def plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percs):
    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title('C4 / C2')

    ax.axhline(moment_pars['k4/k2']['true'], color='r', ls='--')
    ax.fill_between(num_events, event_means['c4/c2'] - event_errs['c4/c2'], event_means['c4/c2'] + event_errs['c4/c2'],
                    color='b', alpha=0.3)
    ax.fill_between(num_events, event_means['k4/k2'] - event_errs['k4/k2'], event_means['k4/k2'] + event_errs['k4/k2'],
                    color='g', alpha=0.3)

    ax.plot(num_events, event_means['c4/c2'], color='b')
    ax.plot(num_events, event_means['k4/k2'], color='g')

    plt.show()


def plot_pmf(dist):
    fig, ax = plt.subplots()
    x_plot = np.arange(dist.ppf(0.001), dist.ppf(0.999))
    ax.plot(x_plot, dist.pmf(x_plot), color='red', lw=4,
            label=fr'Skellam PDF $\mu_1$={dist.args[0]}, $\mu_2$={dist.args[1]}')
    ax.set_xlabel('Distribution Value')
    ax.set_ylabel('Probability Mass')
    ax.legend()
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


def plot_c2_vs_n(dist, num, trials):
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
    ax.axhspan(mean+sem, mean-sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
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


def plot_k2_vs_n(dist, num, trials):
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
    ax.axhspan(mean+sem, mean-sem, color='green', alpha=0.3, label=rf'Mean of Trials: {mean:.4f}$\pm${sem:.4f}')
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

    plt.show()


if __name__ == '__main__':
    main()
