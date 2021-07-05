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


def main():
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
    dist = skellam(mu1, mu2)
    binning = np.arange(-0.5, mu * 10 + 1.5, 1)
    trials = 1000
    moment_pars = {'c2': {'method': lambda x: x.get_cumulant(2).val, 'true': mu},
                   'c3': {'method': lambda x: x.get_cumulant(3).val, 'true': mu},
                   'c4': {'method': lambda x: x.get_cumulant(4).val,
                          'true': mu},
                   'k2': {'method': lambda x: x.get_k_stat(2).val, 'true': mu},
                   'k3': {'method': lambda x: x.get_k_stat(3).val, 'true': mu},
                   'k4': {'method': lambda x: x.get_k_stat(4).val,
                          'true': mu},
                   'c4/c2': {'method': lambda x: x.get_cumulant(4).val / x.get_cumulant(2).val,
                             'true': 1},
                   'k4/k2': {'method': lambda x: x.get_k_stat(4).val / x.get_k_stat(2).val,
                             'true': 1}
                   }

    save_path = '/home/dylan/Desktop/'

    event_means, event_errs, event_percs = simulate(dist, num_events, trials, binning, moment_pars, percentiles)
    plot_moments(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)
    plot_ratios(num_events, event_means, event_errs, event_percs, moment_pars, percentiles)

    print('donzo')


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
            # print(calc_moments(hist, moment_defs))
            for moment, moment_val in calc_moments(hist, moment_pars).items():
                trial_moments[moment].append(moment_val)
        # print(trial_moments)
        trial_mean, trial_err, trial_perc = comb_trials(trial_moments, percentiles)
        # print(trial_mean)
        for moment, mean in trial_mean.items():
            # print(moment)
            # print(trial_mean[moment])
            event_means[moment].append(mean)
            event_errs[moment].append(trial_err[moment])
            for perc, perc_val in trial_perc[moment].items():
                event_percs[moment][perc].append(perc_val)

    # print(event_means)
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
        # print(f'moment: {moment} | method: {method(stats)} | method_val: {method(stats).val}')
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


if __name__ == '__main__':
    main()
