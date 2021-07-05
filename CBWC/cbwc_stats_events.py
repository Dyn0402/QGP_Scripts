#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 28 9:16 PM 2021
Created in PyCharm
Created as QGP_Scripts/cbwc_stats_events

@author: Dylan Neff, Dylan
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
from DistStats import DistStats
import seaborn as sns


def main():
    num_events = np.asarray(np.arange(0, 25, 1))
    binom_n = 20
    binom_p = 0.4
    trials = 10
    # save_path = 'N:/UCLA_Microsoft/OneDrive - personalmicrosoftsoftware.ucla.edu/Research/UCLA/Presentations/1-29-21/'
    save_path = '/home/dylan/Desktop/'

    trial_moments = []
    for trial in range(trials):
        hists = simulate(num_events, binom_n, binom_p)
        # plot_dists(hists)
        trial_moment = calc_moments(hists)
        trial_moments.append(trial_moment)

    # plot_mom_trials(num_events, trial_moments, binom_n, binom_p)

    mom_means, mom_errs = comb_trials(trial_moments)
    plot_moments(num_events, mom_means, mom_errs, trial_moments, binom_n, binom_p, save_path)

    # test_list = [1, 2, 3, 4, 3, 3, 5, 3, 4, 6, 5, 3, 4, 3, 4, 1, 8]
    # hist, bin_edges = np.histogram(test_list, np.arange(-0.5, max(test_list) + 1.5, 1))
    # hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
    #
    # stats = DistStats(hist)
    # print(stats.get_cumulant(2))
    # print(stats.get_sd().val**2)
    # print(np.std(test_list)**2)

    print('donzo')


def simulate(num_events, binom_n, binom_p):
    hists = []
    for events in num_events:
        ref3_binom = binom(binom_n, binom_p)
        hist, bin_edges = np.histogram(ref3_binom.rvs(size=events), np.arange(-0.5, binom_n + 1.5, 1))
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        hists.append(hist)

    return hists


def plot_dists(hists):
    n_plots_y = int(len(hists)**0.5 + 0.5)
    n_plots_x = int(len(hists) / n_plots_y + 0.5)
    fig, axs = plt.subplots(n_plots_y, n_plots_x, sharey=True, sharex=True)

    ref3s = list(hists.keys())
    i = 0
    for y in range(n_plots_y):
        for x in range(n_plots_x):
            if i < len(ref3s):
                norm = sum(hists[ref3s[i]].values())
                axs[y, x].bar(hists[ref3s[i]].keys(), np.asarray(list(hists[ref3s[i]].values())) / norm, 1.0,
                              align='center')
                axs[y, x].text(0, 0.25, norm)
            i += 1

    fig.canvas.manager.window.showMaximized()
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    # plt.show()


def calc_moments(hists):
    moments = {'sd': [], 'skew': [], 'kurt': [], 'c2': [], 'c3': [], 'c4': []}
    for hist in hists:
        stats = DistStats(hist)
        moments['kurt'].append(stats.get_kurtosis().val)
        moments['skew'].append(stats.get_skewness().val)
        moments['sd'].append(stats.get_sd().val)
        moments['c4'].append(stats.get_cumulant(4).val)
        moments['c3'].append(stats.get_cumulant(3).val)
        moments['c2'].append(stats.get_cumulant(2).val)

    return moments


def comb_trials(trial_moments):
    moments = {}
    for mom in trial_moments[0]:
        moments.update({mom: [trial[mom] for trial in trial_moments]})

    trials = len(trial_moments)
    means, errs = {}, {}
    for mom, vals in moments.items():
        means.update({mom: np.mean(vals, axis=0)})
        errs.update({mom: np.std(vals, axis=0) / np.sqrt(trials)})

    return means, errs


def plot_mom_trials(num_events, trial_moments, n, p):
    q = 1 - p
    binom_sd = np.sqrt(n * p * q)
    binom_skew = (q - p) / np.sqrt(n * p * q)
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    binom_c2 = n*p*q
    binom_c3 = n*p*q*(1-2*p)
    binom_c4 = n * p * q * (1 + (3 * n - 6) * p * q) - 3 * binom_c2**2

    mom_bionms = {'sd': binom_sd, 'skew': binom_skew, 'kurt': binom_kurt,
                  'c2': binom_c2, 'c3': binom_c3, 'c4': binom_c4}

    mom_names = {'sd': 'Standard Deviation', 'skew': 'Skewness', 'kurt': 'Kurtosis',
                 'c2': 'C2', 'c3': 'C3', 'c4': 'C4'}

    trials = len(trial_moments)
    moments = {}
    for mom in trial_moments[0]:
        moments.update({mom: [trial[mom] for trial in trial_moments]})

    for mom, values in moments.items():
        values = np.asarray(values).T
        events = []
        mom_vals = []
        for index, val_trials in enumerate(values):
            for val in val_trials:
                events.append(num_events[index])
                mom_vals.append(val)
        fig, ax = plt.subplots()
        fig.canvas.set_window_title(f'{mom_names[mom]} All')
        ax.scatter(events, mom_vals, alpha=0.02, marker='s', label='Simulation')
        # sns.kdeplot(x=events, y=mom_vals, shade=True, cmap='PuBu')
        # sns.scatterplot(x=events, y=mom_vals, alpha=np.exp(-(trials - 1) / 600), label='Simulation')
        ax.plot(num_events, [mom_bionms[mom] for event in num_events], 'r', label='Binomial Expectation')
        ax.legend()
        ax.set_title(f'Binomial {mom_names[mom]} over {trials} Trials')
        ax.set_xlabel('Number of Events')
        ax.set_ylabel(mom_names[mom])
        ax.tight_layout()

    # plt.show()


def plot_moments(num_events, mom_means, mom_errs, trial_moments, n, p, save_path):
    q = 1 - p
    binom_sd = np.sqrt(n * p * q)
    binom_skew = (q - p) / np.sqrt(n * p * q)
    binom_kurt = (1 - 6 * p * q) / (n * p * q)
    binom_c2 = n*p*q
    binom_c3 = n*p*q*(1-2*p)
    binom_c4 = n * p * q * (1 + (3 * n - 6) * p * q) - 3 * binom_c2**2

    mom_bionms = {'sd': binom_sd, 'skew': binom_skew, 'kurt': binom_kurt,
                  'c2': binom_c2, 'c3': binom_c3, 'c4': binom_c4}

    mom_names = {'sd': 'Standard Deviation', 'skew': 'Skewness', 'kurt': 'Kurtosis',
                 'c2': 'C2', 'c3': 'C3', 'c4': 'C4'}

    trials = len(trial_moments)
    moments = {}
    for mom in trial_moments[0]:
        moments.update({mom: [trial[mom] for trial in trial_moments]})

    for mom, values in moments.items():
        values = np.asarray(values).T
        events = []
        mom_vals = []
        for index, val_trials in enumerate(values):
            for val in val_trials:
                events.append(num_events[index])
                mom_vals.append(val)
        
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=(19.25, 7.2))
        fig_name = f'{mom_names[mom]}_{trials}t_{max(num_events)}e'
        fig.canvas.manager.set_window_title(fig_name)

        axs[0].plot(num_events, [mom_bionms[mom] for event in num_events], 'r--')
        axs[0].scatter(events, mom_vals, alpha=0.03, marker='_', s=55, label='Simulation Trials')
        # sns.kdeplot(x=events, y=mom_vals, shade=True, cmap='PuBu')
        # sns.scatterplot(x=events, y=mom_vals, alpha=np.exp(-(trials - 1) / 600), label='Simulation')
        leg0 = axs[0].legend()
        for lh in leg0.legendHandles:
            lh.set_alpha(0.7)
        axs[0].set_title(f'Binomial {mom_names[mom]} over {trials} Trials')
        axs[0].set_ylabel(f'{mom_names[mom]}')

        axs[1].errorbar(num_events, mom_means[mom], yerr=mom_errs[mom], label='Simulation Average')
        axs[1].plot(num_events, [mom_bionms[mom] for event in num_events], 'r--', label='Binomial Expectation')
        axs[1].legend()
        axs[1].set_xlabel('Number of Events')
        axs[1].set_ylabel(f'Average {mom_names[mom]}')

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.03)
        plt.savefig(f'{save_path}{fig_name}.png')


    plt.show()


if __name__ == '__main__':
    main()
