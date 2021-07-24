#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 21 8:36 PM 2021
Created in PyCharm
Created as QGP_Scripts/cbwc_data_emulator

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
import ROOT as r


def main():
    plt.rcParams['figure.autolayout'] = True

    threads = 15

    cent_edges = [6, 9, 12, 17, 24, 32, 42, 54, 69, 86, 106, 129, 156, 188, 226, 271]
    percentiles = []

    mu1, mu2 = 20.9, 0.2
    mu_bar = (mu1 + mu2) / 2
    mu_delta = mu1 - mu2
    dist = skellam(mu1, mu2)
    binning = np.arange(-0.5, mu_bar * 2 * 10 + 1.5, 1)
    trials = 10000
    moment_pars = {'c2': {'method': get_c2, 'true': 2 * mu_bar},
                   'c3': {'method': get_c3, 'true': mu_delta},
                   'c4': {'method': get_c4, 'true': 2 * mu_bar},
                   'c5': {'method': get_c5, 'true': mu_delta},
                   'c6': {'method': get_c6, 'true': 2 * mu_bar},
                   'k2': {'method': get_k2, 'true': 2 * mu_bar},
                   'k3': {'method': get_k3, 'true': mu_delta},
                   'k4': {'method': get_k4, 'true': 2 * mu_bar},
                   'k5': {'method': get_k5, 'true': mu_delta},
                   'k6': {'method': get_k6, 'true': 2 * mu_bar},
                   'c4/c2': {'method': get_c4_div_c2, 'true': 1},
                   'k4/k2': {'method': get_k4_div_k2, 'true': 1},
                   'c6/c2': {'method': get_c6_div_c2, 'true': 1},
                   'k6/k2': {'method': get_k6_div_k2, 'true': 1},
                   'c4/c2 - k4/k2': {'method': get_c4_div_c2_sub_k4_div_k2, 'true': 1},
                   'c6/c2 - k6/k2': {'method': get_c6_div_c2_sub_k6_div_k2, 'true': 1},
                   }

    save_path = '/home/dylan/Desktop/'

    emulate_data(cent_edges, dist, binning, trials, moment_pars, percentiles, threads)

    print('donzo')


def emulate_data(cent_edges, dist, binning, trials, moment_pars, percentiles, threads):
    ref3, events = get_events()
    plot_ref3_dist(ref3, events, cent_edges)

    start = time.time()
    simulate(dist, np.asarray(events), trials, binning, moment_pars, percentiles, threads)
    print(f'Simulation time: {time.time() - start}s')


def plot_ref3_dist(ref3, events, cent_edges):
    plt.plot(ref3, events)
    for edge in cent_edges:
        plt.axvline(edge, color='r', ls='--')
    plt.grid()
    plt.show()


def get_events():
    f = r.TFile('/home/dylan/Research/Data/default/rapid05_n1ratios_dca1_nsprx1_m2r6_m2s0_nhfit20_0/7GeV/QA_7GeV.root',
                'READ')
    h = f.Get('pre_refn_rapid05_n1ratios_dca1_nsprx1_m2r6_m2s0_nhfit20_0_7').Clone()
    # h.Draw()
    ref3 = []
    events = []
    # h.GetXaxis().GetFirst()
    for bini in range(270, 406):
        # print(f'bin: {bini}, bin_center: {h.GetBinCenter(bini)}, bin_content: {h.GetBinContent(bini)}')
        if int(h.GetBinContent(bini)) > 10:
            ref3.append(int(h.GetBinCenter(bini)))
            events.append(int(h.GetBinContent(bini)))
    f.Close()

    return ref3, events


def simulate(dist, events, trials, binning, moment_pars, percentiles, threads):
    trial_mean, trial_err, trial_perc = sim_centrality(dist, trials, moment_pars, binning, percentiles, events, threads)

    print("Plotting...")

    for i in range(2, 7):
        fig, ax = plt.subplots()
        ax.errorbar(['c', 'k'], [trial_mean[f'c{i}'], trial_mean[f'k{i}']], marker='o', ls='none',
                    yerr=[trial_err[f'c{i}'], trial_err[f'k{i}']])
        ax.axhline(moment_pars[f'k{i}']['true'], ls='--', color='red')
    plt.show()


def sim_centrality(dist, trials, moment_pars, binning, percentiles, events, threads):
    print(f'Simulating Centrality')
    trial_cbwc_moments = {x: [] for x in moment_pars.keys()}

    # [0] trial number [1] independent state for each trial
    rand_states = list(zip(range(trials), [np.random.RandomState() for i in range(trials)]))
    func = partial(sim_cent_trial, dist, moment_pars, binning, events)
    with Pool(threads) as p:
        trial_cbwc_stats = p.map(func, rand_states)

    for trial_moments in trial_cbwc_stats:
        for moment, cbwc_val in trial_moments.items():
            trial_cbwc_moments[moment].append(cbwc_val)

    return comb_trials(trial_cbwc_moments, percentiles)


def sim_cent_trial(dist, moment_pars, binning, events, state):
    print(f'Trial #{state[0]}')
    cbwc_moments = {x: 0 for x in moment_pars.keys()}
    event_sum = sum(events)
    for n in events:
        hist, bin_edges = np.histogram(dist.rvs(size=n, random_state=state[1]), binning)
        hist = dict(zip((bin_edges[1:] + bin_edges[:-1]) / 2.0, hist))
        for moment, moment_val in calc_moments(hist, moment_pars).items():
            cbwc_moments[moment] += moment_val * n / event_sum

    return cbwc_moments


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
    fig = plt.figure(figsize=(5, 1))
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


# def plot_moments_together(num_events, event_means, event_errs, event_percs, moment_pars, percs):
#     fig, axs = plt.subplots(5, 1, sharex=True)
#     plt.subplots_adjust(hspace=0.0, wspace=0.0)
#     fig.canvas.manager.set_window_title('Cumulant KStat Compare')
#     for i in range(2, 7):
#         c, k = 'c' + str(i), 'k' + str(i)
#
#         if i == 2:
#             axs[i-2].axhline(moment_pars[k]['true'], color='r', ls='--', label='Analytical')
#         else:
#             axs[i - 2].axhline(moment_pars[k]['true'], color='r', ls='--')
#
#         for j in range(int(len(percs) / 2)):
#             axs[i-2].fill_between(num_events, event_percs[c][percs[j]], event_percs[c][percs[len(percs) - 1 - j]],
#                             color='b', alpha=0.3)
#             axs[i-2].fill_between(num_events, event_percs[k][percs[j]], event_percs[k][percs[len(percs) - 1 - j]],
#                             color='g', alpha=0.3)
#
#         axs[i-2].fill_between(num_events, event_means[c] + event_errs[c], event_means[c] - event_errs[c],
#                         color='b', alpha=0.3)
#         axs[i-2].fill_between(num_events, event_means[k] + event_errs[k], event_means[k] - event_errs[k],
#                         color='g', alpha=0.3)
#
#         axs[i-2].plot(num_events, event_means[c], color='b', label=f'{c}')
#         axs[i-2].plot(num_events, event_means[k], color='g', label=f'{k}')
#         axs[i-2].legend()
#     axs[2].set_xlabel('Sample Size n')
#     plt.show()


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
    fig1.canvas.manager.set_window_title('C4 / C2 with K4 / K2')

    axs1[0].axhline(moment_pars['k4/k2']['true'], color='r', ls='--', label='Analytical')
    axs1[0].fill_between(num_events, event_means['c4/c2'] - event_errs['c4/c2'],
                         event_means['c4/c2'] + event_errs['c4/c2'],
                         color='b', alpha=0.3)
    axs1[0].fill_between(num_events, event_means['k4/k2'] - event_errs['k4/k2'],
                         event_means['k4/k2'] + event_errs['k4/k2'],
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
    fig2.canvas.manager.set_window_title('(C4 / C2) - (K4 / K2)')
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


if __name__ == '__main__':
    main()
