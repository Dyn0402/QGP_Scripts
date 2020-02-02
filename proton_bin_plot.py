#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 30 11:32 PM 2020
Created in PyCharm
Created as QGP_Scripts/proton_bin_plot.py

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from scipy.stats import binom
from fractions import Fraction


def main():
    path = '/home/dylan/Research/Data/Single_Ratio0/7GeV/ratios_divisions_3_centrality_8_local.txt'
    divs = 3
    data = read_azbin_data(path)
    # plot_azbin_data(data, range(0, 40), range(0, 20), divs)
    # plot_azbin_data_trans(data, range(0, 40), range(0, 20), divs)
    # ratio_transform(data, divs)
    # diff_transform(data, divs)
    plot_binomial(data, 13, divs)
    print('donzo')


def read_azbin_data(path):
    azbin_data = np.zeros((20, 40))
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split('\t')
            total_protons = int(line[0])
            bin_protons = line[1].split(' ')
            for entry in bin_protons:
                entry = entry.split(':')
                azbin_data[int(entry[0]), total_protons] = int(entry[1])

    return azbin_data


def plot_binomial(data, protons, divs):
    y = np.asarray([ele[protons] for ele in data])
    x = range(len(y))
    y_binom = sum(y)*binom.pmf(x, protons, 1/divs)
    y_err = np.sqrt(y)
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    ax1.bar(x, y, align='center', zorder=0, label=f'{protons} Proton Events')
    ax1.scatter(x, y_binom, color='red', label='Binomial Distribution')
    ax1.set_xticks(range(0, 20, 2))
    ax1.set_title(f'Protons in 3 Division Bin vs Binomial for {protons} Proton Events')
    ax1.set_xlabel('Number of Protons in Bin')
    ax1.set_ylabel('Events')
    ax1.legend()
    y_diff = y - y_binom
    ax2.axhline(0, color='red', ls='--')
    ax2.errorbar(x, y_diff, yerr=y_err, fmt='bo')
    ax2.set_xticks(range(0, 20, 2))
    ax2.set_title(f'Protons in 3 Division Bin Minus Binomial for {protons} Proton Events')
    ax2.set_xlabel('Number of Protons in Bin')
    ax2.set_ylabel('Data Events Minus Binomial')
    print(f'Binomial Difference Sum: {sum(y_diff)}')
    plt.show()


def plot_azbin_data(data, x_range, y_range, divs, x_label='Number of Protons in Bin', y_label='Number of Protons in Event'):
    # x, y = np.meshgrid(x_range, y_range)
    x, y = np.meshgrid(np.asarray(x_range) - float(x_range[1] - x_range[0]) / 2,
                       np.asarray(y_range) - float(y_range[1] - y_range[0]) / 2)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(x, y, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    y2 = 1.0/divs * np.asarray(x_range)
    y3 = 1.0 * np.asarray(x_range)
    plt.plot(x_range, y2, color='red', label='mean')
    plt.plot(x_range, y3, color='blue', label='max')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, 20])
    plt.legend()
    plt.show()


def plot_azbin_data_trans(data, x_range, y_range, divs, x_label='Number of Protons in Bin',
                    y_label='Number of Protons in Event'):
    # x, y = np.meshgrid(x_range, y_range)
    x, y = np.meshgrid(np.asarray(x_range) - float(x_range[1] - x_range[0]) / 2,
                       np.asarray(y_range) - float(y_range[1] - y_range[0]) / 2)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    y2 = float(divs) * np.asarray(y_range)
    y3 = 1.0 * np.asarray(y_range)
    plt.plot(y_range, y2, color='red', label='mean')
    plt.plot(y_range, y3, color='blue', label='max')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, 40])
    plt.xticks(range(0, 20, 2))
    plt.title('Protons in Event vs Protons in Bin')
    plt.legend(loc='lower right')
    plt.show()


def plot_ratio_data(data, x_range, y_range, divs, x_label='Number of Protons in Bin', y_label='Number of Protons in Event'):
    # x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
    #                    np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    x, y = np.meshgrid(x_range, y_range)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(1.0 / divs, color='red', label='mean')
    plt.axvline(1.05, color='blue', label='max')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title('Protons in Event vs Ratios')
    plt.legend(loc='upper right')
    plt.show()


def plot_ratio_dist(y, x, divs):
    plt.bar(x, y, width=x[1]-x[0], align='edge', zorder=0)
    plt.yscale('log')
    plt.axvline(1.0/divs, color='red', label='mean')
    plt.axvline(1.05, color='blue', label='max')
    plt.xlabel('Ratio')
    plt.ylabel('Events')
    plt.title('Ratio Distribution')
    plt.legend(loc='upper right')
    plt.show()


def plot_ratio_kde(data):
    sns.set()
    ax = sns.distplot(data, rug=True, hist=False)
    plt.show()


def plot_diff_data(data, x_range, y_range, divs, x_label='Number of Protons in Bin', y_label='Number of Protons in Event'):
    # x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
    #                    np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    x, y = np.meshgrid(x_range, y_range)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(0.0, color='red', label='mean')
    y2 = np.asarray(y_range) * (divs - 1) / divs
    y3 = -float(divs) * np.asarray(y_range)
    plt.plot(y_range, y2, color='blue', label='max')
    plt.plot(y_range, y3, color='black', label='min')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, 40])
    plt.xlim([-10, 12])
    plt.xticks(range(-10, 12, 2))
    plt.title('Protons in Event vs Differences')
    plt.legend(loc='lower left')
    plt.show()


def plot_diff_dist(y, x, divs):
    plt.bar(x, y, width=x[1]-x[0], align='edge', zorder=0)
    plt.yscale('log')
    plt.axvline(0.0, color='red', label='mean')
    plt.xlabel('Difference')
    plt.ylabel('Events')
    plt.title('Difference Distribution')
    plt.xlim([-10, 12])
    plt.xticks(range(-10, 12, 2))
    plt.legend(loc='upper right')
    plt.show()


def plot_diff_kde(data):
    sns.set()
    ax = sns.distplot(data, rug=True, hist=False)
    plt.show()


def ratio_transform(data, divs):
    ratio_data = np.zeros((21, 40))
    ratio_dist = np.zeros(21)
    ratio_vals = []
    for num_in_bin, bin_data in enumerate(data):
        for total_protons, num_events in enumerate(bin_data):
            if total_protons > 0 and num_events > 0:
                ratio = float(num_in_bin) / total_protons
                ratio_bin = int(ratio * 20)
                ratio_data[ratio_bin][total_protons] += num_events
                ratio_dist[ratio_bin] += num_events
                # ratio_vals.extend([ratio] * int(num_events))

    # plot_ratio_data(ratio_data, range(0, 40), np.linspace(0, 1.05, 22), divs, x_label='Ratio')
    plot_ratio_dist(ratio_dist, np.linspace(0, 1, 21), divs)
    # plot_ratio_kde(ratio_vals)


def diff_transform(data, divs):
    x_bins = 32
    diff_data = np.zeros((x_bins, 40))
    diff_dist = np.zeros(x_bins)
    diff_values = []
    for num_in_bin, bin_data in enumerate(data):
        for total_protons, num_events in enumerate(bin_data):
            if total_protons >= num_in_bin and total_protons > 0 and num_events > 0:
                # diff = float(num_in_bin) - float(total_protons) / divs
                diff = Fraction(num_in_bin * divs - total_protons, divs)
                if diff < Fraction(-total_protons, divs):
                    print('less than min')
                if diff > Fraction(total_protons * (divs - 1), divs):
                    print(f'greater than max: {total_protons * (divs - 1) / divs} | nT: {total_protons} | '
                          f'np: {num_in_bin} | diff: {diff}')
                diff_norm = Fraction((diff + Fraction(40, divs)), 40)
                diff_bin = int(diff_norm * (x_bins - 1))
                diff_data[diff_bin][total_protons] += num_events
                diff_dist[diff_bin] += num_events
                diff_values.extend([diff] * int(num_events))

    plot_diff_data(diff_data, range(0, 40), np.linspace(-40.0 / divs, 40.0 - 40.0 / divs, x_bins),
                   divs, x_label='Difference')
    plot_diff_dist(diff_dist, np.linspace(-40.0 / divs, 40.0 - 40.0 / divs, x_bins), divs)
    plot_diff_kde(diff_values)


if __name__ == '__main__':
    main()
