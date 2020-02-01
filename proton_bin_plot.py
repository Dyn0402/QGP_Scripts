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
from scipy.stats import binom


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
    ax1.bar(x, y, align='center', zorder=0, label="Data")
    ax1.scatter(x, y_binom, color='red', label="Binomial")
    ax1.legend()
    y_diff = y - y_binom
    ax2.axhline(0, color='red', ls='--')
    ax2.errorbar(x, y_diff, yerr=y_err, fmt='bo')
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
    plt.title('Protons in Event vs Protons in Bin')
    plt.legend()
    plt.show()


def plot_ratio_data(data, x_range, y_range, divs, x_label='Number of Protons in Bin', y_label='Number of Protons in Event'):
    x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
                       np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(1.0 / divs, color='red')
    plt.axvline(1.0, color='blue')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title('Protons in Event vs Ratios')
    plt.show()


def plot_diff_data(data, x_range, y_range, divs, x_label='Number of Protons in Bin', y_label='Number of Protons in Event'):
    x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
                       np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(0.0, color='red')
    y2 = np.asarray(y_range) / (1.0 - 1.0 / divs)
    y3 = -float(divs) * np.asarray(y_range)
    plt.plot(y_range, y2, color='blue')
    plt.plot(y_range, y3, color='blue')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, 40])
    plt.title('Protons in Event vs Differences')
    plt.show()


def ratio_transform(data, divs):
    ratio_data = np.zeros((20, 40))
    for num_in_bin, bin_data in enumerate(data):
        for total_protons, num_events in enumerate(bin_data):
            if total_protons > 0 and num_events > 0:
                ratio = float(num_in_bin) / total_protons
                ratio_bin = int(ratio * 19)
                ratio_data[ratio_bin][total_protons] += num_events

    plot_ratio_data(ratio_data, range(0, 40), np.linspace(0, 1.05, 21), divs, x_label='Ratio')


def diff_transform(data, divs):
    diff_data = np.zeros((20, 40))
    for num_in_bin, bin_data in enumerate(data):
        for total_protons, num_events in enumerate(bin_data):
            if total_protons > 0 and num_events > 0:
                diff = float(num_in_bin) - float(total_protons) / divs
                diff_norm = (float(diff) + 40.0 / divs) / (2 * 40.0 / divs)
                diff_bin = int(diff_norm * 19 + 0.5)
                diff_data[diff_bin][total_protons] += num_events

    plot_diff_data(diff_data, range(0, 40), np.linspace(-40.0 / divs, 40.0 / divs + 40.0 / divs / 20, 21),
                   divs, x_label='Difference')


if __name__ == '__main__':
    main()
