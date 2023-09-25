#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 15 4:56 PM 2021
Created in PyCharm
Created as QGP_Scripts/presentation_plots

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import binom
from fractions import Fraction

from Analyzer.AzimuthBinData import AzimuthBinData


def main():
    divs = 120
    cent = 8
    energy = 39
    samples = 72  # Just for title
    data_name = 'AMPT'  # Just for title
    set_group = 'default_single'
    set_name = 'Ampt_rapid05_single_'
    # set_group = 'default_resample_epbins1'
    # set_name = 'Ampt_rapid05_resample_norotate_epbins1_'
    set_num = 0
    data_set = '_Ampt_New_Coal'
    # base_path = '/home/dylan/Research/Data'
    base_path = 'F:/Research/Data'
    # base_path = 'D:/Transfer/Research/Data'
    set_path = f'{set_group}/{set_name}{set_num}/{energy}GeV/ratios_divisions_{divs}_centrality_{cent}_local.txt'
    path = f'{base_path}{data_set}/{set_path}'
    path_mix = f'{base_path}{data_set}_Mix/{set_path}'
    raw = AzimuthBinData(path=path, div=divs)
    mix = AzimuthBinData(path=path_mix, div=divs)

    # title_sufx = f'\n{energy}GeV, 0-5% Centrality, {divs}° Partitions, {samples} Sample per Event'
    title_sufx = f'\n{energy}GeV, 0-5% Centrality, {divs}° Partitions'

    # plot_2d(raw.get_dist(), raw.max_particle, raw.get_max_bin(), divs, data_name, title_sufx)
    plot_binomial(raw.get_dist(), 20, divs, data_name, title_sufx=title_sufx)
    # plot_ratio(raw.get_dist(), raw.max_particle, divs, x_bins=20, title_sufx=title_sufx)
    # plot_pull(raw.get_dist(), raw.max_particle, divs, x_bins=40, title_sufx=title_sufx)

    print('donzo')


def plot_2d(file_data, max_particle, max_bin, divs, data_name, title_sufx):
    data = np.zeros((max_particle + 1, max_particle + 1))
    for total_particles, bins in file_data.items():
        for bin_particles, counts in enumerate(bins):
            data[total_particles][bin_particles] += counts

    data = np.ma.masked_where(data <= 0, data)
    edges = np.arange(-0.5, max_particle + 1)
    x_edges, y_edges = np.meshgrid(edges, edges)
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=144)
    plt.pcolormesh(x_edges, y_edges, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.plot(edges, 360 / float(divs) * edges, color='red', label='mean')
    plt.plot(edges, edges, color='blue', label='max')
    plt.title(f'{data_name} Protons in Event vs Protons in Partition{title_sufx}')
    plt.xlabel('Number of Protons in Partition')
    plt.ylabel('Number of Protons in Event')
    plt.xlim(-0.5, max_bin + 2)
    plt.ylim(-0.5, max_particle + 3)
    plt.colorbar()
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.subplots_adjust(right=1)
    plt.show()


def plot_ratio(file_data, max_particle, divs, x_bins=None, title_sufx=''):
    if x_bins is None:
        x_bins = max_particle
    y_edges = np.arange(-0.5, max_particle + 1)
    x_edges = np.linspace(0.0, 1.0 + 0.01 / x_bins, x_bins + 1)
    x_range = x_edges[-1] - x_edges[0]

    ratio2d, ratio1d = ratio_transform(file_data, max_particle, x_range, x_bins)

    plot_ratio2d(ratio2d, x_edges, y_edges, max_particle, divs, title_sufx)
    plot_ratio1d(ratio1d, x_edges, divs, title_sufx)


def ratio_transform(file_data, max_particle, x_range, x_bins):
    ratio2d = np.zeros((max_particle + 1, x_bins))
    ratio1d = np.zeros(x_bins)
    for total_particles, bins in file_data.items():
        for bin_particles, counts in enumerate(bins):
            ratio_val = Fraction(bin_particles, total_particles)
            ratio_bin = int(ratio_val * x_bins / x_range)
            ratio2d[total_particles][ratio_bin] += counts
            ratio1d[ratio_bin] += counts

    return ratio2d, ratio1d


def plot_ratio2d(ratio, x_edges, y_edges, max_particle, divs, title_sufx):
    plt.figure(figsize=(6.4735, 4))
    ratio = np.ma.masked_where(ratio <= 0, ratio)
    x_medges, y_medges = np.meshgrid(x_edges, y_edges)
    plt.pcolormesh(x_medges, y_medges, ratio, norm=colors.LogNorm(vmin=ratio.min(), vmax=ratio.max()))
    plt.axvline(divs / 360, color='red', label='mean')
    plt.axvline(1, color='blue', label='max')
    plt.title(f'Particles in Event vs Ratio{title_sufx}')
    plt.xlabel('Ratio')
    plt.ylabel('Number of Particles in Event')
    # plt.xlim(x_edges[], max_bin + 2)
    plt.ylim(-0.5, max_particle + 3)
    plt.colorbar()
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_ratio1d(ratio, x_edges, divs, title_sufx):
    plt.figure(figsize=(6.4735, 4))
    x_centers = (x_edges[1:] + x_edges[:-1]) / 2
    x_widths = x_edges[1:] - x_edges[:-1]
    plt.bar(x_centers, ratio, width=x_widths, align='center')
    plt.axvline(float(divs) / 360, color='red', label='mean')
    plt.axvline(1, color='blue', label='max')
    plt.yscale('log')
    plt.xlabel('Ratio')
    plt.ylabel('Counts')
    plt.title(f'Ratio Distribution{title_sufx}')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_pull(file_data, max_particle, divs, x_bins=None, title_sufx=''):
    p = Fraction(divs, 360)

    if x_bins is None:
        x_bins = max_particle
    y_edges = np.arange(-0.5, max_particle + 1)
    diff_edges = np.linspace(-p * max_particle, (1 - p) * max_particle + 0.01 * max_particle / x_bins, x_bins + 1)
    pull_edges = np.linspace(-(p * max_particle / (1 - p)) ** 0.5,
                             (max_particle / (p * (1 - p))) ** 0.5 * (1 - p) + 0.01 * max_particle / x_bins, x_bins + 1)

    diff2d, pull2d, diff1d, pull1d = pull_transform(file_data, max_particle, divs, diff_edges, pull_edges, x_bins)

    diff_plt_range = get_min_max(diff1d, diff_edges)
    pull_plt_range = get_min_max(pull1d, pull_edges)

    plot_diff2d(diff2d, diff_edges, y_edges, max_particle, divs, diff_plt_range, title_sufx)
    plot_pull2d(pull2d, pull_edges, y_edges, max_particle, divs, pull_plt_range, title_sufx)
    plot_diff1d(diff1d, diff_edges, diff_plt_range, title_sufx)
    plot_pull1d(pull1d, pull_edges, pull_plt_range, title_sufx)


def pull_transform(file_data, max_particle, divs, diff_edges, pull_edges, x_bins):
    diff2d = np.zeros((max_particle + 1, x_bins))
    pull2d = np.zeros((max_particle + 1, x_bins))
    diff1d = np.zeros(x_bins)
    pull1d = np.zeros(x_bins)

    diff_range = diff_edges[-1] - diff_edges[0]
    pull_range = pull_edges[-1] - pull_edges[0]

    for total_particles, bins in file_data.items():
        for bin_particles, counts in enumerate(bins):
            p = Fraction(divs, 360)
            diff_val = bin_particles - p * total_particles
            pull_val = diff_val / (total_particles * p * (1 - p)) ** 0.5
            diff_bin = int((diff_val - diff_edges[0]) * x_bins / diff_range)
            pull_bin = int((pull_val - pull_edges[0]) * x_bins / pull_range)

            diff2d[total_particles][diff_bin] += counts
            pull2d[total_particles][pull_bin] += counts
            diff1d[diff_bin] += counts
            pull1d[pull_bin] += counts

    return diff2d, pull2d, diff1d, pull1d


def plot_diff2d(diff, x_edges, y_edges, max_particle, divs, plt_range=None, title_sufx=''):
    diff = np.ma.masked_where(diff <= 0, diff)
    p = Fraction(divs, 360)
    x_medges, y_medges = np.meshgrid(x_edges, y_edges)
    plt.pcolormesh(x_medges, y_medges, diff, norm=colors.LogNorm(vmin=diff.min(), vmax=diff.max()))
    plt.axvline(0, color='red', label='mean')
    plt.plot(x_edges, -x_edges / p, color='black', label='min')
    plt.plot(x_edges, x_edges / (1 - p), color='blue', label='max')
    plt.title(f'Particles in Event vs Difference{title_sufx}')
    plt.xlabel('Difference')
    plt.ylabel('Number of Particles in Event')
    if plt_range is not None:
        plt.xlim(plt_range)
    plt.ylim(-0.5, max_particle + 3)
    plt.colorbar()
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()


def plot_diff1d(diff, x_edges, plt_range=None, title_sufx=''):
    x_centers = (x_edges[1:] + x_edges[:-1]) / 2
    x_widths = x_edges[1:] - x_edges[:-1]
    plt.bar(x_centers, diff, width=x_widths, align='center')
    plt.axvline(0, color='red', label='mean')
    if plt_range is not None:
        plt.xlim(plt_range)
    plt.yscale('log')
    plt.xlabel('Difference')
    plt.ylabel('Counts')
    plt.title(f'Difference Distribution{title_sufx}')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_pull2d(pull, x_edges, y_edges, max_particle, divs, plt_range=None, title_sufx=''):
    pull = np.ma.masked_where(pull <= 0, pull)
    p = Fraction(divs, 360)
    x_medges, y_medges = np.meshgrid(x_edges, y_edges)
    plt.pcolormesh(x_medges, y_medges, pull, norm=colors.LogNorm(vmin=pull.min(), vmax=pull.max()))
    plt.axvline(0, color='red', label='mean')
    x_pos = np.linspace(0, x_edges[-1], 100)
    x_neg = np.linspace(x_edges[0], 0, 100)
    plt.plot(x_neg, x_neg ** 2 * (1 - p) / p, color='black', label='min')
    plt.plot(x_pos, x_pos ** 2 * p * (1 - p) / (1 - p) ** 2, color='blue', label='max')
    plt.title(f'Particles in Event vs Pull{title_sufx}')
    plt.xlabel('Pull')
    plt.ylabel('Number of Particles in Event')
    if plt_range is not None:
        plt.xlim(plt_range)
    plt.ylim(-0.5, max_particle + 3)
    plt.colorbar()
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()


def plot_pull1d(pull, x_edges, plt_range=None, title_sufx=''):
    x_centers = (x_edges[1:] + x_edges[:-1]) / 2
    x_widths = x_edges[1:] - x_edges[:-1]
    plt.bar(x_centers, pull, width=x_widths, align='center')
    plt.axvline(0, color='red', label='mean')
    if plt_range is not None:
        plt.xlim(plt_range)
    plt.yscale('log')
    plt.xlabel('Pull')
    plt.ylabel('Counts')
    plt.title(f'Pull Distribution{title_sufx}')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_binomial(data, particles, divs, data_name, title_sufx=''):
    y = np.zeros(shape=particles + 1)
    y[:len(data[particles])] = data[particles]
    x = range(particles + 1)
    y_binom = sum(y) * binom.pmf(x, particles, float(divs) / 360)
    # y_err = np.sqrt(y)  # No longer a good approximation with multi-sampling

    fig0, ax0 = plt.subplots(figsize=(6.6, 3.3), dpi=144)
    ax0.bar(x, y, width=0.9, align='center', zorder=0, label=f'{particles} Proton Events')
    ax0.set_xticks(range(0, len(y), 2))
    ax0.set_title(f'{data_name} Protons in {divs}° Partition for {particles} Proton Events' + title_sufx)
    ax0.set_xlabel('Number of Protons in Azimuthal Partition')
    ax0.set_ylabel('Partitions')
    # ax0.set_xticks(x)
    ax0.set_xlim([-0.5, particles + 0.5])
    ax0.legend()
    fig0.tight_layout()

    fig1, ax1 = plt.subplots(figsize=(6.6, 3.3), dpi=144)
    # fig1, ax1 = plt.subplots(figsize=(6, 2.8), dpi=144)
    ax1.bar(x, y, align='center', zorder=0, label=f'{particles} Proton Events')
    ax1.scatter(x, y_binom, color='red', label='Binomial Distribution')
    ax1.set_xticks(range(0, len(y), 2))
    ax1.set_title(f'{data_name} Protons in {divs}° Partition vs Binomial for {particles} Proton Events' + title_sufx)
    ax1.set_xlabel('Number of Protons in Azimuthal Partition')
    ax1.set_ylabel('Partitions')
    ax1.set_xlim([-0.5, particles + 0.5])
    ax1.legend()
    fig1.tight_layout()

    fig2, ax2 = plt.subplots()
    y_diff = y - y_binom
    ax2.axhline(0, color='red', ls='--')
    ax2.errorbar(x, y_diff, fmt='bo')
    ax2.set_xticks(range(0, len(y), 2))
    ax2.set_title(f'Protons in {divs}° Bin Minus Binomial for {particles} Proton Events' + title_sufx)
    ax2.set_xlabel('Number of Protons in Bin')
    ax2.set_ylabel('Data Events Minus Binomial')
    ax2.set_xlim([-0.5, particles + 0.5])

    fig3, ax3 = plt.subplots()
    y_ratio = y / y_binom
    # y_ratio_err = y_err / y_binom
    ax3.axhline(1, color='red', ls='--')
    ax3.errorbar(x, y_ratio, fmt='bo')
    ax3.set_xticks(range(0, len(y), 2))
    ax3.set_title(f'Protons in {divs}° Bin Divided by Binomial for {particles} Proton Events' + title_sufx)
    ax3.set_xlabel('Number of Protons in Bin')
    ax3.set_ylabel('Data Events Divided by Binomial')
    ax3.set_xlim([-0.5, particles + 0.5])

    print(f'Binomial Difference Sum: {sum(y_diff)}')
    plt.show()


def plot_azbin_data_proj(data, x_lim, y_lim, divs, x_label='Number of Protons in Bin',
                         y_label='Number of Protons in Event', title_sufx=''):
    """Copied from old script, needs tweaks to data format to work"""
    left, width = 0.1, 0.55
    bottom, height = 0.1, 0.55
    bottom_h = left_h = left + width + 0.08
    rect_2d = [left, bottom, width, height]
    rect_x = [left, bottom_h, width, 0.2]
    rect_y = [left_h, bottom, 0.2, height]

    ax_2d = plt.axes(rect_2d)
    ax_x = plt.axes(rect_x)
    ax_y = plt.axes(rect_y)

    ax_2d.set_xlabel('x axis')
    ax_2d.set_ylabel('y axis')
    ax_x.set_title('x projection')
    ax_y.set_title('y projection')

    # x, y = np.meshgrid(x_range, y_range)
    x_range = range(0, len(data))
    y_range = range(0, len(data[0]))
    x, y = np.meshgrid(np.asarray(x_range) - float(x_range[1] - x_range[0]) / 2,
                       np.asarray(y_range) - float(y_range[1] - y_range[0]) / 2)
    data = np.ma.masked_where(data <= 0, data)
    plot_2d = ax_2d.pcolormesh(x, y, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    y2 = 360 / float(divs) * np.asarray(y_range)
    y3 = 1.0 * np.asarray(y_range)
    ax_2d.plot(y_range, y2, color='red', label='mean')
    ax_2d.plot(y_range, y3, color='blue', label='max')
    ax_2d.set_xlabel(x_label)
    ax_2d.set_ylabel(y_label)
    ax_2d.set_ylim(y_lim)
    ax_2d.set_xlim(x_lim)
    ax_2d.set_xticks(range(x_lim[0], x_lim[1], 2))
    # ax_2d.set_title('Particles in Event vs Particles in Bin'+title_sufx)
    ax_2d.legend(loc='lower right')

    ax_x.set_xlim(x_lim)
    ax_y.set_ylim(y_lim)

    x_proj = np.sum(data, axis=0)
    y_proj = np.sum(data, axis=1)

    ax_x.bar(x_range, x_proj, width=1, align='center')
    ax_y.barh(y_range, y_proj, height=1, align='center')

    cbar_ax = plt.axes([left + width + 0.03, bottom_h - 0.035, 0.025, 0.29])
    plt.colorbar(plot_2d, cax=cbar_ax)

    # title_ax = plt.axes([left_h, bottom_h, 0.15, 0.2])
    plt.text(0.77, 0.85, 'Protons in Event vs Protons in Partition' + title_sufx, fontsize=10,
             transform=plt.gcf().transFigure, wrap=True, va='center', fontweight='bold')

    plt.show()


def get_min_max(dist_1d, dist_edges):
    index = 0
    while index < len(dist_1d) and dist_1d[index] <= 0:
        index += 1
    min_edge = dist_edges[index]
    index = len(dist_1d) - 1
    while index >= 0 and dist_1d[index] <= 0:
        index -= 1
    max_edge = dist_edges[index + 1]

    return min_edge, max_edge


if __name__ == '__main__':
    main()
