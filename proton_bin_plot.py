#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 30 11:32 PM 2020
Created in PyCharm
Created as QGP_Scripts/particle_bin_plot.py

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
from scipy.stats import binom
from fractions import Fraction
from scipy.optimize import curve_fit


def main():
    divs = 120
    cent = 8
    energy = 7
    set_name = 'eta05_n1ratios_dca3'
    set_num = 0
    path = f'/home/dylan/Research/Data/{set_name}{set_num}/{energy}GeV/ratios_divisions_{divs}_centrality_{cent}_local.txt'
    path_mix = f'/home/dylan/Research/Data_Mix/{set_name}{set_num}/{energy}GeV/ratios_divisions_{divs}_centrality_{cent}_local.txt'
    title_sufx = f'\n7GeV, 0-5% Centrality, {divs}° Bins'
    data = read_azbin_data(path)
    data_mix = read_azbin_data(path_mix)
    # plot_azbin_data(data, [0, 40], [0, 20], divs)
    plot_azbin_data_trans(data, [0, 22], [0, 45], divs, title_sufx=title_sufx)
    # plot_azbin_data_trans_proj(data, [0, 20], [0, 40], divs, title_sufx=title_sufx)
    # ratio_transform(data, divs, max_particles=46, title_sufx=title_sufx)
    # diff_transform(data, divs, max_particles=46, title_sufx=title_sufx)
    # pull_transform(data, data_mix, divs, max_particles=46, title_sufx=title_sufx)
    plot_binomial(data, 17, divs, title_sufx=title_sufx)
    # plot_data_mixed(data, data_mix, 31, divs, range(10, 26), title_sufx=title_sufx)
    print('donzo')


def read_azbin_data(path):
    azbin_data = np.array([])
    with open(path, 'r') as file:
        lines = file.readlines()
        max_particles = int(lines[-1].split('\t')[0])+1
        azbin_data = np.zeros((max_particles, max_particles))
        for line in lines:
            line = line.strip().split('\t')
            total_particles = int(line[0])
            bin_particles = line[1].split(' ')
            for entry in bin_particles:
                entry = entry.split(':')
                azbin_data[int(entry[0]), total_particles] = int(entry[1])

    return azbin_data


def plot_binomial(data, particles, divs, title_sufx=''):
    y = np.asarray([ele[particles] for ele in data])
    x = range(len(y))
    y_binom = sum(y)*binom.pmf(x, particles, float(divs) / 360)
    y_err = np.sqrt(y)
    fig1, ax1 = plt.subplots()
    ax1.bar(x, y, align='center', zorder=0, label=f'{particles} Particle Events')
    ax1.scatter(x, y_binom, color='red', label='Binomial Distribution')
    ax1.set_xticks(range(0, len(y), 2))
    ax1.set_title(f'Particles in {divs}° Bin vs Binomial for {particles} Particle Events'+title_sufx)
    ax1.set_xlabel('Number of Particles in Bin')
    ax1.set_ylabel('Events')
    ax1.set_xlim([-0.5, particles+0.5])
    ax1.legend()

    fig2, ax2 = plt.subplots()
    y_diff = y - y_binom
    ax2.axhline(0, color='red', ls='--')
    ax2.errorbar(x, y_diff, yerr=y_err, fmt='bo')
    ax2.set_xticks(range(0, len(y), 2))
    ax2.set_title(f'Particles in {divs}° Bin Minus Binomial for {particles} Particle Events'+title_sufx)
    ax2.set_xlabel('Number of Particles in Bin')
    ax2.set_ylabel('Data Events Minus Binomial')
    ax2.set_xlim([-0.5, particles+0.5])

    fig3, ax3 = plt.subplots()
    y_ratio = y / y_binom
    y_ratio_err = y_err / y_binom
    ax3.axhline(1, color='red', ls='--')
    ax3.errorbar(x, y_ratio, yerr=y_ratio_err, fmt='bo')
    ax3.set_xticks(range(0, len(y), 2))
    ax3.set_title(f'Particles in {divs}° Bin Divided by Binomial for {particles} Particle Events' + title_sufx)
    ax3.set_xlabel('Number of Particles in Bin')
    ax3.set_ylabel('Data Events Divided by Binomial')
    ax3.set_xlim([-0.5, particles + 0.5])

    print(f'Binomial Difference Sum: {sum(y_diff)}')
    plt.show()


def plot_data_mixed(data, data_mixed, particles, divs, particle_array=[], title_sufx=''):
    y = np.asarray([ele[particles] for ele in data])
    y_norm = y / sum(y)
    y_norm_err = np.sqrt(y) / sum(y)
    y_mixed = np.asarray([ele[particles] for ele in data_mixed])
    y_mixed_norm = y_mixed / sum(y_mixed)
    y_mixed_norm_err = np.sqrt(y_mixed) / sum(y_mixed)
    x = range(len(y))
    y_binom = binom.pmf(x, particles, 1/divs)

    fig1, ax1 = plt.subplots()
    ax1.errorbar(x, y_norm, yerr=y_norm_err, fmt='ob', zorder=2, label=f'{particles} Proton Events')
    ax1.errorbar(x, y_mixed_norm, yerr=y_mixed_norm_err, fmt='og', zorder=1, label=f'{particles} Proton Mixed Events')
    ax1.scatter(x, y_binom, color='red', zorder=0, label='Binomial Distribution')
    ax1.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    ax1.set_xticks(range(0, len(y), 2))
    ax1.set_title(f'Particles in {divs} Division Bin vs Binomial for {particles} Proton Events'+title_sufx)
    ax1.set_xlabel('Number of Particles in Bin')
    ax1.set_ylabel('Events')
    ax1.legend()

    y_err = np.sqrt(y)
    y_mixed_err = np.sqrt(y_mixed)
    x = range(len(y))
    y_binom2 = binom.pmf(x, particles, 1 / divs) * sum(y)

    fig9, ax9 = plt.subplots()
    ax9.errorbar(x, y, yerr=y_err, fmt='ob', zorder=2, label=f'{particles} Proton Events')
    # ax9.errorbar(x, y_mixed_norm, yerr=y_mixed_norm_err, fmt='og', zorder=1, label=f'{particles} Proton Mixed Events')
    ax9.scatter(x, y_binom2, color='red', zorder=0, label='Binomial Distribution')
    ax9.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    ax9.set_xticks(range(0, len(y), 2))
    ax9.set_title(f'Particles in {divs} Division Bin vs Binomial for {particles} Proton Events'+title_sufx)
    ax9.set_xlabel('Number of Particles in Bin')
    ax9.set_ylabel('Events')
    ax9.legend()
    # fig2, ax2 = plt.subplots()
    # y_diff = y_norm - y_mixed_norm
    # diff_err = np.sqrt(y_norm_err**2 + y_mixed_norm_err**2)
    # ax2.axhline(0, color='red', ls='--')
    # ax2.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    # ax2.errorbar(x, y_diff, yerr=diff_err, fmt='bo')
    # ax2.set_xticks(range(0, len(y), 2))
    # ax2.set_title(f'Particles in {divs} Division Bin Minus Mixed for {particles} Proton Events')
    # ax2.set_xlabel('Number of Particles in Bin')
    # ax2.set_ylabel('Data Events Minus Mixed')
    #
    # fig3, ax3 = plt.subplots()
    # x_div = []
    # y_div = []
    # y_bin_div = []
    # bin_div_err = []
    # div_err = []
    # for i in range(len(y_norm)):
    #     if y_mixed_norm[i] != 0 and y_norm[i] != 0:
    #         x_div.append(x[i])
    #         y_div.append(y_norm[i] / y_mixed_norm[i])
    #         y_bin_div.append(y_norm[i] / y_binom[i])
    #         div_err.append(abs(y_div[-1]) * np.sqrt((y_norm_err[i] / y_norm[i])**2 +
    #                                                 (y_mixed_norm_err[i] / y_mixed_norm[i])**2))
    #         bin_div_err.append(y_norm_err[i] / y_binom[i])
    # ax3.axhline(1, color='red', ls='--')
    # ax3.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    # ax3.errorbar(x_div, y_div, yerr=div_err, zorder=1, fmt='bo', label='Data Divided by Mixed')
    # ax3.errorbar(x_div, y_bin_div, yerr=bin_div_err, zorder=0, fmt='ro', label='Data Divided by Binomial')
    # ax3.set_xticks(range(0, max(x_div), 2))
    # ax3.set_title(f'Particles in {divs} Division Bin Divided by Mixed for {particles} Proton Events')
    # ax3.set_xlabel('Number of Particles in Bin')
    # ax3.set_ylabel('Data Events Divided by Mixed')
    # ax3.legend()
    #
    # fig4, ax4 = plt.subplots()
    # x_data_div = []
    # x_mix_div = []
    # y_data_bin_div = []
    # y_mix_bin_div = []
    # data_bin_div_err = []
    # mix_bin_div_err = []
    # for i in range(len(y_norm)):
    #     if y_norm[i] != 0:
    #         x_data_div.append(x[i])
    #         y_data_bin_div.append(y_norm[i] / y_binom[i])
    #         data_bin_div_err.append(y_norm_err[i] / y_binom[i])
    #     if y_mixed_norm[i] != 0:
    #         x_mix_div.append(x[i])
    #         y_mix_bin_div.append(y_mixed_norm[i] / y_binom[i])
    #         mix_bin_div_err.append(y_mixed_norm_err[i] / y_binom[i])
    # ax4.axhline(1, color='red', ls='--')
    # ax4.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    # ax4.errorbar(x_data_div, y_data_bin_div, yerr=data_bin_div_err, zorder=1, fmt='bo', label='Data Divided by Binomial')
    # ax4.errorbar(x_mix_div, y_mix_bin_div, yerr=mix_bin_div_err, zorder=0, fmt='go', label='Mixed Divided by Binomial')
    # ax4.set_xticks(range(0, max([max(x_data_div),max(x_mix_div)]), 2))
    # ax4.set_title(f'Particles in {divs} Division Bin Divided by Binomial for {particles} Proton Events')
    # ax4.set_xlabel('Number of Particles in Bin')
    # ax4.set_ylabel('Data/Mixed Events Divided by Binomial')
    # ax4.legend()
    #
    # fig5, ax5 = plt.subplots()
    # y_diff_raw = y_norm - y_binom
    # y_diff_mix = y_mixed_norm - y_binom
    # ax5.axhline(0, color='red', ls='--')
    # ax5.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    # ax5.errorbar(x, y_diff_raw, yerr=y_norm_err, fmt='bo', label='Raw - Binomial')
    # ax5.errorbar(x, y_diff_mix, yerr=y_mixed_norm_err, fmt='go', label='Mix - Binomial')
    # ax5.set_xticks(range(0, len(y), 2))
    # ax5.set_title(f'Particles in {divs} Division Bin Minus Binomial for {particles} Proton Events')
    # ax5.set_xlabel('Number of Particles in Bin')
    # ax5.set_ylabel('Data Events Minus Mixed')
    # ax5.legend()
    #
    # fig6, axes1 = plt.subplots(4, 4, sharex='all', sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    # fig7, axes2 = plt.subplots(4, 4, sharex='all', sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    # fig7.suptitle('Data (blue) and Mix (green) Minus Binomial for 10-26 Total Particles')
    # for index, particle in enumerate(particle_array):
    #     yi = np.asarray([ele[particle] for ele in data])
    #     y_normi = yi / sum(yi)
    #     y_norm_erri = np.sqrt(yi) / sum(yi)
    #     y_mixedi = np.asarray([ele[particle] for ele in data_mixed])
    #     y_mixed_normi = y_mixedi / sum(y_mixedi)
    #     y_mixed_norm_erri = np.sqrt(y_mixedi) / sum(y_mixedi)
    #     xi = range(len(yi))
    #     y_binomi = binom.pmf(xi, particle, 1 / divs)
    #     axes1[int(index / 4), int(index % 4)].errorbar(xi, y_normi, yerr=y_norm_erri, fmt='ob', zorder=2,
    #                                                    label=f'{particle} Proton Events')
    #     axes1[int(index / 4), int(index % 4)].errorbar(xi, y_mixed_normi, yerr=y_mixed_norm_erri, fmt='og', zorder=1,
    #                                                    label=f'{particle} Proton Mixed Events')
    #     axes1[int(index / 4), int(index % 4)].scatter(xi, y_binomi, color='red', zorder=0,
    #                                                   label='Binomial Distribution')
    #     axes1[int(index / 4), int(index % 4)].axvline(float(particle) / divs, color='red', ls='--', label='Mean')
    #     axes1[int(index / 4), int(index % 4)].set_xticks(range(0, len(yi), 4))
    #     axes1[int(index / 4), int(index % 4)].set_xlim(0, 22)
    #
    #     y_diffi = y_normi - y_binomi
    #     y_diff_mixi = y_mixed_normi - y_binomi
    #     axes2[int(index / 4), int(index % 4)].errorbar(xi, y_diffi, yerr=y_norm_erri, fmt='ob', zorder=2,
    #                                                    label=f'{particle} Proton Events')
    #     axes2[int(index / 4), int(index % 4)].errorbar(xi, y_diff_mixi, yerr=y_mixed_norm_erri, fmt='og', zorder=1,
    #                                                    label=f'{particle} Proton Mixed Events')
    #     axes2[int(index / 4), int(index % 4)].axvline(float(particle) / divs, color='red', ls='--', label='Mean')
    #     axes2[int(index / 4), int(index % 4)].axhline(0, color='red', ls='--')
    #     axes2[int(index / 4), int(index % 4)].set_xticks(range(0, len(yi), 4))
    #     axes2[int(index / 4), int(index % 4)].set_xlim(-1, 17)
    #     # axes[int(index / 4), int(index % 4)].set_title(f'{particle} Total Particles')
    #     # axes[int(index / 4), int(index % 4)].set_xlabel('Number of Particles in Bin')
    #     # axes[int(index / 4), int(index % 4)].set_ylabel('Events')
    #     # axes[int(index / 4), int(index % 4)].legend()
    # axes1[0, 0].legend()
    # axes2[0, 0].legend()

    # fig8, ax8 = plt.subplots()
    # x_raw = []
    # raw_sd = []
    # x_mix = []
    # mix_sd = []
    # for total_particles in enumerate(data):
    #     x_raw.append(total_particles)
    #     raw_sd.append(hist_sd(range(len(data)), y))
    # mix_sd = hist_sd(x, y_mixed)
    # ax5.axhline(0, color='red', ls='--')
    # ax5.axvline(float(particles) / divs, color='red', ls='--', label='Mean')
    # ax5.errorbar(x, y_diff_raw, yerr=y_norm_err, fmt='bo', label='Raw - Binomial')
    # ax5.errorbar(x, y_diff_mix, yerr=y_mixed_norm_err, fmt='go', label='Mix - Binomial')
    # ax5.set_xticks(range(0, len(y), 2))
    # ax5.set_title(f'Particles in {divs} Division Bin Minus Binomial for {particles} Proton Events')
    # ax5.set_xlabel('Number of Particles in Bin')
    # ax5.set_ylabel('Data Events Minus Mixed')
    # ax5.legend()

    plt.show()


def plot_azbin_data(data, x_lim, y_lim, divs, x_label='Number of Particles in Bin', y_label='Number of Particles in Event',
                    title_sufx=''):
    # x, y = np.meshgrid(x_range, y_range)
    # x, y = np.meshgrid(np.asarray(x_range) - float(x_range[1] - x_range[0]) / 2,
    #                    np.asarray(y_range) - float(y_range[1] - y_range[0]) / 2)
    x_range = range(0, len(data))
    y_range = range(0, len(data[0]))
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
    plt.ylim(y_lim)
    plt.xlim(x_lim)
    plt.legend()
    plt.show()


def plot_azbin_data_trans(data, x_lim, y_lim, divs, x_label='Number of Particles in Bin',
                               y_label='Number of Particles in Event', title_sufx=''):
    ax_2d = plt.axes()

    # x, y = np.meshgrid(x_range, y_range)
    x_range = range(0, len(data))
    y_range = range(0, len(data[0]))
    x, y = np.meshgrid(np.asarray(x_range) - float(x_range[1] - x_range[0]) / 2,
                       np.asarray(y_range) - float(y_range[1] - y_range[0]) / 2)
    data = np.ma.masked_where(data <= 0, data)
    plot_2d = ax_2d.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    y2 = 360 / float(divs) * np.asarray(y_range)
    y3 = 1.0 * np.asarray(y_range)
    ax_2d.plot(y_range, y2, color='red', label='mean')
    ax_2d.plot(y_range, y3, color='blue', label='max')
    ax_2d.set_xlabel(x_label)
    ax_2d.set_ylabel(y_label)
    ax_2d.set_ylim(y_lim)
    ax_2d.set_xlim(x_lim)
    ax_2d.set_xticks(range(x_lim[0], x_lim[1], 2))
    ax_2d.set_title('Particles in Event vs Particles in Bin'+title_sufx)
    ax_2d.legend(loc='lower right')
    plt.colorbar(plot_2d)

    plt.show()


def plot_azbin_data_trans_proj(data, x_lim, y_lim, divs, x_label='Number of Particles in Bin',
                          y_label='Number of Particles in Event', title_sufx=''):
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
    plot_2d = ax_2d.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
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

    x_proj = np.sum(data, axis=1)
    y_proj = np.sum(data, axis=0)

    ax_x.bar(x_range, x_proj, width=1, align='center')
    ax_y.barh(y_range, y_proj, height=1, align='center')

    cbar_ax = plt.axes([left+width+0.03, bottom_h-0.035, 0.025, 0.29])
    plt.colorbar(plot_2d, cax=cbar_ax)

    # title_ax = plt.axes([left_h, bottom_h, 0.15, 0.2])
    plt.text(0.78, 0.85, 'Particles in Event vs Particles in Bin'+title_sufx, fontsize=12, transform=plt.gcf().transFigure,
             wrap=True, va='center', fontweight='bold')

    plt.show()


def plot_ratio_data(data, x_range, y_range, divs, x_label='Number of Particles in Bin',
                    y_label='Number of Particles in Event', title_sufx=''):
    # x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
    #                    np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    x, y = np.meshgrid(x_range, y_range)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(float(divs) / 360, color='red', label='mean')
    plt.axvline(1.05, color='blue', label='max')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title('Particles in Event vs Ratios'+title_sufx)
    plt.legend(loc='upper right')
    plt.show()


def plot_ratio_dist(y, x, divs, title_sufx=''):
    plt.bar(x, y, width=x[1]-x[0], align='edge', zorder=0)
    plt.yscale('log')
    plt.axvline(float(divs) / 360, color='red', label='mean')
    plt.axvline(1.05, color='blue', label='max')
    plt.xlabel('Ratio')
    plt.ylabel('Events')
    plt.title('Ratio Distribution'+title_sufx)
    plt.legend(loc='upper right')
    plt.show()


def plot_ratio_kde(data):
    sns.set()
    ax = sns.distplot(data, rug=True, hist=False)
    plt.show()


def plot_diff_data(data, x_range, y_range, divs, max_particles=40, x_label='Number of Particles in Bin',
                   y_label='Number of Particles in Event', title_sufx=''):
    # x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
    #                    np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    x, y = np.meshgrid(x_range, y_range)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(0.0, color='red', label='mean')
    print(divs)
    y2 = np.asarray(y_range) * (1 - divs / 360)
    y3 = -float(divs) / 360 * np.asarray(y_range)
    plt.plot(y2, y_range, color='blue', label='max')
    plt.plot(y3, y_range, color='black', label='min')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, max_particles])
    plt.xlim([-10, 12])
    plt.xticks(range(-10, 12, 2))
    plt.title('Particles in Event vs Differences'+title_sufx)
    plt.legend(loc='lower left')
    plt.show()


def plot_pull_data(data, x_range, y_range, divs, max_particles=40, x_label='Number of Particles in Bin',
                   y_label='Number of Particles in Event', title_sufx=''):
    # x, y = np.meshgrid(np.asarray(x_range)-float(x_range[1]-x_range[0])/2,
    #                    np.asarray(y_range)-float(y_range[1]-y_range[0])/2)
    x, y = np.meshgrid(x_range, y_range)
    data = np.ma.masked_where(data <= 0, data)
    plt.pcolormesh(y, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
    plt.axvline(0.0, color='red', label='mean')
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([0, max_particles])
    plt.xlim([-10, 12])
    plt.xticks(range(-10, 12, 2))
    plt.title('Particles in Event vs Pull'+title_sufx)
    plt.legend(loc='lower left')
    plt.show()


def plot_diff_dist(y, x, divs, title_sufx=''):
    plt.bar(x, y, width=x[1]-x[0], align='edge', zorder=0)
    plt.yscale('log')
    plt.axvline(0.0, color='red', label='mean')
    plt.xlabel('Difference')
    plt.ylabel('Events')
    plt.title('Difference Distribution'+title_sufx)
    plt.xlim([-10, 12])
    plt.xticks(range(-10, 12, 2))
    plt.legend(loc='upper right')
    plt.show()


def plot_diff_kde(data):
    sns.set()
    ax = sns.distplot(data, rug=True, hist=False)
    plt.show()


def ratio_transform(data, divs, max_particles=40, title_sufx=''):
    ratio_data = np.zeros((21, max_particles))
    ratio_dist = np.zeros(21)
    ratio_vals = []
    for num_in_bin, bin_data in enumerate(data):
        for total_particles, num_events in enumerate(bin_data):
            if total_particles > 0 and num_events > 0:
                ratio = float(num_in_bin) / total_particles
                ratio_bin = int(ratio * 20)
                ratio_data[ratio_bin][total_particles] += num_events
                ratio_dist[ratio_bin] += num_events
                # ratio_vals.extend([ratio] * int(num_events))

    plot_ratio_data(ratio_data, range(0, max_particles), np.linspace(0, 1.05, 22), divs, x_label='Ratio', title_sufx=title_sufx)
    plot_ratio_dist(ratio_dist, np.linspace(0, 1, 21), divs, title_sufx=title_sufx)
    # plot_ratio_kde(ratio_vals)


def diff_transform(data, divs, max_particles=40, title_sufx=''):
    x_bins = 32
    diff_data = np.zeros((x_bins, max_particles))
    diff_dist = np.zeros(x_bins)
    diff_values = []
    bin_width = divs
    divs = int(360 / divs)  # Convert from divs in bin width to number of divisions
    for num_in_bin, bin_data in enumerate(data):
        for total_particles, num_events in enumerate(bin_data):
            if total_particles >= num_in_bin and total_particles > 0 and num_events > 0:
                # diff = float(num_in_bin) - float(total_particles) / divs
                diff = Fraction(num_in_bin * divs - total_particles, divs)
                if diff < Fraction(-total_particles, divs):
                    print('less than min')
                if diff > Fraction(total_particles * (divs - 1), divs):
                    print(f'greater than max: {total_particles * (divs - 1) / divs} | nT: {total_particles} | '
                          f'np: {num_in_bin} | diff: {diff}')
                diff_norm = Fraction((diff + Fraction(max_particles, divs)), max_particles)
                diff_bin = int(diff_norm * (x_bins - 1))
                diff_data[diff_bin][total_particles] += num_events
                diff_dist[diff_bin] += num_events
                diff_values.extend([diff] * int(num_events))


    plot_diff_data(diff_data, range(0, max_particles), np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins),
                   bin_width, x_label='Difference', title_sufx=title_sufx, max_particles=max_particles)
    plot_diff_dist(diff_dist, np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins), bin_width, title_sufx=title_sufx)
    # plot_diff_kde(diff_values)


def pull_transform(data, data_mix, divs, max_particles=40, title_sufx=''):
    x_bins = 32
    diff_mix_data = np.zeros((x_bins, max_particles))
    diff_mix_dist = np.zeros(x_bins)
    diff_mix_values = []

    bin_width = divs
    divs = int(360 / divs)  # Convert from divs in bin width to number of divisions
    for num_in_bin, bin_data in enumerate(data_mix):
        for total_particles, num_events in enumerate(bin_data):
            diff_mix_values.append([])
            if total_particles >= num_in_bin and total_particles > 0 and num_events > 0:
                # diff = float(num_in_bin) - float(total_particles) / divs
                diff = Fraction(num_in_bin * divs - total_particles, divs)
                if diff < Fraction(-total_particles, divs):
                    print('less than min')
                if diff > Fraction(total_particles * (divs - 1), divs):
                    print(f'greater than max: {total_particles * (divs - 1) / divs} | nT: {total_particles} | '
                          f'np: {num_in_bin} | diff: {diff}')
                diff_norm = Fraction((diff + Fraction(max_particles, divs)), max_particles)
                diff_bin = int(diff_norm * (x_bins - 1))
                diff_mix_data[diff_bin][total_particles] += num_events
                diff_mix_dist[diff_bin] += num_events
                diff_mix_values[total_particles].extend([diff] * int(num_events))

    # diff_mix_sds = np.zeros(max_particles)
    # for total_particles in range(len(diff_mix_sds)):
    #     if len(diff_mix_values[total_particles]) > 2:
    #         diff_mix_sds[total_particles] = fraction_sd(diff_mix_values[total_particles])
    #     else:
    #         diff_mix_sds[total_particles] = 1
    #     print(f'{total_particles} total particles, mix diff slice sd: {diff_mix_sds[total_particles]}')
    diff_binom_sds = np.zeros(max_particles)
    for total_particles in range(len(diff_binom_sds)):
        p = float(bin_width) / 360
        diff_binom_sds[total_particles] = (total_particles * p * (1 - p))**0.5

    diff_data = np.zeros((x_bins, max_particles))
    diff_dist = np.zeros(x_bins)
    pull_data = np.zeros((x_bins, max_particles))
    pull_dist = np.zeros(x_bins)
    diff_values = []
    pull_values = []
    for num_in_bin, bin_data in enumerate(data):
        for total_particles, num_events in enumerate(bin_data):
            pull_values.append([])
            if total_particles >= num_in_bin and total_particles > 0 and num_events > 0:
                # diff = float(num_in_bin) - float(total_particles) / divs
                diff = Fraction(num_in_bin * divs - total_particles, divs)
                if diff < Fraction(-total_particles, divs):
                    print('less than min')
                if diff > Fraction(total_particles * (divs - 1), divs):
                    print(f'greater than max: {total_particles * (divs - 1) / divs} | nT: {total_particles} | '
                          f'np: {num_in_bin} | diff: {diff}')
                diff_norm = Fraction((diff + Fraction(max_particles, divs)), max_particles)
                diff_bin = int(diff_norm * (x_bins - 1))
                diff_data[diff_bin][total_particles] += num_events
                diff_dist[diff_bin] += num_events
                pull = diff / diff_binom_sds[total_particles]
                pull_norm = (pull + 20 / divs) / 20
                pull_bin = int(pull_norm * (x_bins - 1))
                if pull_bin == x_bins:
                    print('Over max')
                    pull_bin -= 1
                pull_data[pull_bin][total_particles] += num_events
                pull_dist[pull_bin] += num_events
                diff_values.extend([diff] * int(num_events))
                pull_values[total_particles].extend([pull] * int(num_events))

    # for total_particles in range(len(pull_values)):
    #     if len(pull_values[total_particles]) > 2:
    #         sd = fraction_sd(diff_mix_values[total_particles])
    #     else:
    #         sd = 1
    #     print(f'{total_particles} total particles, pull slice sd: {sd}')

    plot_diff_data(diff_mix_data, range(0, max_particles), np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins),
                   bin_width, x_label='Difference', title_sufx=title_sufx + ' Mix', max_particles=max_particles)

    plot_diff_data(diff_data, range(0, max_particles), np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins),
                   bin_width, x_label='Difference', title_sufx=title_sufx, max_particles=max_particles)
    plot_diff_dist(diff_dist, np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins), bin_width, title_sufx=title_sufx)

    plot_pull_data(pull_data, range(0, max_particles), np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins),
                   bin_width, x_label='Pull', title_sufx=title_sufx, max_particles=max_particles)
    plot_diff_dist(pull_dist, np.linspace(-float(max_particles) / divs, float(max_particles) - float(max_particles) / divs, x_bins), bin_width, title_sufx=title_sufx)
    # plot_diff_kde(diff_values)


def fraction_sd(x):
    mean = 0
    for xi in x:
        mean += xi
    mean /= len(x)

    variance = 0
    for xi in x:
        variance += (xi - mean)**2
    variance = Fraction(variance, len(x) - 1)

    return variance ** Fraction(1, 2)


def hist_sd(x, y):
    mean = hist_mean(x, y)
    variance = sum((np.asarray(x) - mean) * np.asarray(y)) / (sum(y) - 1)  # Think I'm missing **2?
    return variance**0.5


def hist_mean(x, y):
    mean = sum(np.asarray(x) * np.asarray(y)) / sum(y)
    return mean


if __name__ == '__main__':
    main()
