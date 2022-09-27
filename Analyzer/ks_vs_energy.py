#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 06 8:34 PM 2020
Created in PyCharm
Created as QGP_Scripts/ks_vs_energy.py

@author: Dylan Neff, dylan
"""

from AzimuthBinData import AzimuthBinData as AzData
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Measure import Measure


def main():
    # from_hist_files()
    from_dataframe()


def from_dataframe():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    divs = 120
    energy = 39
    data_set_name = 'bes_resample_def'
    stat = 'standard deviation'

    df_path = base_path + df_name
    df = pd.read_csv(df_path)
    df = df.dropna()
    print(df.head())
    df = df[(df['name'] == data_set_name) & (df['divs'] == divs) & (df['energy'] == energy) & (df['stat'] == stat)]
    df.sort_values(by=['total_protons'])
    df_raw = df[df['data_type'] == 'raw']
    df_mix = df[df['data_type'] == 'mix']

    p = float(divs) / 360
    y_binom_mix = (np.asarray(df_mix['total_protons']) * p * (1 - p)) ** 0.5
    y_binom_raw = (np.asarray(df_raw['total_protons']) * p * (1 - p)) ** 0.5

    fig1, ax1 = plt.subplots()
    ax1.errorbar(df_raw['total_protons'], df_raw['val'], df_raw['err'], alpha=0.8, zorder=2, color='blue',
                 ls='', marker='o', label='Raw')
    ax1.errorbar(df_mix['total_protons'], df_mix['val'], df_mix['err'], alpha=0.8, zorder=1, color='green',
                 ls='', marker='o', label='Mix')
    ax1.plot(df_mix['total_protons'], y_binom_mix, color='red', alpha=0.8, zorder=0, label='Binomial')
    ax1.set_xlabel('Total Protons in Event')
    ax1.set_ylabel('Standard Deviation of Slice')
    ax1.set_title(f'{energy}GeV, 0-5% Centrality, {divs}° Partitions')
    ax1.legend()
    fig1.tight_layout()

    fig2, ax2 = plt.subplots()
    raw_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
    mix_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
    ax2.errorbar(df_raw['total_protons'], [x.val for x in raw_ratio], [x.err for x in raw_ratio], ls='', marker='o',
                 zorder=2, color='blue', alpha=0.8, label='Raw / Binomial')
    ax2.errorbar(df_mix['total_protons'], [x.val for x in mix_ratio], [x.err for x in mix_ratio], ls='', marker='o',
                 zorder=1, color='green', alpha=0.8, label='Mix / Binomial')
    ax2.axhline(1, zorder=0, color='red', ls='--')
    ax2.set_xlabel('Total Particles')
    ax2.set_ylabel('Standard Deviation Ratio')
    ax2.set_title(f'{energy}GeV, 0-5% Centrality, {divs}° Partitions')
    ax2.legend()
    fig2.tight_layout()

    plt.show()


def from_hist_files():
    energies = [7, 11, 19, 27, 39, 62]
    divisions = [120]
    centralities = [8]
    # path = {'raw': ['F:/Research/Data_Ampt_Old/default/Ampt_rapid05_n1ratios_0/', 'GeV/ratios_divisions_', '_centrality_', '_local.txt'],
    #         'mix': ['F:/Research/Data_Ampt_Old_Mix/default/Ampt_rapid05_n1ratios_0/', 'GeV/ratios_divisions_', '_centrality_',
    #                 '_local.txt']}
    path = {
        'raw': ['F:/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/',
                'GeV/ratios_divisions_', '_centrality_', '_local.txt'],
        'mix': ['F:/Research/Data_Mix/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/',
                'GeV/ratios_divisions_', '_centrality_', '_local.txt']}
    colors = {7: 'red', 11: 'blue', 19: 'green', 27: 'cyan', 39: 'magenta', 62: 'black'}
    cdf_plot = {'energy': 39, 'division': 120, 'centrality': 8, 'total particles': 20}
    sd_plot = {'energy': 39, 'division': 120, 'centrality': 8}
    title_sufx = f'\n{sd_plot["energy"]}GeV, 0-5% Centrality, {sd_plot["division"]}° Bins'

    data = {'raw': {}, 'mix': {}}
    for source in data.keys():
        for energy in energies:
            data[source][energy] = {}
            for division in divisions:
                data[source][energy][division] = {}
                for centrality in centralities:
                    current_path = path[source][0] + str(energy)+path[source][1] + str(division) + path[source][2] +\
                                   str(centrality) + path[source][3]
                    data[source][energy][division][centrality] = AzData(path=current_path, div=division)

    plot_data_ks2 = {}
    plot_data_ks2_test = {}
    plot_data_raw = {}
    plot_data_mix = {}
    sd_data_raw = {}
    sd_data_mix = {}
    # plot_data_cdf = {}
    for energy in energies:
        print(f'Working on {energy}GeV')
        plot_data_ks2[energy] = [[], []]
        plot_data_ks2_test[energy] = [[], []]
        plot_data_raw[energy] = [[], []]
        plot_data_mix[energy] = [[], []]
        sd_data_raw[energy] = [[], []]
        sd_data_mix[energy] = [[], []]
        for division in divisions:
            for centrality in centralities:
                raw = data['raw'][energy][division][centrality]
                mix = data['mix'][energy][division][centrality]
                for total_particles in raw.data:
                    if total_particles in mix.data:
                        plot_data_ks2_test[energy][0].append(total_particles)
                        plot_data_ks2_test[energy][1].append(ks2_test(raw.data[total_particles], mix.data[total_particles]))
                        plot_data_raw[energy][0].append(total_particles)
                        plot_data_raw[energy][1].append(ks2_test(raw.data[total_particles], stats.binom.pmf(
                            range(0, total_particles+1), total_particles, 1.0/raw.div)))
                        plot_data_mix[energy][0].append(total_particles)
                        plot_data_mix[energy][1].append(ks2_test(mix.data[total_particles], stats.binom.pmf(
                            range(0, total_particles+1), total_particles, 1.0/mix.div)))
                        sd_data_raw[energy][0].append(total_particles)
                        sd_data_raw[energy][1].append(hist_sd(range(len(raw.data[total_particles])), raw.data[total_particles]))
                        sd_data_mix[energy][0].append(total_particles)
                        sd_data_mix[energy][1].append(
                            hist_sd(range(len(mix.data[total_particles])), mix.data[total_particles]))
                        if sd_data_raw[energy][1][-1] is None or sd_data_mix[energy][1][-1] is None:
                            sd_data_raw[energy][0].pop()
                            sd_data_raw[energy][1].pop()
                            sd_data_mix[energy][0].pop()
                            sd_data_mix[energy][1].pop()
                        # if energy == cdf_plot['energy'] and total_particles == cdf_plot['total particles'] and \
                        #         division == cdf_plot['division'] and centrality == cdf_plot['centrality']:
                        #     plot_data_cdf['raw']
                        # raw_dist = [bin_particles for bin_particles, events in enumerate(raw.data[total_particles])
                        #             for l in range(events)]
                        # mix_dist = [bin_particles for bin_particles, events in enumerate(mix.data[total_particles])
                        #             for l in range(events)]
                        # ks2 = stats.ks_2samp(raw_dist, mix_dist)
                        # ks_raw = stats.kstest(raw_dist, 'binom',
                        #                       args=(total_particles, 1.0/raw.div))
                        # ks_mix = stats.kstest(mix_dist, 'binom',
                        #                       args=(total_particles, 1.0/mix.div))
                        # # print(total_particles)
                        # # plt.hist(mix_dist, bins=np.linspace(-0.5, total_particles+0.5, total_particles+2), align='mid', density=True, zorder=0)
                        # # plt.scatter(range(0, total_particles+1), stats.binom.pmf(range(0, total_particles+1), total_particles, 1.0/mix.div),
                        # #          zorder=1, color='red')
                        # # plt.show()
                        # plot_data_ks2[energy][0].append(total_particles)
                        # plot_data_ks2[energy][1].append(ks2[0])
                        # plot_data_raw[energy][0].append(total_particles)
                        # plot_data_raw[energy][1].append(ks_raw[0])
                        # plot_data_mix[energy][0].append(total_particles)
                        # plot_data_mix[energy][1].append(ks_mix[0])

    raw_cdf = get_cdf(data['raw'][cdf_plot['energy']][cdf_plot['division']][cdf_plot['centrality']]
                      .data[cdf_plot['total particles']])
    mix_cdf = get_cdf(data['mix'][cdf_plot['energy']][cdf_plot['division']][cdf_plot['centrality']]
                      .data[cdf_plot['total particles']])

    # fig1, ax1 = plt.subplots()
    # fig2, ax2 = plt.subplots()
    # fig3, ax3 = plt.subplots()
    # fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()

    # for energy in energies:
    #     # ax1.plot(plot_data_ks2[energy][0], plot_data_ks2[energy][1], color=colors[energy], label=f'{energy}GeV')
    #     ax2.plot(plot_data_raw[energy][0], plot_data_raw[energy][1], color=colors[energy], label=f'{energy}GeV')
    #     ax3.plot(plot_data_mix[energy][0], plot_data_mix[energy][1], color=colors[energy], label=f'{energy}GeV')
    #     ax4.plot(plot_data_ks2_test[energy][0], plot_data_ks2_test[energy][1], color=colors[energy],
    #              label=f'{energy}GeV')
    #
    # ax1.scatter(range(0, len(raw_cdf)), raw_cdf, label='Raw CDF')
    # ax1.scatter(range(0, len(mix_cdf)), mix_cdf, label='Mix CDF')
    # ax1.legend()
    # ax1.set_title(f'{cdf_plot["energy"]}GeV, {cdf_plot["division"]} div, {cdf_plot["centrality"]} cent,'
    #               f' {cdf_plot["total particles"]} particles, CDF Comparison')
    # ax1.set_xlabel('Particles in Bin')
    # ax1.set_ylabel('Integrated Probability')
    #
    # ax2.legend()
    # ax2.set_title('KS Statistics for Data to Binomial | 3 Divisions Most Central')
    # ax2.set_xlabel('Total Particles in Event')
    # ax2.set_ylabel('KS Statistic')
    #
    # ax3.legend()
    # ax3.set_title('KS Statistics for Mixed to Binomial | 3 Divisions Most Central')
    # ax3.set_xlabel('Total Particles in Event')
    # ax3.set_ylabel('KS Statistic')
    #
    # ax4.legend()
    # ax4.set_title('Test Two Sample KS Statistics for Data to Mixed | 3 Divisions Most Central')
    # ax4.set_xlabel('Total Particles in Event')
    # ax4.set_ylabel('KS Statistic')

    raw_x = sd_data_raw[sd_plot['energy']][0]
    raw_y = sd_data_raw[sd_plot['energy']][1]
    mix_x = sd_data_mix[sd_plot['energy']][0]
    mix_y = sd_data_mix[sd_plot['energy']][1]
    p = float(sd_plot['division']) / 360
    y_mix = (np.asarray(mix_x) * p * (1 - p)) ** 0.5

    raw_x, raw_y = list(zip(*[(x, y) for x, y in zip(raw_x, raw_y) if y > 0]))
    y_raw = (np.asarray(raw_x) * p * (1 - p)) ** 0.5

    # fig5.set_size_inches(7, 7)
    ax5.scatter(raw_x, raw_y, zorder=2, color='blue', label='Raw SD')
    ax5.scatter(mix_x, mix_y, zorder=1, color='green', label='Mix SD')
    ax5.plot(mix_x, y_mix, color='red', zorder=0, label='Binomial SD')
    ax5.set_xlabel('Total Particles')
    ax5.set_ylabel('Standard Deviation of Slice')
    ax5.set_title(f'Standard Deviation of Total Particle Slices for {sd_plot["energy"]}GeV'+title_sufx)
    ax5.legend()

    # fig6.set_size_inches(10, 7)
    raw_ratio = raw_y / y_raw
    mix_ratio = mix_y / y_mix
    ax6.scatter(raw_x, raw_ratio, zorder=2, color='blue', label='Raw SD / Binomial SD')
    ax6.scatter(mix_x, mix_ratio, zorder=1, color='green', label='Mix SD / Binomial SD')
    ax6.axhline(1, zorder=0, color='red', ls='--')
    # ax6.axhline(np.average(raw_ratio), zorder=0, color='blue', ls='--', label='Raw Avg')
    # ax6.axhline(np.average(mix_ratio), zorder=0, color='green', ls='--', label='Mix Avg')
    ax6.set_xlabel('Total Particles')
    ax6.set_ylabel('Standard Deviation of Slice Divided by Binomial')
    ax6.set_title(f'SD Divided by Binomial of Total Particle Slices for {sd_plot["energy"]}GeV'+title_sufx)
    ax6.legend()

    plt.show()

    print('donzo')


def ks2_test(dist1, dist2):
    ks = -1
    cdf1 = get_cdf(dist1)
    cdf2 = get_cdf(dist2)

    for i, j in zip(cdf1, cdf2):
        if abs(i-j) > ks:
            ks = abs(i-j)

    return ks


def get_cdf(y):
    cdf = [0]
    norm = float(sum(y))
    for i in y:
        cdf.append(i / norm + cdf[-1])

    return cdf[1:]


def hist_sd(x, y):
    mean = hist_mean(x, y)
    y_sum = sum(y)
    if y_sum == 1:
        print('Only one event what is this amateur hour')
        return None
    else:
        variance = sum((np.asarray(x) - mean)**2 * np.asarray(y)) / (sum(y) - 1)
        return variance**0.5


def hist_mean(x, y):
    mean = sum(np.asarray(x) * np.asarray(y)) / sum(y)
    return mean


if __name__ == '__main__':
    main()
