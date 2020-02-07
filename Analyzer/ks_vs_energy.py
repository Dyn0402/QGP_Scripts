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


def main():
    energies = [7, 11, 19, 27, 39, 62]
    divisions = [3]
    centralities = [8]
    path = {'raw': ['/home/dylan/Research/Data/Single_Ratio0/', 'GeV/ratios_divisions_', '_centrality_', '_local.txt'],
            'mix': ['/home/dylan/Research/Data_Mix/Single_Ratio0/', 'GeV/ratios_divisions_', '_centrality_',
                    '_local.txt']}
    colors = {7: 'red', 11: 'blue', 19: 'green', 27: 'cyan', 39: 'magenta', 62: 'black'}
    cdf_plot = {'energy': 7, 'division': 3, 'centrality': 8, 'total protons': 35}

    data = {'raw': {}, 'mix': {}}
    for source in data.keys():
        for energy in energies:
            data[source][energy] = {}
            for division in divisions:
                data[source][energy][division] = {}
                for centrality in centralities:
                    current_path = path[source][0]+str(energy)+path[source][1]+str(division)+path[source][2]+\
                                   str(centrality)+path[source][3]
                    data[source][energy][division][centrality] = AzData(path=current_path, div=division)

    plot_data_ks2 = {}
    plot_data_ks2_test = {}
    plot_data_raw = {}
    plot_data_mix = {}
    # plot_data_cdf = {}
    for energy in energies:
        print(f'Working on {energy}GeV')
        plot_data_ks2[energy] = [[], []]
        plot_data_ks2_test[energy] = [[], []]
        plot_data_raw[energy] = [[], []]
        plot_data_mix[energy] = [[], []]
        for division in divisions:
            for centrality in centralities:
                raw = data['raw'][energy][division][centrality]
                mix = data['mix'][energy][division][centrality]
                for total_protons in raw.data:
                    if total_protons in mix.data:
                        plot_data_ks2_test[energy][0].append(total_protons)
                        plot_data_ks2_test[energy][1].append(ks2_test(raw.data[total_protons], mix.data[total_protons]))
                        plot_data_raw[energy][0].append(total_protons)
                        plot_data_raw[energy][1].append(ks2_test(raw.data[total_protons], stats.binom.pmf(
                            range(0, total_protons+1), total_protons, 1.0/raw.div)))
                        plot_data_mix[energy][0].append(total_protons)
                        plot_data_mix[energy][1].append(ks2_test(mix.data[total_protons], stats.binom.pmf(
                            range(0, total_protons+1), total_protons, 1.0/mix.div)))
                        # if energy == cdf_plot['energy'] and total_protons == cdf_plot['total protons'] and \
                        #         division == cdf_plot['division'] and centrality == cdf_plot['centrality']:
                        #     plot_data_cdf['raw']
                        # raw_dist = [bin_protons for bin_protons, events in enumerate(raw.data[total_protons])
                        #             for l in range(events)]
                        # mix_dist = [bin_protons for bin_protons, events in enumerate(mix.data[total_protons])
                        #             for l in range(events)]
                        # ks2 = stats.ks_2samp(raw_dist, mix_dist)
                        # ks_raw = stats.kstest(raw_dist, 'binom',
                        #                       args=(total_protons, 1.0/raw.div))
                        # ks_mix = stats.kstest(mix_dist, 'binom',
                        #                       args=(total_protons, 1.0/mix.div))
                        # # print(total_protons)
                        # # plt.hist(mix_dist, bins=np.linspace(-0.5, total_protons+0.5, total_protons+2), align='mid', density=True, zorder=0)
                        # # plt.scatter(range(0, total_protons+1), stats.binom.pmf(range(0, total_protons+1), total_protons, 1.0/mix.div),
                        # #          zorder=1, color='red')
                        # # plt.show()
                        # plot_data_ks2[energy][0].append(total_protons)
                        # plot_data_ks2[energy][1].append(ks2[0])
                        # plot_data_raw[energy][0].append(total_protons)
                        # plot_data_raw[energy][1].append(ks_raw[0])
                        # plot_data_mix[energy][0].append(total_protons)
                        # plot_data_mix[energy][1].append(ks_mix[0])

    raw_cdf = get_cdf(data['raw'][cdf_plot['energy']][cdf_plot['division']][cdf_plot['centrality']]
                      .data[cdf_plot['total protons']])
    mix_cdf = get_cdf(data['mix'][cdf_plot['energy']][cdf_plot['division']][cdf_plot['centrality']]
                      .data[cdf_plot['total protons']])

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()

    for energy in energies:
        # ax1.plot(plot_data_ks2[energy][0], plot_data_ks2[energy][1], color=colors[energy], label=f'{energy}GeV')
        ax2.plot(plot_data_raw[energy][0], plot_data_raw[energy][1], color=colors[energy], label=f'{energy}GeV')
        ax3.plot(plot_data_mix[energy][0], plot_data_mix[energy][1], color=colors[energy], label=f'{energy}GeV')
        ax4.plot(plot_data_ks2_test[energy][0], plot_data_ks2_test[energy][1], color=colors[energy],
                 label=f'{energy}GeV')

    ax1.scatter(range(0, len(raw_cdf)), raw_cdf, label='Raw CDF')
    ax1.scatter(range(0, len(mix_cdf)), mix_cdf, label='Mix CDF')
    ax1.legend()
    ax1.set_title(f'{cdf_plot["energy"]}GeV, {cdf_plot["division"]} div, {cdf_plot["centrality"]} cent,'
                  f' {cdf_plot["total protons"]} protons, CDF Comparison')
    ax1.set_xlabel('Protons in Bin')
    ax1.set_ylabel('Integrated Probability')

    ax2.legend()
    ax2.set_title('KS Statistics for Data to Binomial | 3 Divisions Most Central')
    ax2.set_xlabel('Total Protons in Event')
    ax2.set_ylabel('KS Statistic')

    ax3.legend()
    ax3.set_title('KS Statistics for Mixed to Binomial | 3 Divisions Most Central')
    ax3.set_xlabel('Total Protons in Event')
    ax3.set_ylabel('KS Statistic')

    ax4.legend()
    ax4.set_title('Test Two Sample KS Statistics for Data to Mixed | 3 Divisions Most Central')
    ax4.set_xlabel('Total Protons in Event')
    ax4.set_ylabel('KS Statistic')

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


if __name__ == '__main__':
    main()
