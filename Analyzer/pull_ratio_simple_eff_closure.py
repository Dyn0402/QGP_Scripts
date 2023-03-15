#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 15 11:09 AM 2023
Created in PyCharm
Created as QGP_Scripts/pull_ratio_simple_eff_closure.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats
from Measure import Measure


def main():
    base_path = 'F:/Research/'
    dir_names = ['Data_Sim_2source_Tests',
                 'Data_Sim_2source_Tests']
    group_names = ['flat80_simpleclust_spread05_amp2_resample_epbins1_test',
                   'flat80_Eff_simpleclust_spread05_amp2_resample_epbins1_test']
    set_names = ['Sim_flat80_simpleclust_spread05_amp2_norotate_resample_epbins1_0',
                 'Sim_flat80_Eff_simpleclust_spread05_amp2_norotate_resample_epbins1_0']
    set_titles = ['Simple Clustering', 'Simple Clustering + Efficiency']
    # dir_names = ['Data_Sim_Flow',
    #              'Data_Sim_2source_Tests']
    # group_names = ['flow_flat80_res15_v207_resample',
    #                'flat80_Eff_flow_res15_v207_resample_epbins1_test']
    # set_names = ['Sim_flow_flat80_res15_v207_resample_norotate_0',
    #              'Sim_flat80_Eff_flow_res15_v207_norotate_resample_epbins1_0']
    # set_titles = ['Flow', 'Flow + Efficiency']
    energy = 62
    divs = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    cent = 8
    stats = {'kurtosis': lambda x: x.get_kurtosis(),
             'non-excess kurtosis': lambda x: x.get_non_excess_kurtosis(),
             'variance': lambda x: x.get_variance(),
             'standard deviation': lambda x: x.get_sd(),
             'skewness': lambda x: x.get_skewness()}

    for stat, stat_method in stats.items():
        fig_ratio, ax_ratio = plt.subplots(dpi=144)
        fig_pull, ax_pull = plt.subplots(dpi=144)

        for set_name, group_name, dir_name, set_title in zip(set_names, group_names, dir_names, set_titles):
            ratio_raw, ratio_mix, pull_raw, pull_mix = [], [], [], []
            for div in divs:
                set_dir_path = f'{group_name}/{set_name}/'
                file_path = f'{set_dir_path}{energy}GeV/ratios_divisions_{div}_centrality_{cent}_local.txt'
                path = f'{base_path}/{dir_name}/{file_path}'
                path_mix = f'{base_path}/{dir_name}_Mix/{file_path}'

                # print(path)
                # print(path_mix)

                az_data_raw = AzimuthBinData(div=div, path=path)
                az_data_mix = AzimuthBinData(div=div, path=path_mix)

                ratio_dist_raw = az_data_raw.get_ratio_dist()
                az_ratio_stats_raw = DistStats(dist=ratio_dist_raw)
                ratio_raw.append(stat_method(az_ratio_stats_raw))

                ratio_dist_mix = az_data_mix.get_ratio_dist()
                az_ratio_stats_mix = DistStats(dist=ratio_dist_mix)
                ratio_mix.append(stat_method(az_ratio_stats_mix))

                pull_dist_raw = az_data_raw.get_pull_dist()
                az_pull_stats_raw = DistStats(dist=pull_dist_raw)
                pull_raw.append(stat_method(az_pull_stats_raw))

                pull_dist_mix = az_data_mix.get_pull_dist()
                az_pull_stats_mix = DistStats(dist=pull_dist_mix)
                pull_mix.append(stat_method(az_pull_stats_mix))

            for ax, raw, mix in [(ax_ratio, ratio_raw, ratio_mix), (ax_pull, pull_raw, pull_mix)]:
                if 'Efficiency' in set_title:
                    ax.errorbar(divs, [x.val for x in raw], [x.err for x in raw], label=f'{set_title} Raw', marker='o',
                                alpha=0.6, color='blue')
                    ax.errorbar(divs, [x.val for x in mix], [x.err for x in mix], label=f'{set_title} Mix', marker='o',
                                alpha=0.6, color='green')
                    div = np.array(raw) / np.array(mix)
                    ax.errorbar(divs, [x.val for x in div], [x.err for x in div], label=f'{set_title} Div', marker='o',
                                alpha=0.6, color='red')
                    sub = np.array(raw) - np.array(mix)
                    ax.errorbar(divs, [x.val for x in sub], [x.err for x in sub], label=f'{set_title} Sub', marker='o',
                                alpha=0.6, color='purple')
                else:
                    ax.errorbar(divs, [x.val for x in raw], [x.err for x in raw], label=f'{set_title} Raw', marker='o',
                                alpha=0.6, color='orange')

        for fig, ax, title in [(fig_ratio, ax_ratio, 'Ratio'), (fig_pull, ax_pull, 'Pull')]:
            ax.grid()
            ax.set_ylabel(f'{stat}')
            ax.set_xlabel('Azimuthal Partition Width (w)')
            ax.set_title(title)
            ax.legend()
            fig.canvas.manager.set_window_title(title)
            fig.tight_layout()

    plt.show()

    print('dono')


if __name__ == '__main__':
    main()
