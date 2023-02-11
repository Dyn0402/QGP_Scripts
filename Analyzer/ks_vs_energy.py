#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 06 8:34 PM 2020
Created in PyCharm
Created as QGP_Scripts/ks_vs_energy.py

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit as cf

from Measure import Measure
from AzimuthBinData import AzimuthBinData as AzData
from BootstrapAzBin import BootstrapAzBin as BsData
from DistStats import DistStats


def main():
    # from_hist_files()
    # from_hist_files_simple()
    from_dataframe()
    # flow_plots()
    print('donzo')


def from_dataframe():
    base_path = '/media/ucla/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'

    divs = 120
    energy = 39
    cent = 4
    samples = 72  # For title only
    data_set = 'AMPT_noprerot_rp'
    if data_set == 'AMPT':
        df_name = 'binom_slice_stats_cents.csv'
        data_set_name = 'ampt_new_coal_resample_def'
    elif data_set == 'AMPT_noprerot':
        df_name = 'binom_slice_stats_ampt_noprerot.csv'
        data_set_name = 'ampt_new_coal_resample_def_noprerot'
    elif data_set == 'AMPT_noprerot_rp':
        df_name = 'binom_slice_stats_ampt_noprerot_rp.csv'
        data_set_name = 'ampt_new_coal_resample_def_noprerot_rp'
    elif data_set == 'BES':
        df_name = 'binom_slice_stats_cents.csv'
        data_set_name = 'bes_resample_def'
    elif data_set == 'BES_epbins1':
        df_name = 'binom_slice_stats_bes_epbins1.csv'
        data_set_name = 'bes_resample_epbins1'
    elif data_set == 'Flow':
        df_name = 'binom_slice_stats_flow.csv'
        data_set_name = 'flow_resample_res15_v205'
    elif data_set == 'CF':
        df_name = 'binom_slice_stats_cents.csv'
        data_set_name = 'cf_resample_def'
    stat = 'standard deviation'

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '60-70%', 1: '70-80%'}

    df_path = base_path + df_name
    df = pd.read_csv(df_path)
    df = df.dropna()
    df = df[(df['name'] == data_set_name) & (df['divs'] == divs) & (df['energy'] == energy) & (df['stat'] == stat)]
    if cent:
        df = df[df['cent'] == cent]
    df = df.sort_values(by=['total_protons'])
    df_raw = df[df['data_type'] == 'raw']
    df_mix = df[df['data_type'] == 'mix']

    p = float(divs) / 360
    y_binom_mix = (np.asarray(df_mix['total_protons']) * p * (1 - p)) ** 0.5
    y_binom_raw = (np.asarray(df_raw['total_protons']) * p * (1 - p)) ** 0.5

    fig1, ax1 = plt.subplots(figsize=(6.66, 5), dpi=144)
    ax1.errorbar(df_raw['total_protons'], df_raw['val'], df_raw['err'], alpha=0.8, zorder=2, color='blue',
                 ls='', marker='o', label='Raw')
    ax1.errorbar(df_mix['total_protons'], df_mix['val'], df_mix['err'], alpha=0.8, zorder=1, color='green',
                 ls='', marker='o', label='Mix')
    ax1.plot(df_mix['total_protons'], y_binom_mix, color='red', alpha=0.8, zorder=0, label='Binomial')
    ax1.set_xlabel('Total Protons in Event')
    ax1.set_ylabel('Standard Deviation of Slice')
    ax1.set_title(f'{data_set} {energy}GeV, {cent_map[cent]} Centrality, {divs}° Partitions, '
                  f'{samples} Samples per Event')
    ax1.legend()
    fig1.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(6.66, 5), dpi=144)
    raw_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
    mix_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
    raw_ratio_vals, raw_ratio_errs = [x.val for x in raw_ratio], [x.err for x in raw_ratio]
    mix_ratio_vals, mix_ratio_errs = [x.val for x in mix_ratio], [x.err for x in mix_ratio]
    ax2.axhline(1, zorder=0, color='red', ls='--')
    ax2.errorbar(df_raw['total_protons'], raw_ratio_vals, raw_ratio_errs, ls='', marker='o',
                 zorder=2, color='blue', alpha=0.8, label='Raw / Binomial')
    ax2.errorbar(df_mix['total_protons'], mix_ratio_vals, mix_ratio_errs, ls='', marker='o',
                 zorder=1, color='green', alpha=0.8, label='Mix / Binomial')
    popt_raw, pcov_raw = cf(line_xy11, df_raw['total_protons'], raw_ratio_vals, sigma=raw_ratio_errs,
                            absolute_sigma=True)
    popt_mix, pcov_mix = cf(line_xy11, df_mix['total_protons'], mix_ratio_vals, sigma=mix_ratio_errs,
                            absolute_sigma=True)
    x_plot = np.array([0, max(max(df_raw['total_protons']), max(df_mix['total_protons']))])
    ax2.plot(x_plot, line_xy11(x_plot, *popt_raw), ls=':', color='blue', alpha=0.6)
    ax2.plot(x_plot, line_xy11(x_plot, *popt_mix), ls=':', color='green', alpha=0.6)
    ax2.annotate(f'raw slope: {popt_raw[0]:.6f}+-{np.sqrt(np.diag(pcov_raw))[0]:.6f}\n'
                 f'mix slope: {popt_mix[0]:.6f}+-{np.sqrt(np.diag(pcov_mix))[0]:.6f}',
                 xy=(0.05, 0.1), xycoords='axes fraction')
    ax2.set_xlabel('Total Protons in Event')
    ax2.set_ylabel('Standard Deviation Ratio')
    ax2.set_title(
        f'{data_set} {energy}GeV, {cent_map[cent]} Centrality, {divs}° Partitions, {samples} Samples per Event')
    ax2.legend()
    fig2.tight_layout()

    fig3, ax3 = plt.subplots(figsize=(6.66, 5), dpi=144)
    raw_diff = [Measure(val, err) - binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
    mix_diff = [Measure(val, err) - binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
    ax3.errorbar(df_raw['total_protons'], [x.val for x in raw_diff], [x.err for x in raw_diff], alpha=0.8, zorder=2,
                 color='blue', ls='', marker='o', label='Raw - Binomial')
    ax3.errorbar(df_mix['total_protons'], [x.val for x in mix_diff], [x.err for x in mix_diff], alpha=0.8, zorder=1,
                 color='green', ls='', marker='o', label='Mix - Binomial')
    ax3.axhline(0, color='red', ls='--')
    ax3.set_xlabel('Total Protons in Event')
    ax3.set_ylabel('Standard Deviation Difference')
    ax3.set_title(
        f'{data_set} {energy}GeV, {cent_map[cent]} Centrality, {divs}° Partitions, {samples} Samples per Event')
    ax3.legend()
    fig3.tight_layout()

    plt.show()


def flow_plots():
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_stats_flow.csv'
    divs = 120
    energy = 62
    samples = 72  # For title only
    stat = 'standard deviation'

    df_path = base_path + df_name
    df_all = pd.read_csv(df_path)
    df_all = df_all.dropna()
    df_all = df_all[(df_all['divs'] == divs) & (df_all['energy'] == energy) & (df_all['stat'] == stat)]

    data_set_names = ['flow_resample_res99_v207', 'flow_resample_res9_v207', 'flow_resample_res75_v207',
                      'flow_resample_res5_v207', 'flow_resample_res3_v207', 'flow_resample_res15_v207']
    data_set_name_labels = ['v2=0.07, resolution=0.99', 'v2=0.07, resolution=0.9', 'v2=0.07, resolution=0.75',
                            'v2=0.07, resolution=0.5', 'v2=0.07, resolution=0.3', 'v2=0.07, resolution=0.15']
    data_set_name_labels = dict(zip(data_set_names, data_set_name_labels))

    fig, axs = plt.subplots(2, 3, figsize=(13.33, 6), dpi=144, sharex=True, sharey=True)
    for ax_i, (data_set_name, ax) in enumerate(zip(data_set_names, axs.flat)):
        df = df_all[(df_all['name'] == data_set_name)]
        df = df.sort_values(by=['total_protons'])
        df_raw = df[df['data_type'] == 'raw']
        df_mix = df[df['data_type'] == 'mix']

        p = float(divs) / 360
        y_binom_mix = (np.asarray(df_mix['total_protons']) * p * (1 - p)) ** 0.5
        y_binom_raw = (np.asarray(df_raw['total_protons']) * p * (1 - p)) ** 0.5

        raw_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
        mix_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
        raw_ratio_vals, raw_ratio_errs = [x.val for x in raw_ratio], [x.err for x in raw_ratio]
        mix_ratio_vals, mix_ratio_errs = [x.val for x in mix_ratio], [x.err for x in mix_ratio]
        ax.axhline(1, zorder=0, color='black', ls='-')
        ax.errorbar(df_raw['total_protons'], raw_ratio_vals, raw_ratio_errs, ls='', marker='o',
                    zorder=2, color='blue', alpha=0.8, label='Raw / Binomial')
        ax.errorbar(df_mix['total_protons'], mix_ratio_vals, mix_ratio_errs, ls='', marker='o',
                    zorder=1, color='green', alpha=0.8, label='Mix / Binomial')
        popt_raw, pcov_raw = cf(line_xy11, df_raw['total_protons'], raw_ratio_vals, sigma=raw_ratio_errs,
                                absolute_sigma=True)
        popt_mix, pcov_mix = cf(line_xy11, df_mix['total_protons'], mix_ratio_vals, sigma=mix_ratio_errs,
                                absolute_sigma=True)
        x_plot = np.array([0, max(max(df_raw['total_protons']), max(df_mix['total_protons']))])
        ax.plot(x_plot, line_xy11(x_plot, *popt_raw), ls=':', color='blue', alpha=0.6)
        ax.plot(x_plot, line_xy11(x_plot, *popt_mix), ls=':', color='green', alpha=0.6)
        ax.annotate(f'raw slope: {popt_raw[0]:.6f}+-{np.sqrt(np.diag(pcov_raw))[0]:.6f}\n'
                    f'mix slope: {popt_mix[0]:.6f}+-{np.sqrt(np.diag(pcov_mix))[0]:.6f}',
                    xy=(0.01, 0.6), xycoords='axes fraction')
        # ax.text(0, 1.023, data_set_name_labels[data_set_name])
        ax.annotate(data_set_name_labels[data_set_name], xy=(0.01, 0.7), xycoords='axes fraction')
        if ax_i >= 3:
            ax.set_xlabel('Total Protons in Event')
        if ax_i in [0, 3]:
            ax.set_ylabel('Standard Deviation Ratio')
        if ax_i == 0:
            ax.legend()
    fig.suptitle(f'Flow Simulation, {divs}° Partitions, {samples} Samples per Event')
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)

    fig2, axs2 = plt.subplots(2, 3, figsize=(13.33, 6), dpi=144, sharex=True, sharey=True)
    for ax_i, (data_set_name, ax) in enumerate(zip(data_set_names, axs2.flat)):
        df = df_all[(df_all['name'] == data_set_name)]
        df = df.sort_values(by=['total_protons'])
        df_raw = df[df['data_type'] == 'raw']
        df_mix = df[df['data_type'] == 'mix']

        p = float(divs) / 360
        y_binom_mix = (np.asarray(df_mix['total_protons']) * p * (1 - p)) ** 0.5
        y_binom_raw = (np.asarray(df_raw['total_protons']) * p * (1 - p)) ** 0.5

        raw_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
        mix_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
        ratio_diff = [raw - mix + 1 for raw, mix in zip(raw_ratio, mix_ratio)]
        ax.errorbar(df_raw['total_protons'], [x.val for x in ratio_diff], [x.err for x in raw_ratio], ls='', marker='o',
                    zorder=2, color='blue', alpha=0.8, label='Raw')
        ax.axhline(1, zorder=0, color='black', ls='-')
        # ax.text(0, 1.023, data_set_name_labels[data_set_name])
        ax.annotate(data_set_name_labels[data_set_name], xy=(0.01, 0.73), xycoords='axes fraction')
        if ax_i >= 3:
            ax.set_xlabel('Total Protons in Event')
        if ax_i in [0, 3]:
            ax.set_ylabel('Standard Deviation Ratio')
        if ax_i == 0:
            ax.legend()
    fig2.suptitle(f'Flow Simulation, {divs}° Partitions, {samples} Samples per Event')
    fig2.tight_layout()
    fig2.subplots_adjust(wspace=0, hspace=0)

    data_set_names = ['flow_resample_res15_v207', 'flow_resample_res15_v206', 'flow_resample_res15_v205',
                      'flow_resample_res15_v204', 'flow_resample_res15_v203', 'flow_resample_res15_v202']
    data_set_name_labels = ['v2=0.07, resolution=0.15', 'v2=0.06, resolution=0.15', 'v2=0.05, resolution=0.15',
                            'v2=0.04, resolution=0.15', 'v2=0.03, resolution=0.15', 'v2=0.02, resolution=0.15']
    data_set_name_labels = dict(zip(data_set_names, data_set_name_labels))

    fig, axs = plt.subplots(2, 3, figsize=(13.33, 6), dpi=144, sharex=True, sharey=True)
    for ax_i, (data_set_name, ax) in enumerate(zip(data_set_names, axs.flat)):
        df = df_all[(df_all['name'] == data_set_name)]
        df = df.sort_values(by=['total_protons'])
        df_raw = df[df['data_type'] == 'raw']
        df_mix = df[df['data_type'] == 'mix']

        p = float(divs) / 360
        y_binom_mix = (np.asarray(df_mix['total_protons']) * p * (1 - p)) ** 0.5
        y_binom_raw = (np.asarray(df_raw['total_protons']) * p * (1 - p)) ** 0.5

        raw_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_raw['val'], df_raw['err'], y_binom_raw)]
        mix_ratio = [Measure(val, err) / binom for val, err, binom in zip(df_mix['val'], df_mix['err'], y_binom_mix)]
        raw_ratio_vals, raw_ratio_errs = [x.val for x in raw_ratio], [x.err for x in raw_ratio]
        mix_ratio_vals, mix_ratio_errs = [x.val for x in mix_ratio], [x.err for x in mix_ratio]
        ax.axhline(1, zorder=0, color='black', ls='-')
        ax.errorbar(df_raw['total_protons'], raw_ratio_vals, raw_ratio_errs, ls='', marker='o',
                    zorder=2, color='blue', alpha=0.8, label='Raw / Binomial')
        ax.errorbar(df_mix['total_protons'], mix_ratio_vals, mix_ratio_errs, ls='', marker='o',
                    zorder=1, color='green', alpha=0.8, label='Mix / Binomial')
        popt_raw, pcov_raw = cf(line_xy11, df_raw['total_protons'], raw_ratio_vals, sigma=raw_ratio_errs,
                                absolute_sigma=True)
        popt_mix, pcov_mix = cf(line_xy11, df_mix['total_protons'], mix_ratio_vals, sigma=mix_ratio_errs,
                                absolute_sigma=True)
        x_plot = np.array([0, max(max(df_raw['total_protons']), max(df_mix['total_protons']))])
        ax.plot(x_plot, line_xy11(x_plot, *popt_raw), ls=':', color='blue', alpha=0.6)
        ax.plot(x_plot, line_xy11(x_plot, *popt_mix), ls=':', color='green', alpha=0.6)
        ax.annotate(f'raw slope: {popt_raw[0]:.6f}+-{np.sqrt(np.diag(pcov_raw))[0]:.6f}\n'
                    f'mix slope: {popt_mix[0]:.6f}+-{np.sqrt(np.diag(pcov_mix))[0]:.6f}',
                    xy=(0.01, 0.6), xycoords='axes fraction')
        # ax.text(0, 1.023, data_set_name_labels[data_set_name])
        ax.annotate(data_set_name_labels[data_set_name], xy=(0.01, 0.73), xycoords='axes fraction')
        if ax_i >= 3:
            ax.set_xlabel('Total Protons in Event')
        if ax_i in [0, 3]:
            ax.set_ylabel('Standard Deviation Ratio')
        if ax_i == 0:
            ax.legend()
    fig.suptitle(f'Flow Simulation, {divs}° Partitions, {samples} Samples per Event')
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)

    plt.show()


def from_hist_files_simple():
    base_path = 'F:/Research/'
    energy = 39
    divs = [60, 180]
    cent = 8
    min_bs = 200  # of 250
    # data_name = 'Ampt_New_Coal'
    # set_group = 'default_resample'
    # set_name = 'Ampt_rapid05_resample_norotate_'
    data_name = 'CF'
    set_group = 'default_resample'
    set_name = 'CF_rapid05_resample_norotate_'
    set_number = 0

    plot_data = []
    for div in divs:
        file_name = f'ratios_divisions_{div}_centrality_{cent}_local.txt'

        raw_path = f'{base_path}Data_{data_name}/{set_group}/{set_name}{set_number}/{energy}GeV/{file_name}'
        mix_path = f'{base_path}Data_{data_name}_Mix/{set_group}/{set_name}{set_number}/{energy}GeV/{file_name}'
        raw = BsData(div, raw_path)
        mix = BsData(div, mix_path)

        raw_stats = {tp: DistStats(x) for tp, x in raw.data.data.items()}
        mix_stats = {tp: DistStats(x) for tp, x in mix.data.data.items()}
        raw_bs_stats = [{tp: DistStats(y) for tp, y in x.data.items()} for x in raw.data_bs]
        mix_bs_stats = [{tp: DistStats(y) for tp, y in x.data.items()} for x in mix.data_bs]

        raw_tps, raw_c2_div_c1_vals, raw_c2_div_c1_errs = [], [], []
        for tp, raw_stat in raw_stats.items():
            val = raw_stat.get_cumulant(2).val / raw_stat.get_mean().val
            bs_c2_div_c1s = []
            for bs in raw_bs_stats:
                if tp in bs:
                    bs_c2_div_c1s.append(bs[tp].get_cumulant(2).val / bs[tp].get_cumulant(1).val)
            if len(bs_c2_div_c1s) > min_bs:
                raw_tps.append(tp)
                raw_c2_div_c1_vals.append(val)
                raw_c2_div_c1_errs.append(np.std(bs_c2_div_c1s))

        mix_tps, mix_c2_div_c1_vals, mix_c2_div_c1_errs = [], [], []
        for tp, mix_stat in mix_stats.items():
            val = mix_stat.get_cumulant(2).val / mix_stat.get_mean().val
            bs_c2_div_c1s = []
            for bs in mix_bs_stats:
                if tp in bs:
                    bs_c2_div_c1s.append(bs[tp].get_cumulant(2).val / bs[tp].get_cumulant(1).val)
            if len(bs_c2_div_c1s) > min_bs:
                mix_tps.append(tp)
                mix_c2_div_c1_vals.append(val)
                mix_c2_div_c1_errs.append(np.std(bs_c2_div_c1s))

        q = 1 - (float(div) / 360)
        raw_y_vals = (np.array(raw_c2_div_c1_vals) - q) / np.array(raw_tps)
        raw_y_errs = np.array(raw_c2_div_c1_errs) / np.array(raw_tps)
        mix_y_vals = (np.array(mix_c2_div_c1_vals) - q) / np.array(mix_tps)
        mix_y_errs = np.array(mix_c2_div_c1_errs) / np.array(mix_tps)

        raw_mix_tps, raw_mix_diff_vals, raw_mix_diff_errs = [], [], []
        raw_mix_diff_mean = 0
        for raw_index, raw_tp in enumerate(raw_tps):
            if raw_tp in mix_tps:
                mix_index = mix_tps.index(raw_tp)
                raw_mix_tps.append(raw_tp)
                raw_mix_diff = Measure(raw_y_vals[raw_index], raw_y_errs[raw_index]) - \
                               Measure(mix_y_vals[mix_index], mix_y_errs[mix_index])
                raw_mix_diff_mean += raw_mix_diff.val * 1 / raw_mix_diff.err ** 2
                raw_mix_diff_vals.append(raw_mix_diff.val)
                raw_mix_diff_errs.append(raw_mix_diff.err)

        raw_mix_diff_mean /= np.sum(1 / np.array(raw_mix_diff_errs) ** 2)

        plot_data.append((q, raw_tps, raw_c2_div_c1_vals, raw_c2_div_c1_errs, mix_tps, mix_c2_div_c1_vals,
                          mix_c2_div_c1_errs, raw_y_vals, raw_y_errs, mix_y_vals, mix_y_errs,
                          raw_mix_tps, raw_mix_diff_vals, raw_mix_diff_errs, raw_mix_diff_mean))

    fig1, ax1 = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig2, ax2 = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig3, ax3 = plt.subplots(figsize=(6.66, 5), dpi=144)
    markers = iter(['o', 's', 'p'])
    colors = iter(['b', 'g', 'r', 'm', 'c', 'k'])
    for div, plot_data in zip(divs, plot_data):
        marker, color = next(markers), next(colors)
        q, raw_tps, raw_c2_div_c1_vals, raw_c2_div_c1_errs, mix_tps, mix_c2_div_c1_vals, mix_c2_div_c1_errs, \
        raw_y_vals, raw_y_errs, mix_y_vals, mix_y_errs, raw_mix_tps, raw_mix_diff_vals, \
        raw_mix_diff_errs, raw_mix_diff_mean = plot_data

        ax1.errorbar(raw_tps, raw_c2_div_c1_vals, raw_c2_div_c1_errs, marker=marker, ls='none', color='blue', alpha=0.8,
                     label=f'Raw {div}°')
        ax1.errorbar(mix_tps, mix_c2_div_c1_vals, mix_c2_div_c1_errs, marker=marker, ls='none', color='green',
                     alpha=0.8, label=f'Mix {div}°')
        ax1.axhline(q, color='red', ls='--', label=f'Binomial {div}°')

        ax2.errorbar(raw_tps, raw_y_vals, raw_y_errs, marker=marker, ls='none', color='blue', alpha=0.8,
                     label=f'Raw {div}°')
        ax2.errorbar(mix_tps, mix_y_vals, mix_y_errs, marker=marker, ls='none', color='green', alpha=0.8,
                     label=f'Mix {div}°')

        ax3.errorbar(raw_mix_tps, raw_mix_diff_vals, raw_mix_diff_errs, alpha=0.8, ls='none', marker=marker,
                     color=color, label=f'{div}°')
        ax3.axhline(raw_mix_diff_mean, color=color, alpha=0.73)
        # ax3.axhspan(raw_mix_diff_mean.val - raw_mix_diff_mean.err, raw_mix_diff_mean.val + raw_mix_diff_mean.err,
        #             color=color, alpha=0.4)
        print(f'{div}° Mean: {raw_mix_diff_mean}')

    ax1.set_xlabel('Total Protons in Event')
    ax1.set_ylabel('C2 / C1')
    ax1.legend()
    fig1.tight_layout()

    ax2.axhline(0, color='black', ls='--')
    ax2.set_xlabel('Total Protons in Event')
    ax2.set_ylabel('(C2 / C1 - q) / n')
    ax2.legend()
    fig2.tight_layout()

    ax3.axhline(0, color='black', ls='--')
    ax3.set_xlabel('Total Protons in Event')
    ax3.set_ylabel('(C2 / C1 - q) / n - Mix')
    ax3.legend()
    fig3.tight_layout()

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
                    current_path = path[source][0] + str(energy) + path[source][1] + str(division) + path[source][2] + \
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
                        plot_data_ks2_test[energy][1].append(
                            ks2_test(raw.data[total_particles], mix.data[total_particles]))
                        plot_data_raw[energy][0].append(total_particles)
                        plot_data_raw[energy][1].append(ks2_test(raw.data[total_particles], stats.binom.pmf(
                            range(0, total_particles + 1), total_particles, 1.0 / raw.div)))
                        plot_data_mix[energy][0].append(total_particles)
                        plot_data_mix[energy][1].append(ks2_test(mix.data[total_particles], stats.binom.pmf(
                            range(0, total_particles + 1), total_particles, 1.0 / mix.div)))
                        sd_data_raw[energy][0].append(total_particles)
                        sd_data_raw[energy][1].append(
                            hist_sd(range(len(raw.data[total_particles])), raw.data[total_particles]))
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
    fig5, ax5 = plt.subplots(figsize=(6.66, 5), dpi=144)
    fig6, ax6 = plt.subplots(figsize=(6.66, 5), dpi=144)

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
    ax5.set_title(f'Standard Deviation of Total Particle Slices for {sd_plot["energy"]}GeV' + title_sufx)
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
    ax6.set_title(f'SD Divided by Binomial of Total Particle Slices for {sd_plot["energy"]}GeV' + title_sufx)
    ax6.legend()

    fig5.tight_layout()
    fig6.tight_layout()

    plt.show()

    print('donzo')


def ks2_test(dist1, dist2):
    ks = -1
    cdf1 = get_cdf(dist1)
    cdf2 = get_cdf(dist2)

    for i, j in zip(cdf1, cdf2):
        if abs(i - j) > ks:
            ks = abs(i - j)

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
        variance = sum((np.asarray(x) - mean) ** 2 * np.asarray(y)) / (sum(y) - 1)
        return variance ** 0.5


def hist_mean(x, y):
    mean = sum(np.asarray(x) * np.asarray(y)) / sum(y)
    return mean


def line_xy11(x, m):
    return m * (x - 1) + 1


if __name__ == '__main__':
    main()
