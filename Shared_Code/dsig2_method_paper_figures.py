#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 28 8:06 PM 2024
Created in PyCharm
Created as QGP_Scripts/dsig2_method_paper_figures.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf


def main():
    plot_method_paper_figs()
    print('donzo')


def plot_method_paper_figs():
    plt.rcParams["figure.dpi"] = 144

    # Gaussian correlation model simulations
    df_sim_width_fits = pd.read_csv(f'binom_slice_stats_sim_demos_width_fits.csv')

    plot_b_vs_amp(df_sim_width_fits, alpha=0.8, ylim=(-0.00046, 0.00053))

    # amp_shifts = {0.002: -0.01, 0.006: 0.0, 0.01: +0.01}
    # amp_colors = {0.002: 'red', 0.006: 'black', 0.01: 'blue'}
    # amp_markers = {0.002: 's', 0.006: 'o', 0.01: '^'}
    # plot_z_vs_spread(df_sim_width_fits, amps=list(amp_colors.keys()), amps_colors=amp_colors, amps_markers=amp_markers,
    #                  amps_x_shifts=amp_shifts, alpha=0.7)

    plt.show()


def plot_z_vs_spread(df_fits, amps=None, spreads=None, amps_colors=None, amps_markers=None, amps_x_shifts=None,
                     alpha=1.0):
    fig_zeros_spread, ax_zeros_spread = plt.subplots(figsize=(8, 4), dpi=144)
    ax_zeros_spread.set_xlabel(r'Simulation Range ($\sigma$)')
    ax_zeros_spread.set_ylabel(r'$z$')
    fig_zeros_spread.canvas.manager.set_window_title('Quad Fit Zeros vs Spread')

    if amps is None:
        amps = pd.unique(df_fits['amp'])
    if spreads is None:
        spreads = pd.unique(df_fits['spread'])

    print(f'Amps: {amps}\nSpreads: {spreads}')

    colors_amps = iter(plt.cm.rainbow(np.linspace(0, 1, len(amps))))
    for amp in amps:
        df_amp = df_fits[df_fits['amp'] == amp]
        df_cl_amp = df_amp[df_amp['data_set'].str.contains('_clmul_')]
        df_acl_amp = df_amp[df_amp['data_set'].str.contains('_aclmul_')]
        if amps_colors is None or amp not in amps_colors:
            color = next(colors_amps)
        else:
            color = amps_colors[amp]
        if amps_markers is None or amp not in amps_markers:
            marker = 'o'
        else:
            marker = amps_markers[amp]
        if amps_x_shifts is None or amp not in amps_x_shifts:
            x_shift = 0
        else:
            x_shift = amps_x_shifts[amp]
        ax_zeros_spread.errorbar(df_cl_amp['spread'] + x_shift, df_cl_amp['zero_mag'],
                                 yerr=df_cl_amp['zero_mag_err'], ls='none', alpha=alpha,
                                 marker=marker, label=rf'$A=+{amp}$', color=color, fillstyle='none')
        ax_zeros_spread.errorbar(df_acl_amp['spread'] + x_shift, df_acl_amp['zero_mag'],
                                 yerr=df_acl_amp['zero_mag_err'], ls='none', alpha=alpha,
                                 marker=marker, label=rf'$A=-{amp}$', color=color, fillstyle='full')

    ax_zeros_spread.legend()
    fig_zeros_spread.tight_layout()
    fig_zeros_spread.subplots_adjust(left=0.07, top=0.99, bottom=0.11, right=0.995)


def plot_b_vs_amp(df_fits, spreads=None, alpha=1.0, ylim=None):
    fig_base_amp, ax_base_amp = plt.subplots(figsize=(8, 4), dpi=144)
    ax_base_amp.set_xlabel(r'Simulation Amplitude Magnitude ($\left| A \right|$)')
    ax_base_amp.set_ylabel(r'$b$')
    ax_base_amp.axhline(0, color='black')
    ax_base_amp.axvline(0, color='black')
    fig_base_amp.canvas.manager.set_window_title('Quad Fit Baseline vs Amplitude')

    if spreads is None:
        spreads = pd.unique(df_fits['spread'])
    cl_type_data = {'_clmul_': {'name': 'Attractive', 'line': '--', 'fill': 'none'},
                    '_aclmul_': {'name': 'Repulsive', 'line': ':', 'fill': 'full'}}
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(spreads))))
    for spread in spreads:
        color = next(colors)
        df_spread = df_fits[df_fits['spread'] == spread]
        for cl_type, cl_data in cl_type_data.items():
            df_set = df_spread[df_spread['data_set'].str.contains(cl_type)]
            if df_set.size == 0 or any(np.isnan(df_set['baseline'])) or any(np.isinf(df_set['baseline'])) or \
                    any(np.isnan(df_set['base_err'])) or any(np.isinf(df_set['base_err'])):
                continue
            if cl_type == '_aclmul_':
                lab = f'σ={round(spread, 2)}'
            else:
                lab = None
            ax_base_amp.errorbar(df_set['amp'], df_set['baseline'], yerr=df_set['base_err'], ls='none',
                                 marker='o', label=lab, color=color, fillstyle=cl_data['fill'], alpha=alpha)
            popt, pcov = cf(line, df_set['amp'], df_set['baseline'], sigma=df_set['base_err'],
                            absolute_sigma=True)
            x_vals = np.array([0, max(df_set['amp'])])
            ax_base_amp.plot(x_vals, line(x_vals, *popt), ls=cl_data['line'], color=color)
    if ylim is not None:
        ax_base_amp.set_ylim(ylim)
    ax_base_amp.legend()
    fig_base_amp.tight_layout()
    fig_base_amp.subplots_adjust(left=0.115, top=0.99, bottom=0.11, right=0.995)


def line(x, a, b):
    return a * x + b


if __name__ == '__main__':
    main()
