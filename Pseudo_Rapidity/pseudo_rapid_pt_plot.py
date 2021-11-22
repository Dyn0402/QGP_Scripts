#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 29 9:01 PM 2021
Created in PyCharm
Created as QGP_Scripts/BES_Pseudo_Rapid_Plot

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os

import uproot
import awkward as ak
import vector
import hist

from rapidity_vs_pseudo import rapidity


def main():
    m = 0.93827  # GeV, proton mass
    energy = '7GeV'
    vector.register_awkward()
    bes1(m, energy)
    ampt(m, energy)
    cf(m, energy)
    plt.show()

    print('donzo')


def bes1(m, energy):
    path = '/home/dylan/Research/BES1_Trees/' + energy + '/'
    tree_name = 'tree'
    track_atts = {'eta': 'proton.eta', 'pt': 'proton.pt'}
    files = 100

    h_rap_eta, h_rap_pt, h_eta_pt = make_hists()

    file_num = 0
    for file_name in os.listdir(path):
        print(file_name)
        if '.root' not in file_name:
            continue
        root_path = os.path.join(path, file_name)
        with uproot.open(root_path) as file:
            tracks = file[tree_name].arrays(track_atts.values())
            m_array = tracks[track_atts['pt']] * 0 + m
            rapid = rapidity(tracks[track_atts['eta']], tracks[track_atts['pt']], m_array)
            h_rap_eta.fill(rapid=ak.flatten(rapid), eta=ak.flatten(tracks[track_atts['eta']]))
            h_rap_pt.fill(rapid=ak.flatten(rapid), pt=ak.flatten(tracks[track_atts['pt']]))
            h_eta_pt.fill(eta=ak.flatten(tracks[track_atts['eta']]), pt=ak.flatten(tracks[track_atts['pt']]))
        file_num += 1
        if file_num >= files:
            break

    plot_hists(h_rap_eta, h_rap_pt, h_eta_pt, m, 'STAR ' + energy)


def ampt(m, energy):
    path = '/home/dylan/Research/AMPT_Trees/slim_most_central/string_melting/' + energy + '/'
    tree_name = 'tree'
    good_pid = 2212
    track_atts = ['pid', 'px', 'py', 'pz']
    files = 80

    h_rap_eta, h_rap_pt, h_eta_pt = make_hists()

    file_num = 0
    for file_name in os.listdir(path):
        print(file_name)
        if '.root' not in file_name:
            continue
        root_path = os.path.join(path, file_name)
        with uproot.open(root_path) as file:
            tracks = file[tree_name].arrays(track_atts)
            tracks = tracks[tracks.pid == good_pid]
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz'],
                             'M': tracks['px'] * 0 + m}, with_name='Momentum4D')
            h_rap_eta.fill(rapid=ak.flatten(tracks.rapidity), eta=ak.flatten(tracks.eta))
            h_rap_pt.fill(rapid=ak.flatten(tracks.rapidity), pt=ak.flatten(tracks.pt))
            h_eta_pt.fill(eta=ak.flatten(tracks.eta), pt=ak.flatten(tracks.pt))
        file_num += 1
        if file_num >= files:
            break

    plot_hists(h_rap_eta, h_rap_pt, h_eta_pt, m, 'AMPT ' + energy)


def cf(m, energy):
    path = '/home/dylan/Research/Cooper_Frye_Trees/' + energy + '/'
    tree_name = 'tree'
    good_pid = 2212
    track_atts = ['pid', 'px', 'py', 'pz']
    files = 1
    max_eta = 10
    min_pt = 0.3

    h_rap_eta, h_rap_pt, h_eta_pt = make_hists()

    file_num = 0
    for file_name in os.listdir(path):
        print(file_name)
        if '.root' not in file_name:
            continue
        root_path = os.path.join(path, file_name)
        with uproot.open(root_path) as file:
            tracks = file[tree_name].arrays(track_atts)
            tracks = tracks[tracks.pid == good_pid]
            tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz'],
                             'M': tracks['px'] * 0 + m}, with_name='Momentum4D')
            tracks = tracks[(abs(tracks.pt) >= min_pt) & (abs(tracks.eta) <= max_eta)]
            h_rap_eta.fill(rapid=ak.flatten(tracks.rapidity), eta=ak.flatten(tracks.eta))
            h_rap_pt.fill(rapid=ak.flatten(tracks.rapidity), pt=ak.flatten(tracks.pt))
            h_eta_pt.fill(eta=ak.flatten(tracks.eta), pt=ak.flatten(tracks.pt))
        file_num += 1
        if file_num >= files:
            break

    plot_hists(h_rap_eta, h_rap_pt, h_eta_pt, m, 'CF ' + energy)


def make_hists():
    h_rap_eta = (
        hist.Hist.new
            .Reg(100, -2.2, 2.2, name='rapid', label='Rapidity')
            .Reg(100, -2.2, 2.2, name='eta', label='Pseudo-Rapidity')
            .Int64()
    )

    h_rap_pt = (
        hist.Hist.new
            .Reg(100, -2.2, 2.2, name='rapid', label='Rapidity')
            .Reg(100, 0, 2.4, name='pt', label='Transverse Momentum')
            .Int64()
    )

    h_eta_pt = (
        hist.Hist.new
            .Reg(100, -2.2, 2.2, name='eta', label='Pseudo-Rapidity')
            .Reg(100, 0, 2.4, name='pt', label='Transverse Momentum')
            .Int64()
    )

    return h_rap_eta, h_rap_pt, h_eta_pt


def plot_hists(h_rap_eta, h_rap_pt, h_eta_pt, m, title=''):
    etas = [1.0, 1.4, 1.8, 2.1]
    pt_cut_low = 0.4
    pt_cut_high = 2.0

    shade = True
    to_plot = ['h_eta_pt']  # ['h_rap_eta', 'h_rap_pt', 'h_eta_pt']

    if 'h_rap_eta' in to_plot:
        fig1, ax1 = plt.subplots()
        h_rap_eta.plot2d(cmap='jet', norm=colors.LogNorm(vmin=1))
        ax1.axhline(pt_cut_low, color='red', ls='--')
        ax1.axhline(pt_cut_high, color='red', ls='--')
        ax1.grid()
        ax1.set_title(title + ' Rapidity vs Eta')
        fig1.canvas.manager.set_window_title(title + ' Rapidity vs Eta')
        fig1.tight_layout()

    if 'h_rap_pt' in to_plot:
        fig2, ax2 = plt.subplots()
        h_rap_pt.plot2d(cmap='jet', norm=colors.LogNorm(vmin=1))
        y = np.linspace(min(h_rap_pt.axes[1].edges), max(h_rap_pt.axes[1].edges), 1000)
        m_array = np.array([m] * len(y))
        for eta in etas:
            eta_array = np.array([eta]*len(y))
            x_plus = np.asarray(rapidity(eta_array, y, m_array))
            x_minus = np.asarray(rapidity(-eta_array, y, m_array))
            ax2.plot(x_plus, y, color='black', ls='--', label=f'|eta| = {eta}')
            ax2.plot(x_minus, y, color='black', ls='--')
            if eta == 1.4 and shade:
                ax2.fill_between(x_plus, y, pt_cut_low, where=(y > pt_cut_low) & (x_plus < 1), color='gray', alpha=0.6)
                ax2.fill_between(x_minus, y, pt_cut_low, where=(y > pt_cut_low) & (x_minus > -1), color='gray', alpha=0.6)
            elif eta == 1.8 and shade:
                ax2.fill_between(x_plus, y, pt_cut_low, where=(y > pt_cut_low) & (x_plus < 1), color='black', alpha=0.7)
                ax2.fill_between(x_minus, y, pt_cut_low, where=(y > pt_cut_low) & (x_minus > -1), color='black', alpha=0.7)
        ax2.axhline(pt_cut_low, color='blue', ls='--')
        ax2.axhline(pt_cut_high, color='red', ls='--')
        ax2.grid()
        ax2.legend()
        ax2.set_title(title + ' Rapidity vs Pt')
        fig2.canvas.manager.set_window_title(title + ' Rapidity vs Pt')
        fig2.tight_layout()

    if 'h_eta_pt' in to_plot:
        fig3, ax3 = plt.subplots()
        h_eta_pt.plot2d(cmap='jet', norm=colors.LogNorm(vmin=1))
        for eta in etas:
            ax3.axvline(eta, color='black', ls='--', label=f'|eta| = {eta}')
            ax3.axvline(-eta, color='black', ls='--')
        ax3.axhline(pt_cut_low, color='blue', ls='--')
        ax3.axhline(pt_cut_high, color='red', ls='--')
        ax3.grid()
        ax3.legend()
        ax3.set_title(title + ' Eta vs Pt')
        fig3.canvas.manager.set_window_title(title + ' Eta vs Pt')
        fig3.tight_layout()

    # plt.show()


if __name__ == '__main__':
    main()
