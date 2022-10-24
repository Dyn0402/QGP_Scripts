#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 24 8:27 AM 2022
Created in PyCharm
Created as QGP_Scripts/make_cent_ref_file

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import os
import time

import ROOT
import awkward as ak
import vector
import uproot

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from DistStats import DistStats


def main():
    df_out_path = '/home/dylan/Research/Results/Azimuth_Analysis/mean_cent_ref.csv'
    df_bes = get_bes_cent_ref()
    df_ampt = get_ampt_cent_ref()
    df = pd.DataFrame(df_bes + df_ampt)
    df.to_csv(df_out_path, index=False)
    print('donzo')


def get_bes_cent_ref():
    energies = [7, 11, 19, 27, 39, 62]
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    path = '/home/dylan/Research/Data/qaonly/'
    set_name = 'rapid05_resample_dca1_nsprx1_m2r6_m2s0_nhfit20_qaonly_0'

    print('Get BES Refmult Averages per Centrality')
    df = []
    for energy in energies:
        qa_path = f'{path}{set_name}/{energy}GeV/QA_{energy}GeV.root'
        file = ROOT.TFile(qa_path, "READ")
        for cent in cents:
            hist_ref = file.Get(f'RefMult_{set_name}_{energy}_{cent}')
            mean_ref = hist_ref.GetMean()
            sd_ref = hist_ref.GetRMS()

            hist_refn = file.Get(f'RefMultn_{set_name}_{energy}_{cent}')
            mean_refn = hist_refn.GetMean()
            sd_refn = hist_refn.GetRMS()

            df.append({'data_set': 'bes_resample_def', 'energy': energy, 'cent': cent, 'mean_ref_val': mean_ref,
                       'mean_ref_sd': sd_ref, 'mean_refn_val': mean_refn, 'mean_refn_sd': sd_refn})
    return df


def get_ampt_cent_ref():
    energies = [7, 11, 19, 27, 39, 62]
    # energies = [27]
    cents = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    bin_edges = np.arange(-0.5, 801.5, 1)
    plot = False

    threads = 15
    ampt_tree_path = '/home/dylan/Research/AMPT_Trees_New_Coalescence/min_bias/string_melting/'
    ampt_cent_path = '/home/dylan/Research/Ampt_Centralities_New_Coalescence/string_melting/'

    print('Get Ampt Refmult Averages per Centrality')
    df = []
    for energy in energies:
        print(f'Starting {energy}GeV:')
        with open(f'{ampt_cent_path}{energy}GeV_Ampt_ref_bin_edge.txt', 'r') as file:
            edges = file.readlines()[2].strip().split()
            ref3_edges = ampt_str_edges_to_9(edges)

        print(ref3_edges)
        tree_names = os.listdir(f'{ampt_tree_path}{energy}GeV')
        jobs = [(f'{ampt_tree_path}{energy}GeV/{tree_name}', ref3_edges, bin_edges) for tree_name in tree_names]

        refmult_hists, refmult3_hists = [{cent: np.zeros(len(bin_edges) - 1) for cent in cents} for i in range(2)]
        with Pool(threads) as pool:
            for refmult_h, refmult3_h in tqdm.tqdm(pool.istarmap(read_ampt_tree, jobs), total=len(jobs)):
                for cent in cents:
                    refmult_hists[cent] += refmult_h[cent]
                    refmult3_hists[cent] += refmult3_h[cent]

        bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
        for cent in cents:
            ref_stats = DistStats(refmult_hists[cent])
            ref3_stats = DistStats(refmult3_hists[cent])
            ref_mean, ref_sd = ref_stats.get_mean().val, ref_stats.get_sd().val
            ref3_mean, ref3_sd = ref3_stats.get_mean().val, ref3_stats.get_sd().val
            if plot:
                print(f'Cent {cent}:\nref={ref_mean} +- {ref_sd}\n ref3={ref3_mean} +- {ref3_sd}')
                fig_ref, ax_ref = plt.subplots()
                ax_ref.bar(bin_centers, refmult_hists[cent], width=1, align='center')
                ax_ref.axvline(ref_mean, ls='--', color='red')
                ax_ref.axvspan(ref_mean - ref_sd, ref_mean + ref_sd, color='red', alpha=0.3)
                ax_ref.set_title(f'{energy}GeV Centrality {cent} Refmult')
                fig_ref.canvas.manager.set_window_title(f'{energy}GeV Centrality {cent} Refmult')

                fig_ref3, ax_ref3 = plt.subplots()
                ax_ref3.bar(bin_centers, refmult3_hists[cent], width=1, align='center')
                ax_ref3.axvline(ref3_mean, ls='--', color='red')
                ax_ref3.axvspan(ref3_mean - ref3_sd, ref3_mean + ref3_sd, color='red', alpha=0.3)
                ax_ref3.set_title(f'{energy}GeV Centrality {cent} Refmult3')
                fig_ref3.canvas.manager.set_window_title(f'{energy}GeV Centrality {cent} Refmult3')
            df.append({'data_set': 'ampt_new_coal_resample_def', 'energy': energy, 'cent': cent,
                       'mean_ref_val': ref_stats.get_mean().val, 'mean_ref_sd': ref_stats.get_sd().val,
                       'mean_refn_val': ref3_stats.get_mean().val, 'mean_refn_sd': ref3_stats.get_sd().val})
    if plot:
        plt.show()

    return df


def ampt_str_edges_to_9(ref_str_edges):
    edges = [int(edge) for edge in ref_str_edges]
    edges = sorted(edges, reverse=True)
    edges_9bin = {8: [1000, edges[0]],  # Set upper edge of most central to very high value
                  7: [edges[0], edges[1]]}  # Do 0-5% and 5-10% manually

    cent_bin, edge_index = 6, 1
    while edge_index + 2 < len(edges):
        edges_9bin.update({cent_bin: [edges[edge_index], edges[edge_index + 2]]})
        cent_bin -= 1
        edge_index += 2

    return edges_9bin


def read_ampt_tree(tree_path, ref3_edges, bin_edges):
    refmult_hists, refmult3_hists = [{cent: np.zeros(len(bin_edges) - 1) for cent in ref3_edges.keys()}
                                     for i in range(2)]
    with uproot.open(tree_path) as file:
        tracks = file['tree'].arrays(['refmult', 'refmult3'])
        for cent_bin, ref3_edges in ref3_edges.items():
            cent_bin_events = tracks[(tracks['refmult3'] >= ref3_edges[1]) & (tracks['refmult3'] < ref3_edges[0])]
            refmult_hists[cent_bin] += np.histogram(cent_bin_events['refmult'], bin_edges)[0]
            refmult3_hists[cent_bin] += np.histogram(cent_bin_events['refmult3'], bin_edges)[0]
    return refmult_hists, refmult3_hists


if __name__ == '__main__':
    main()
