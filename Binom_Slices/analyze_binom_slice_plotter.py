#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 17 4:27 PM 2022
Created in PyCharm
Created as QGP_Scripts/analyze_binom_slice_plotter

@author: Dylan Neff, Dyn04
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from scipy.optimize import basinhopping
import warnings

import itertools
import inspect

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

from analyze_binom_slices import *
from presentation_plots import method_paper_plot
from pdf_examples import cluster_voids_example_plots
from sub_event_resample_algorithm import plot_method_paper_event
from momentum_conservation_model import plot_full_test_from_file

from integrate_pdf_var import base_gaus_pdf_wrap, get_partition_variance


def main():
    # plot_paper_figs()
    # plot_qm_figs()
    # plot_method_paper_figs()

    # plot_star_model_var()  # This doesn't do bkg subtraction error correction
    # plot_vs_cent_var()
    # plot_sims_var()
    # get_sim_mapping_var()
    # plot_all_zero_base()

    plot_star_var_sys()
    # make_models_csv()
    # plot_star_var_rand_sys()

    # plot_vs_cent_var_fits()
    # plot_vs_cent_var_fit_tests()
    # plot_vs_cent_var_fit_tests2()
    # plot_vs_cent_var_fit_tests3()

    # plot_sims()
    # get_sim_mapping()
    # get_sim_mapping_pm()
    # plot_star_model()
    # plot_star_model_onediv()
    # plot_vs_cent()
    # plot_closest_sims()
    # plot_vs_cent_nofit()
    # plot_vs_cent_fittest()
    # plot_ampt_efficiency()
    # plot_ampt_efficiency_var()
    # plot_flow()
    # plot_flow_k2()
    # plot_ampt_v2_closure()
    # plot_ampt_v2_closure_var()
    # plot_flow_v2_closure()
    # plot_flow_v2_closure_raw()
    # plot_flow_eff_test()
    # plot_anticl_flow_closure_test()
    # plot_anticl_flow_closure_test_simple()
    # plot_efficiency_closure_tests()
    # plot_closure_tests()

    # plot_raw_mix_sys_comp()
    print('donzo')


def plot_qm_figs():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144

    presentation_mode = True
    if presentation_mode:
        # plt.rcParams['axes.labelsize'] = 14  # Adjust the value as needed
        # plt.rcParams['axes.titlesize'] = 16  # Adjust the value as needed
        # plt.rcParams['legend.fontsize'] = 14  # Adjust the value as needed
        plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1

    # base_path = 'F:/Research/Results/Azimuth_Analysis/'
    base_path = 'C:/Users/Dyn04/OneDrive - personalmicrosoftsoftware.ucla.edu/OneDrive - UCLA IT Services/Research/UCLA/Results/Azimuth_Analysis/'
    df_dsigma_v2sub_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma_v2sub.csv'
    df_dsigma_v2sub_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma_v2sub.csv'
    df_def_avgs_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_bes.csv'
    df_def_avgs_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_model.csv'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [89, 356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_def', 'ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))
    data_sets_markers = dict(zip(data_sets_plt, [dict(zip(['raw', 'mix', 'diff'], [x, x, x]))
                                                 for x in ['o', 's', '^', '*']]))
    data_sets_bands = ['ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    legend_order = ['STAR', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$', 'AMPT']

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    df_dsigma_v2sub = pd.read_csv(f'{base_path}{df_dsigma_v2sub_name}')
    df_dsigma_v2sub_model = pd.read_csv(f'{base_path}{df_dsigma_v2sub_model_name}')
    df_dsigma_v2sub_model = df_dsigma_v2sub_model[df_dsigma_v2sub_model['err'] < 0.0001]
    df_dsigma_v2sub = pd.concat([df_dsigma_v2sub, df_dsigma_v2sub_model])

    dvar_vs_protons(df_dsigma_v2sub, div_plt, cent_plt, [39], ['diff'], data_sets_plt, ylabel=r'$\Delta\sigma^2$',
                    data_sets_colors=data_sets_colors, plot=True, avg=False, alpha=1.0, y_ranges=[-0.00124, 0.0009],
                    data_sets_labels=data_sets_labels, marker_map=data_sets_markers, legend_pos='lower right',
                    kin_info_loc=(0.22, 0.13), star_prelim_loc=(0.68, 0.35), data_sets_bands=data_sets_bands,
                    legend_order=legend_order, title=f'{cent_map[8]} Centrality, {div_plt}° Partitions',
                    xlim=(0, 43))  # <---

    dsig_avgs_all = pd.read_csv(f'{base_path}{df_def_avgs_out_name}')
    dsig_avgs_all_model = pd.read_csv(f'{base_path}{df_def_avgs_out_model_name}')
    dsig_avgs_all = pd.concat([dsig_avgs_all, dsig_avgs_all_model])

    dsig_avgs_v2sub = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_name}')
    dsig_avgs_v2sub_model = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_model_name}')
    dsig_avgs_v2sub = pd.concat([dsig_avgs_v2sub, dsig_avgs_v2sub_model])

    dsig_avgs_v2sub['data_type'] = 'diff'

    plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.2
    dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], data_sets_plt,
                             plot=True, avg=False, plot_avg=False, data_sets_colors=data_sets_colors, no_hydro_label=1,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.0003], avgs_df=dsig_avgs_v2sub,
                             ylabel=r'$\Delta\sigma^2$', ylim=(-0.0014, 0.00035),
                             # kin_loc=(0.65, 0.2), star_prelim_loc=(1, 0.54, 0.78),
                             kin_loc=(0.58, 0.45), star_prelim_loc=(2, 0.55, 0.2),
                             marker_map=data_sets_markers, data_sets_bands=data_sets_bands, legend_order=legend_order,
                             # title=f'{cent_map[cent_plt]} Centrality, 120° Partitions'
                             title=f''
                             )  # <---

    dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], data_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors, no_hydro_label=1,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.0003], avgs_df=dsig_avgs_v2sub,
                             ylabel=r'$\Delta\sigma^2$', ylim=(-0.0014, 0.00035),
                             # kin_loc=(0.65, 0.2), star_prelim_loc=(1, 0.54, 0.78),
                             kin_loc=(0.58, 0.45), star_prelim_loc=(2, 0.55, 0.2),
                             marker_map=data_sets_markers, data_sets_bands=data_sets_bands, legend_order=legend_order,
                             # title=f'{cent_map[cent_plt]} Centrality, 120° Partitions'
                             title=f''
                             )  # <---

    dsig_avgs_v2_sub_cent8_div120 = dsig_avgs_v2sub[(dsig_avgs_v2sub['cent'] == 8) & (dsig_avgs_v2sub['divs'] == 120)]
    plot_protons_avgs_vs_energy(dsig_avgs_v2_sub_cent8_div120, data_sets_plt, data_sets_colors=data_sets_colors,  # <---
                                data_sets_labels=data_sets_labels, alpha=1, kin_info_loc=(0.73, 0.45),
                                star_prelim_loc=(0.02, 0.75), marker_map=data_sets_markers,
                                data_sets_bands=data_sets_bands, legend_order=legend_order,
                                # title=f'{cent_map[8]} Centrality, {div_plt}° Partitions, {samples} Samples per Event')
                                title=f'{cent_map[8]} Centrality, {div_plt}° Partitions')
    plot_protons_avgs_vs_energy(dsig_avgs_v2_sub_cent8_div120, ['bes_def'], data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, alpha=1, kin_info_loc=(0.123, 0.68),
                                star_prelim_loc=(0.65, 0.9), leg_loc='lower right',
                                title=f'{cent_map[8]} Centrality, {div_plt}° Partitions, {samples} Samples per Event')

    dsig_avgs_v2_sub_div120 = dsig_avgs_v2sub[dsig_avgs_v2sub['divs'] == 120]
    data_sets_energies_colors = \
        {'bes_def': {7: 'red', 11: 'blue', 19: 'green', 27: 'orange', 39: 'purple', 62: 'black'},
         'ampt_new_coal_epbins1': {7: 'red', 11: 'blue', 19: 'green', 27: 'orange', 39: 'purple', 62: 'black'}}
    vs_cent_sets = list(data_sets_energies_colors.keys())
    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, ['bes_def'], data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.6, 0.5),
                              data_sets_energies_colors=data_sets_energies_colors)
    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, ['bes_def'], data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type='npart',  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.6, 0.5),
                              data_sets_energies_colors=data_sets_energies_colors)

    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, vs_cent_sets, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.3, 0.5), marker_map=data_sets_markers,
                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, vs_cent_sets, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type='npart',  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.3, 0.5), marker_map=data_sets_markers,
                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)

    data_sets_cent = ['ampt_new_coal_epbins1', 'bes_def']
    legend_order = ['7.7 GeV', '11.5 GeV', '19.6 GeV', '27 GeV', '39 GeV', '62.4 GeV', 'AMPT Fit']
    plot_dsig_avg_vs_cent_2panel(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
                                 cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
                                 # title=f'{div_plt}° Partitions, {samples} Samples per Event',
                                 title='',
                                 errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
                                 kin_info_loc=(0.45, 0.1), star_prelim_loc=(0.4, 0.5), marker_map=data_sets_markers,
                                 data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    plot_dsig_avg_vs_cent_2panel2(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
                                  cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
                                  # title=f'{div_plt}° Partitions, {samples} Samples per Event',
                                  title='',
                                  errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
                                  kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.62, 0.4), marker_map=data_sets_markers,
                                  data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.3
    plot_dsig_avg_vs_cent_2panel3(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
                                  cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
                                  # title=f'{div_plt}° Partitions, {samples} Samples per Event',
                                  title='',
                                  errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
                                  kin_info_loc=(0.3, 0.1), star_prelim_loc=(0.25, 0.3), marker_map=data_sets_markers,
                                  data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    # plt.savefig('C:/Users/Dyn04/Desktop/test.pdf', format='pdf')
    plt.show()


def plot_paper_figs():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144

    presentation_mode = False
    if presentation_mode:
        # plt.rcParams['axes.labelsize'] = 14  # Adjust the value as needed
        # plt.rcParams['axes.titlesize'] = 16  # Adjust the value as needed
        # plt.rcParams['legend.fontsize'] = 14  # Adjust the value as needed
        plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.2
    # else:
    #     plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.1

    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'C:/Users/Dyn04/OneDrive - personalmicrosoftsoftware.ucla.edu/OneDrive - UCLA IT Services/Research/UCLA/Results/Azimuth_Analysis/'
    # base_path = 'C:/Users/Dyn04/Research/'
    df_name = 'Bes_with_Sys/binom_slice_vars_bes.csv'
    df_model_name = 'Bes_with_Sys/binom_slice_vars_model.csv'
    df_dsigma_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma.csv'
    df_dsigma_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma.csv'
    df_dsigma_v2sub_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma_v2sub.csv'
    df_dsigma_v2sub_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma_v2sub.csv'
    df_def_avgs_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_bes.csv'
    df_def_avgs_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_model.csv'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'
    df_partitions_fits_name = 'Bes_with_Sys/partition_width_fits_bes.csv'
    df_partitions_fits_model_name = 'Bes_with_Sys/partition_width_fits_model.csv'

    save_dir = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Presentations/STAR_Paper/Autogen_Figures/'

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [89, 356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_def', 'ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))
    data_sets_markers = dict(zip(data_sets_plt, [dict(zip(['raw', 'mix', 'diff'], [x, x, x]))
                                                 for x in ['o', 's', '^', '*']]))
    data_sets_bands = ['ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    legend_order = ['STAR', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$', 'AMPT']

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    df = pd.read_csv(f'{base_path}{df_name}')
    df_model = pd.read_csv(f'{base_path}{df_model_name}')
    df = pd.concat([df, df_model])
    df = df[df['stat'] == stat_plot]

    df_dsigma = pd.read_csv(f'{base_path}{df_dsigma_name}')
    df_dsigma_model = pd.read_csv(f'{base_path}{df_dsigma_model_name}')
    df_dsigma = pd.concat([df_dsigma, df_dsigma_model])

    dvar_vs_protons(df_dsigma, div_plt, cent_plt, [39], ['raw', 'mix', 'diff'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00085, 0.00055], print_data=True,
                    data_sets_labels={'bes_def': {'raw': 'Single Event', 'mix': 'Mixed Event',
                                                  'diff': 'Single Event - Mixed Event'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^'}})  # Mix subtraction demo
    dvar_vs_protons(df_dsigma, div_plt, 4, [39], ['raw', 'mix', 'diff'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00085, 0.00055],
                    data_sets_labels={'bes_def': {'raw': 'Single Event', 'mix': 'Mixed Event',
                                                  'diff': 'Single Event - Mixed Event'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^'}})  # Mix subtraction demo

    dvar_vs_protons(df_dsigma, div_plt, cent_plt, [62], ['raw', 'mix', 'diff'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00085, 0.00055],
                    data_sets_labels={'bes_def': {'raw': 'Single Event', 'mix': 'Mixed Event',
                                                  'diff': 'Single Event - Mixed Event'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^'}})  # Mix subtraction demo

    # dvar_vs_protons(df_dsigma, div_plt, cent_plt, [39], ['diff'], data_sets_plt, data_sets_colors=data_sets_colors,
    #                 plot=True, avg=False, alpha=1.0, y_ranges=[-0.00124, 0.0009], ylabel=r'$\Delta\sigma^2$',
    #                 data_sets_labels=data_sets_labels, marker_map=data_sets_markers, legend_pos='lower right',
    #                 data_sets_bands=data_sets_bands)

    df_dsigma_v2sub = pd.read_csv(f'{base_path}{df_dsigma_v2sub_name}')
    df_dsigma_v2sub_model = pd.read_csv(f'{base_path}{df_dsigma_v2sub_model_name}')
    df_dsigma_v2sub_model = df_dsigma_v2sub_model[df_dsigma_v2sub_model['err'] < 0.0001]
    df_dsigma_v2sub = pd.concat([df_dsigma_v2sub, df_dsigma_v2sub_model])

    dvar_vs_protons(df_dsigma_v2sub, div_plt, cent_plt, [39], ['diff'], data_sets_plt, ylabel=r'$\Delta\sigma^2$',
                    data_sets_colors=data_sets_colors, plot=True, avg=False, alpha=1.0, y_ranges=[-0.00124, 0.0009],
                    data_sets_labels=data_sets_labels, marker_map=data_sets_markers, legend_pos='lower right',
                    kin_info_loc=(0.22, 0.13), star_prelim_loc=(0.65, 0.96), data_sets_bands=data_sets_bands,
                    legend_order=legend_order)  # 39 GeV mix and v2 subtract dsig2

    df_dsigma_v2sub_diffs = df_dsigma_v2sub[df_dsigma_v2sub['data_type'] == 'diff'].assign(data_type='v2_sub')
    df_dsigma_with_v2sub = pd.concat([df_dsigma, df_dsigma_v2sub_diffs])
    dvar_vs_protons(df_dsigma_with_v2sub, div_plt, cent_plt, [39], ['raw', 'mix', 'diff', 'v2_sub'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00124, 0.0009],
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^', 'v2_sub': '*'}})

    dvar_vs_protons(df_dsigma_with_v2sub, 90, 4, [62], ['raw', 'diff', 'v2_sub'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, ylabel=r'$\Delta\sigma^2$',
                    data_sets_labels={'bes_def': {'raw': 'Uncorrected', 'diff': 'Mixed Corrected',
                                                  'v2_sub': 'Mixed and Flow Corrected'}},
                    data_sets_colors={'bes_def': {'raw': 'blue', 'diff': 'red', 'v2_sub': 'black'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^', 'v2_sub': '*'}},
                    # y_ranges=[-0.00124, 0.0009])  # v2 sub demo
                    y_ranges=[-0.0039, 0.0009], kin_info_loc=(0.26, 0.94))  # v2 sub demo
    dvar_vs_protons(df_dsigma_with_v2sub, div_plt, 4, [39], ['raw', 'diff', 'v2_sub'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, ylabel=r'$\Delta\sigma^2$', print_data=True,
                    data_sets_labels={'bes_def': {'raw': 'Uncorrected', 'diff': 'Mixed Corrected',
                                                  'v2_sub': 'Mixed and Flow Corrected'}},
                    data_sets_colors={'bes_def': {'raw': 'blue', 'diff': 'red', 'v2_sub': 'black'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^', 'v2_sub': '*'}},
                    # y_ranges=[-0.00124, 0.0009])  # v2 sub demo
                    y_ranges=[-0.0039, 0.0009], kin_info_loc=(0.26, 0.94))  # v2 sub demo

    # plt.show()

    dsig_avgs_all = pd.read_csv(f'{base_path}{df_def_avgs_out_name}')
    dsig_avgs_all_model = pd.read_csv(f'{base_path}{df_def_avgs_out_model_name}')
    dsig_avgs_all = pd.concat([dsig_avgs_all, dsig_avgs_all_model])

    # dvar_vs_protons_energies(df_dsigma, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], ['bes_def'],
    #                          plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.0005], avgs_df=dsig_avgs_all,
    #                          ylabel=r'$\Delta\sigma^2$', kin_loc=(0.65, 0.2), legend_order=legend_order)

    dsig_avgs_v2sub = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_name}')
    dsig_avgs_v2sub_model = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_model_name}')
    dsig_avgs_v2sub = pd.concat([dsig_avgs_v2sub, dsig_avgs_v2sub_model])

    dsig_avgs_v2sub['data_type'] = 'diff'
    print(df_dsigma_v2sub.columns)
    df_dsigma_v2sub = df_dsigma_v2sub[df_dsigma_v2sub['err'] < 0.0001]
    # dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], ['bes_def'],
    #                          plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.0005], avgs_df=dsig_avgs_v2sub,
    #                          ylabel=r'$\Delta\sigma^2$',
    #                          # kin_loc=(0.55, 0.6), star_prelim_loc=(1, 0.54, 0.53)
    #                          kin_loc=(0.58, 0.45), star_prelim_loc=(1, 0.54, 0.4)
    #                          )

    # dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], data_sets_plt,
    #                          plot=True, avg=False, plot_avg=False, data_sets_colors=data_sets_colors, no_hydro_label=1,
    #                          data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.0005], avgs_df=dsig_avgs_v2sub,
    #                          ylabel=r'$\Delta\sigma^2$',
    #                          # kin_loc=(0.65, 0.2), star_prelim_loc=(1, 0.54, 0.78),
    #                          kin_loc=(0.58, 0.45), star_prelim_loc=(1, 0.54, 0.4),
    #                          marker_map=data_sets_markers, data_sets_bands=data_sets_bands, legend_order=legend_order,
    #                          # title=f'{cent_map[cent_plt]} Centrality, 120° Partitions'
    #                          title=f''
    #                          )  # <---

    dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], data_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors, no_hydro_label=1,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.00124, 0.00033], avgs_df=dsig_avgs_v2sub,
                             ylabel=r'$\Delta\sigma^2$', print_data=True,
                             # kin_loc=(0.65, 0.2), star_prelim_loc=(1, 0.54, 0.78),
                             kin_loc=(0.58, 0.43), star_prelim_loc=(1, 0.54, 0.37),
                             marker_map=data_sets_markers, data_sets_bands=data_sets_bands, legend_order=legend_order,
                             # title=f'{cent_map[cent_plt]} Centrality, 120° Partitions'
                             title=f''
                             )  # <---

    # plot_protons_avgs_vs_energy(dsig_avg, ['bes_def'], data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, title=f'{cent_map[cent_plt]} Centrality, {div_plt}° '
    #                                                                      f'Partitions, {samples} Samples per Event')

    dsig_avgs_v2_sub_cent8 = dsig_avgs_v2sub[dsig_avgs_v2sub['cent'] == 8]
    # plot_dvar_avgs_divs(dsig_avgs_v2_sub_cent8, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,  # <---
    #                     data_sets_labels=data_sets_labels, plot_energy_panels=True, legend_order=legend_order,
    #                     ylab=r'$\langle\Delta\sigma^2\rangle$', data_sets_bands=data_sets_bands,
    #                     plot_indiv=False, ylim=(-0.00079, 0.0001), leg_panel=5, no_hydro_label=True,
    #                     # star_prelim_loc=(1, 0.3, 0.7),
    #                     star_prelim_loc=(1, 0.5, 0.8),
    #                     # xlim=(-10, 370), title=f'0-5% Centrality, {samples} Samples per Event',
    #                     xlim=(-10, 370), title=f'',
    #                     exclude_divs=exclude_divs)
    plot_dvar_avgs_divs(dsig_avgs_v2_sub_cent8, data_sets_plt, data_sets_colors=data_sets_colors, fit=True,  # <---
                        data_sets_labels=data_sets_labels, plot_energy_panels=True, legend_order=legend_order,
                        ylab=r'$\langle\Delta\sigma^2\rangle$', data_sets_bands=data_sets_bands, print_data=True,
                        plot_indiv=False, ylim=(-0.0009, 0.0001), leg_panel=5, no_hydro_label=True,
                        # star_prelim_loc=(1, 0.3, 0.7),
                        star_prelim_loc=(1, 0.5, 0.8),
                        # xlim=(-10, 370), title=f'0-5% Centrality, {samples} Samples per Event',
                        xlim=(-10, 370), title=f'',
                        exclude_divs=exclude_divs)
    # plt.show()

    # plot_dvar_avgs_divs(dsig_avgs_v2_sub_cent8, ['bes_def'], data_sets_colors=data_sets_colors, fit=False,
    #                     data_sets_labels=data_sets_labels, plot_energy_panels=True,
    #                     ylab=r'$\langle\Delta\sigma^2\rangle$',
    #                     plot_indiv=False, ylim=(-0.00079, 0.00019), leg_panel=5,
    #                     # star_prelim_loc=(1, 0.3, 0.7),
    #                     star_prelim_loc=(1, 0.5, 0.8),
    #                     xlim=(-10, 370), title='',
    #                     # title=f'0-5% Centrality, {samples} Samples per Event',
    #                     exclude_divs=exclude_divs)

    dsig_avgs_v2_sub_cent8_div120 = dsig_avgs_v2sub[(dsig_avgs_v2sub['cent'] == 8) & (dsig_avgs_v2sub['divs'] == 120)]
    plot_protons_avgs_vs_energy(dsig_avgs_v2_sub_cent8_div120, data_sets_plt, data_sets_colors=data_sets_colors,  # <---
                                data_sets_labels=data_sets_labels, alpha=1, kin_info_loc=(0.02, 0.595),
                                star_prelim_loc=(0.25, 0.32), marker_map=data_sets_markers, ylim=(-0.00066, 0.00002),
                                data_sets_bands=data_sets_bands, legend_order=legend_order, leg_loc='lower right',
                                # title=f'{cent_map[8]} Centrality, {div_plt}° Partitions, {samples} Samples per Event')
                                title=f'{cent_map[8]} Centrality, {div_plt}° Partitions')
    # plot_protons_avgs_vs_energy(dsig_avgs_v2_sub_cent8_div120, ['bes_def'], data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, alpha=1, kin_info_loc=(0.123, 0.68),
    #                             star_prelim_loc=(0.65, 0.9), leg_loc='lower right',
    #                             title=f'{cent_map[8]} Centrality, {div_plt}° Partitions, {samples} Samples per Event')

    # dsig_avgs_v2_sub_div120 = dsig_avgs_v2sub[dsig_avgs_v2sub['divs'] == 120]
    dsig_avgs_v2_sub_div120 = dsig_avgs_v2sub[(dsig_avgs_v2sub['divs'] == 120) & (dsig_avgs_v2sub['cent'] > -1)]
    data_sets_energies_colors = \
        {'bes_def': {7: 'red', 11: 'blue', 19: 'green', 27: 'orange', 39: 'purple', 62: 'black'},
         'ampt_new_coal_epbins1': {7: 'red', 11: 'blue', 19: 'green', 27: 'orange', 39: 'purple', 62: 'black'}}
    vs_cent_sets = list(data_sets_energies_colors.keys())
    # plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, ['bes_def'], data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,  # <---
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
    #                           kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.6, 0.5),
    #                           data_sets_energies_colors=data_sets_energies_colors)
    # plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, ['bes_def'], data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type='npart',  # <---
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
    #                           kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.6, 0.5),
    #                           data_sets_energies_colors=data_sets_energies_colors)

    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, vs_cent_sets, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.3, 0.5), marker_map=data_sets_markers,
                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands,
                              print_data=True)
    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, vs_cent_sets, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type='npart',  # <---
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', alpha=0.8, errbar_alpha=0.3,
                              kin_info_loc=(0.2, 0.1), star_prelim_loc=(0.3, 0.5), marker_map=data_sets_markers,
                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)

    data_sets_cent = ['ampt_new_coal_epbins1', 'bes_def']
    legend_order = ['7.7 GeV', '11.5 GeV', '19.6 GeV', '27 GeV', '39 GeV', '62.4 GeV', 'AMPT Fit']
    plot_dsig_avg_vs_cent_2panel(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
                                 cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
                                 # title=f'{div_plt}° Partitions, {samples} Samples per Event',
                                 title='',
                                 errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
                                 kin_info_loc=(0.45, 0.1), star_prelim_loc=(0.4, 0.5), marker_map=data_sets_markers,
                                 data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    # plot_dsig_avg_vs_cent_2panel2(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
    #                               cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
    #                               # title=f'{div_plt}° Partitions, {samples} Samples per Event',
    #                               title='',
    #                               errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
    #                               kin_info_loc=(0.45, 0.1), star_prelim_loc=(0.4, 0.5), marker_map=data_sets_markers,
    #                               data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)
    plot_dsig_avg_vs_cent_2panel62ref(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors,
                                      fit=False, cent_ref=cent_ref_df, ref_type=ref_type, legend_order=None,
                                      title='', errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8, kin_info_loc=(0.2, 0.75),
                                      star_prelim_loc=None, marker_map=data_sets_markers,
                                      data_sets_energies_colors=data_sets_energies_colors,
                                      data_sets_bands=data_sets_bands)

    cent_fits = plot_dsig_avg_vs_cent_fit(dsig_avgs_v2_sub_div120, data_sets_cent, data_sets_colors=data_sets_colors,
                                          fit=True, cent_ref=cent_ref_df, ref_type=ref_type, legend_order=None,
                                          title='', errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8, kin_info_loc=(0.2, 0.75),
                                          star_prelim_loc=None, marker_map=data_sets_markers,
                                          data_sets_energies_colors=data_sets_energies_colors,
                                          data_sets_bands=None)
    print(cent_fits)
    plt.show()

    dsig_avgs_62ref = dsig_avgs_v2_sub_div120.rename(columns={'name': 'data_set'})
    plot_div_fits_vs_cent_62res(dsig_avgs_62ref, data_sets_cent, data_sets_colors, data_sets_labels, ref_type=ref_type,
                                cent_ref=cent_ref_df, val_col='avg', err_col='avg_err')

    df_fits = pd.read_csv(f'{base_path}{df_partitions_fits_name}')
    df_fits_model = pd.read_csv(f'{base_path}{df_partitions_fits_model_name}')
    df_fits = pd.concat([df_fits, df_fits_model])

    # data_sets_energies_cmaps = dict(zip(data_sets_cent, ['winter', 'copper']))
    # data_sets_markers2 = dict(zip(data_sets_cent, ['s', 'o']))
    # plot_div_fits_vs_cent(df_fits, data_sets_cent,  # data_sets_energies_cmaps=data_sets_energies_cmaps,
    #                       data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
    #                       ref_type=ref_type, data_sets_colors=data_sets_energies_colors,
    #                       data_sets_markers=data_sets_markers2)

    plot_div_fits_vs_cent_2panel(df_fits, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
                                 cent_ref=cent_ref_df, ref_type=ref_type,  # <---
                                 # title=f'{div_plt}° Partitions, {samples} Samples per Event',
                                 title='', print_data=True,
                                 errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
                                 kin_info_loc=(0.45, 0.03), marker_map=data_sets_markers,
                                 data_sets_energies_colors=data_sets_energies_colors,
                                 # data_sets_bands=data_sets_bands
                                 )
    plt.show()

    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels, data_sets=data_sets_plt)
    # plot_slope_div_fits_simpars(df_fits)

    # Plot avg dsig2 vs refmult for mixed events. Wierd stuff at most peripheral bin or two

    # dvar_vs_protons(df_dsigma_with_v2sub, div_plt, 0, [7], ['raw', 'mix', 'diff', 'v2_sub'], ['bes_def'],
    #                 plot=True, avg=True, alpha=1.0, ylabel=r'$\Delta\sigma^2$',
    #                 data_sets_labels={'bes_def': {'raw': 'Uncorrected', 'diff': 'Mixed Corrected', 'mix': 'Mixed',
    #                                               'v2_sub': 'Mixed and Flow Corrected'}},
    #                 data_sets_colors={'bes_def': {'raw': 'blue', 'diff': 'red', 'v2_sub': 'black', 'mix': 'purple'}},
    #                 marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^', 'v2_sub': '*'}},
    #                 # y_ranges=[-0.00124, 0.0009])  # v2 sub demo
    #                 y_ranges=[-0.017, 0.007], kin_info_loc=(0.26, 0.94))  # v2 sub demo
    #
    # dvar_vs_protons(df_dsigma_with_v2sub, div_plt, 2, [7], ['raw', 'mix', 'diff', 'v2_sub'], ['bes_def'],
    #                 plot=True, avg=True, alpha=1.0, ylabel=r'$\Delta\sigma^2$',
    #                 data_sets_labels={'bes_def': {'raw': 'Uncorrected', 'diff': 'Mixed Corrected', 'mix': 'Mixed',
    #                                               'v2_sub': 'Mixed and Flow Corrected'}},
    #                 data_sets_colors={'bes_def': {'raw': 'blue', 'diff': 'red', 'v2_sub': 'black', 'mix': 'purple'}},
    #                 marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^', 'v2_sub': '*'}},
    #                 # y_ranges=[-0.00124, 0.0009])  # v2 sub demo
    #                 y_ranges=[-0.017, 0.007], kin_info_loc=(0.26, 0.94))  # v2 sub demo
    #
    # df_mix = []
    # for cent in range(0, 9):
    #     df_mix_cent = dvar_vs_protons(df_dsigma, div_plt, cent, [7, 11, 19, 27, 39, 62],
    #                                   ['mix'], data_sets_plt, plot=False, avg=True)
    #     print(df_mix_cent)
    #     df_mix.append(df_mix_cent)
    # print(df_mix)
    # df_mix = pd.concat(df_mix)
    # print(df_mix)
    # print(df_mix.columns)
    # df_mix = df_mix[(df_mix['divs'] == 120) & (df_mix['data_type'] == 'mix')]
    # plot_dsig_avg_vs_cent_2panel(df_mix, data_sets_cent, data_sets_colors=data_sets_colors, fit=False,
    #                              cent_ref=cent_ref_df, ref_type=ref_type, legend_order=legend_order,  # <---
    #                              # title=f'{div_plt}° Partitions, {samples} Samples per Event',
    #                              title='',
    #                              errbar_alpha=0.3, xlim=(-20, 720), alpha=0.8,
    #                              kin_info_loc=(0.45, 0.1), star_prelim_loc=(0.4, 0.5), marker_map=data_sets_markers,
    #                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)

    # Save all open figures
    for i in plt.get_fignums():
        plt.figure(i)
        window_title = plt.gcf().canvas.manager.get_window_title().replace(' ', '_')
        plt.savefig(f'{save_dir}{window_title}_{i}.png')
        plt.savefig(f'{save_dir}{window_title}_{i}.pdf')

    plt.show()


def plot_method_paper_figs():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144

    presentation_mode = False
    if presentation_mode:
        # plt.rcParams['axes.labelsize'] = 14  # Adjust the value as needed
        # plt.rcParams['axes.titlesize'] = 16  # Adjust the value as needed
        # plt.rcParams['legend.fontsize'] = 14  # Adjust the value as needed
        plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.2

    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'C:/Users/Dyn04/Research/'
    df_model_name = 'Bes_with_Sys/binom_slice_vars_model.csv'
    df_dsigma_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma.csv'
    df_dsigma_v2sub_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma_v2sub.csv'
    # df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_raw_model.csv'
    # df_partitions_fits_model_name = 'Bes_with_Sys/partition_width_fits_model.csv'
    df_partitions_fits_model_name = 'Bes_with_Sys/partition_width_fits_raw_model.csv'
    df_partitions_fits_model_nov2sub_name = 'Bes_with_Sys/partition_width_fits_nov2sub_model.csv'

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [89, 356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    data_type_plt = 'diff'  # 'raw' 'diff'
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))
    data_sets_markers = dict(zip(data_sets_plt, [dict(zip(['raw', 'mix', 'diff'], [x, x, x]))
                                                 for x in ['o', 's', '^', '*']]))
    data_sets_bands = []
    legend_order = ['MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$', 'AMPT']

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    df = pd.read_csv(f'{base_path}{df_model_name}')
    df = df[df['stat'] == stat_plot]

    df_dsigma = pd.read_csv(f'{base_path}{df_dsigma_model_name}')

    plot_method_paper_event()

    # plt.rcParams.update({'font.size': 12})
    method_paper_plot()  # From Plotter/presentation_plots
    # cluster_voids_example_plots()  # From Anti_Clustering/pdf_examples. Take time to generate.

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex='all', figsize=(8, 4), dpi=144)
    stat_binom_vs_protons(df, stat_plot, div_plt, cent_plt, 39, ['raw', 'mix'], 'ampt_new_coal_epbins1',
                          data_sets_labels=data_sets_labels, ax_in=ax1)
    # ax1.set_title(f'AMPT 39 GeV, 0-5% Centrality, 120° Partitions, 72 Samples/Event', pad=-20)
    ax1.text(0.5, 0.9, 'AMPT 39 GeV, 0-5% Centrality, 120° Partitions', ha='center',
             va='center', transform=ax1.transAxes, fontsize=16)
    ax1.set_ylim(top=ax1.get_ylim()[1] * 1.1)
    ax1.set_ylabel(r'$\sigma^2$', fontsize=15)
    ax1.text(0.005, 0.9, '(a)', ha='left', va='center', transform=ax1.transAxes, fontsize=12)
    ax1.legend(loc='lower right')

    markers = {'ampt_new_coal_epbins1': {'raw': 'o', 'mix': 's', 'diff': '^'}}
    dvar_vs_protons(df_dsigma, div_plt, cent_plt, [39], ['raw', 'mix'], ['ampt_new_coal_epbins1'],
                    data_sets_colors=None, data_sets_labels=data_sets_labels, ax_in=ax2, title='', marker_map=markers,
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00124, 0.0009], ylabel=r'$\Delta\sigma^2$',
                    kin_info_loc=(0.35, 0.35), leg=False)
    ax2.set_xlabel('Total Protons in Event', fontsize=15)
    ax2.set_ylabel(r'$\Delta\sigma^2$', fontsize=15)
    ax2.text(0.005, 0.9, '(b)', ha='left', va='center', transform=ax2.transAxes, fontsize=12)
    fig.subplots_adjust(hspace=0.0, left=0.126, right=0.995, top=0.995, bottom=0.115)
    fig.canvas.manager.set_window_title('variance_dsig2_vs_tprotons_ampt_39gev')

    df_dsigma_v2sub = pd.read_csv(f'{base_path}{df_dsigma_v2sub_model_name}')

    dsig_avgs_v2sub = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_model_name}')

    dsig_avgs_v2sub['data_type'] = data_type_plt
    subplot_adjust = {'wspace': 0.0, 'hspace': 0.0, 'left': 0.075, 'right': 0.995, 'top': 0.995, 'bottom': 0.09}
    dvar_vs_protons_energies(df_dsigma_v2sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], [data_type_plt], data_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors, no_hydro_label=1,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.00053], avgs_df=dsig_avgs_v2sub,
                             ylabel=r'$\Delta\sigma^2$', kin_loc=(0.65, 0.41), legend_panel=3,
                             marker_map=data_sets_markers, data_sets_bands=data_sets_bands, legend_order=legend_order,
                             title=f'', alpha=0.8, fig_splt_adjust=subplot_adjust, plot_letters=True)

    dsig_avgs_v2_sub_cent8 = dsig_avgs_v2sub[dsig_avgs_v2sub['cent'] == 8]
    subplot_adjust = {'wspace': 0.0, 'hspace': 0.0, 'left': 0.075, 'right': 0.995, 'top': 0.995, 'bottom': 0.095}
    plot_dvar_avgs_divs(dsig_avgs_v2_sub_cent8, data_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                        data_sets_labels=data_sets_labels, plot_energy_panels=True, legend_order=legend_order,
                        ylab=r'$\langle\Delta\sigma^2\rangle$', data_sets_bands=data_sets_bands,
                        plot_indiv=False, ylim=(-0.00054, 0.00009), leg_panel=5, no_hydro_label=True,
                        xlim=(-10, 370), title=f'', leg_frameon=True, exclude_divs=exclude_divs, alpha=0.8,
                        data_sets_markers=data_sets_markers, fig_splt_adj=subplot_adjust, panel_letters=True)

    plt.rcParams["figure.figsize"] = (7, 3.5)
    subplot_adjust = {'left': 0.142, 'right': 0.995, 'top': 0.995, 'bottom': 0.140}
    dsig_avgs_v2_sub_cent8_div120 = dsig_avgs_v2sub[(dsig_avgs_v2sub['cent'] == 8) & (dsig_avgs_v2sub['divs'] == 120)]
    plot_protons_avgs_vs_energy(dsig_avgs_v2_sub_cent8_div120, data_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, alpha=0.6, kin_info_loc=(0.02, 0.45),
                                marker_map=data_sets_markers, fig_splt_adj=subplot_adjust,
                                data_sets_bands=data_sets_bands, legend_order=legend_order, title=f'')

    dsig_avgs_v2_sub_div120 = dsig_avgs_v2sub[dsig_avgs_v2sub['divs'] == 120]
    data_sets_energies_colors = \
        {'ampt_new_coal_epbins1': {7: 'red', 11: 'blue', 19: 'green', 27: 'orange', 39: 'purple', 62: 'black'}}
    data_sets_energies_markers = {'ampt_new_coal_epbins1': {7: 'o', 11: 's', 19: '^', 27: 'p', 39: 'D', 62: '*'}}
    vs_cent_sets = list(data_sets_energies_colors.keys())

    splt_adjust = {'left': 0.142, 'right': 0.995, 'top': 0.995, 'bottom': 0.11}
    plot_protons_avgs_vs_cent(dsig_avgs_v2_sub_div120, vs_cent_sets, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'', alpha=0.8, errbar_alpha=0.3, xerr=False, xlim=(-10, 699), figsize=(7, 4.5),
                              kin_info_loc=(0.2, 0.1), marker_map=data_sets_energies_markers, fig_splt_adj=splt_adjust,
                              data_sets_energies_colors=data_sets_energies_colors, data_sets_bands=data_sets_bands)

    df_fits = pd.read_csv(f'{base_path}{df_partitions_fits_model_name}')
    data_sets_cent = ['ampt_new_coal_epbins1']

    # plot_div_fits_vs_cent(df_fits, data_sets_cent,  # data_sets_energies_cmaps=data_sets_energies_cmaps,
    #                       data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df, ls='-',
    #                       ref_type=ref_type, data_sets_colors=data_sets_energies_colors,
    #                       data_sets_markers=data_sets_energies_markers)

    plot_div_fits_vs_cent_single_set(df_fits, data_sets_cent, data_sets_labels=data_sets_labels, title=None, ls='-',
                                     data_sets_colors=data_sets_energies_colors, cent_ref=cent_ref_df, show_xerr=False,
                                     data_sets_markers=data_sets_energies_markers, ref_type=ref_type)

    df_fits_nov2sub = pd.read_csv(f'{base_path}{df_partitions_fits_model_nov2sub_name}')
    plot_div_fits_vs_cent_single_set(df_fits_nov2sub, data_sets_cent, data_sets_labels=data_sets_labels, title=None,
                                     data_sets_colors=data_sets_energies_colors, cent_ref=cent_ref_df, show_xerr=False,
                                     data_sets_markers=data_sets_energies_markers, ref_type=ref_type, ls='-')

    df_fits_amptv2rp = df_fits[df_fits['data_set'] == 'ampt_new_coal_epbins1_v2rp']
    df_fits_amptv2rp['data_set'] = 'ampt_new_coal_epbins1'
    plot_div_fits_vs_cent_single_set(df_fits_amptv2rp, data_sets_cent, data_sets_labels=data_sets_labels, title=None,
                                     data_sets_colors=data_sets_energies_colors, cent_ref=cent_ref_df, show_xerr=False,
                                     data_sets_markers=data_sets_energies_markers, ref_type=ref_type, ls='-')

    df_fits_amptv2cum = df_fits[df_fits['data_set'] == 'ampt_new_coal_epbins1_v2cum']
    df_fits_amptv2cum['data_set'] = 'ampt_new_coal_epbins1'
    plot_div_fits_vs_cent_single_set(df_fits_amptv2cum, data_sets_cent, data_sets_labels=data_sets_labels, title=None,
                                     data_sets_colors=data_sets_energies_colors, cent_ref=cent_ref_df, show_xerr=False,
                                     data_sets_markers=data_sets_energies_markers, ref_type=ref_type, ls='-')

    plt.show()

    # Gaussian correlation model simulations
    sim_base_path = f'{base_path}Binomial_Slice_Moments/'
    df_sim_dsig = pd.read_csv(f'{sim_base_path}binom_slice_stats_sim_demos_dsigma.csv')
    df_sim_dsig_avgs = pd.read_csv(f'{sim_base_path}binom_slice_stats_sim_demos_dsigma_avg.csv')
    df_sim_width_fits = pd.read_csv(f'{sim_base_path}binom_slice_stats_sim_demos_width_fits.csv')

    amps = ['002', '006', '008', '01']
    spreads = ['04', '05', '06', '07', '08', '09', '1', '11', '12']
    amp_spreads_all = [amps, spreads]

    amps = ['002', '006', '01']
    spreads = ['04', '05', '06', '07', '08', '09', '1', '11', '12']
    amp_spreads_z_w = [amps, spreads]

    amps = ['002', '006', '01']
    spreads = ['1']
    amp_spreads_dsig_n = [amps, spreads]

    amps = ['002', '008']
    spreads = ['04', '08']
    amp_spreads_avg_w = [amps, spreads]

    sim_set_defs = [amp_spreads_all, amp_spreads_dsig_n, amp_spreads_avg_w, amp_spreads_z_w]
    sim_sets_out = []
    for (amps, spreads) in sim_set_defs:
        sim_sets = []
        for amp in amps:
            for spread in spreads:
                sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
                sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
        sim_sets = sorted(sim_sets, reverse=True)
        sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])
        sim_sets_out.append(sim_sets)
    sim_sets_all, sim_sets_dsig_n, sim_sets_avg_w, sim_sets_z_w = sim_sets_out

    amp_spread_markers = {0.002: {0.4: 'o', 0.8: 's'}, 0.006: {0.4: '^', 0.8: 'v'}, 0.008: {0.4: '^', 0.8: 'v'},
                          0.01: {0.4: 'D', 0.8: 'x'}}
    amp_spread_colors = {0.002: {0.4: 'black', 0.8: 'red'}, 0.006: {0.4: 'green', 0.8: 'orange'},
                         0.008: {0.4: 'black', 0.8: 'red'}, 0.01: {0.4: 'purple', 0.8: 'black'}}
    amp_spread_ls = {0.002: {0.4: '-', 0.8: '--'}, 0.006: {0.4: '-', 0.8: '--'}, 0.008: {0.4: '-', 0.8: '--'},
                     0.01: {0.4: '-', 0.8: '--'}}
    data_sets_labels_sim, data_sets_colors_sim, data_sets_fill_sim, data_sets_markers_sim = {}, {}, {}, {}
    data_sets_ls_sim = {}
    for sim_set in sim_sets_all:
        label, fillstyle, marker = '', '', ''
        amp, spread = get_name_amp_spread(sim_set)
        if '_clmul_' in sim_set:
            label += 'Attractive ' + rf'$A=+{amp:.3f}$ $Σ={spread}$'
            fillstyle = 'none'
        elif '_aclmul_' in sim_set:
            label += 'Repulsive ' + rf'$A=-{amp:.3f}$ $Σ={spread}$'
            fillstyle = 'full'
        data_sets_labels_sim.update({sim_set: label})
        data_sets_fill_sim.update({sim_set: fillstyle})
        if amp in amp_spread_markers and spread in amp_spread_markers[amp]:
            marker = amp_spread_markers[amp][spread]
        else:
            marker = 'o'
        if amp in amp_spread_colors and spread in amp_spread_colors[amp]:
            color = amp_spread_colors[amp][spread]
        else:
            color = 'black'
        if amp in amp_spread_ls and spread in amp_spread_ls[amp]:
            ls = amp_spread_ls[amp][spread]
        else:
            ls = '-'
        data_sets_markers_sim.update({sim_set: marker})
        data_sets_colors_sim.update({sim_set: color})
        data_sets_ls_sim.update({sim_set: ls})

    fig, ax = plt.subplots(figsize=(8, 4), dpi=144)
    dvar_vs_protons(df_sim_dsig, div_plt, cent_plt, ['sim'], ['raw'], sim_sets_dsig_n, plot=True, avg=True,
                    data_sets_labels=data_sets_labels_sim, ylabel=r'$\Delta\sigma^2$', data_sets_bands=sim_sets_dsig_n,
                    title=None, ax_in=ax, ylim=(-0.00099, 0.00124), kin_info_loc=(0.4, 0.17), leg=False)
    title = f'Gaussian Correlation Model - {div_plt}° Partitions'
    ax.text(0.5, 0.98, title, ha='center', va='top', transform=ax.transAxes, fontsize=16)
    ax.set_xlabel('Total Protons in Event', fontsize=14)
    ax.set_ylabel(r'$\Delta\sigma^2$', fontsize=14)
    ax.axhline(0, color='black', zorder=1)

    handles, labels = ax.get_legend_handles_labels()
    attractive_handles = [handle for handle, label in zip(handles, labels) if 'Attractive' in label]
    repulsive_handles = [handle for handle, label in zip(handles, labels) if 'Repulsive' in label]
    attractive_label = ax.legend(attractive_handles,
                                 [label.replace('Attractive ', '') for label in labels if 'Attractive' in label],
                                 loc='upper right', title='Attractive', bbox_to_anchor=(0.99, 0.935), title_fontsize=11)
    repulsive_label = ax.legend(repulsive_handles,
                                [label.replace('Repulsive ', '') for label in labels if 'Repulsive' in label],
                                loc='upper right', title='Repulsive', bbox_to_anchor=(0.99, 0.27), title_fontsize=11)
    ax.add_artist(attractive_label)

    fig.subplots_adjust(left=0.135, right=0.995, top=0.995, bottom=0.115)
    fig.canvas.manager.set_window_title('dsigma_vs_Total_Protons_simGeV_120_example')

    plt.rcParams["figure.figsize"] = (9, 4.5)
    # plt.rcParams['figure.subplot.left'], plt.rcParams['figure.subplot.right'] = 0.125, 0.995
    # plt.rcParams['figure.subplot.bottom'], plt.rcParams['figure.subplot.top'] = 0.1, 0.94
    subplot_adjust = {'left': 0.111, 'right': 0.995, 'top': 0.995, 'bottom': 0.105}
    plot_dvar_avgs_divs(df_sim_dsig_avgs, sim_sets_avg_w, data_sets_colors=data_sets_colors_sim, fit=True,
                        data_sets_labels=data_sets_labels_sim, plot_energy_panels=False, ylim=(-0.00045, 0.00045),
                        data_sets_markers=data_sets_markers_sim, data_sets_ls=data_sets_ls_sim, xlim=(-15, 485),
                        ylab=r'$\langle\Delta\sigma^2\rangle$', data_sets_fills=data_sets_fill_sim, alpha=0.8,
                        exclude_divs=exclude_divs, fig_splt_adj=subplot_adjust, title='', leg=False)

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    attractive_handles = [handle for handle, label in zip(handles, labels) if 'Attractive' in label]
    repulsive_handles = [handle for handle, label in zip(handles, labels) if 'Repulsive' in label]
    attractive_label = ax.legend(attractive_handles,
                                 [label.replace('Attractive ', '') for label in labels if 'Attractive' in label],
                                 loc='center right', title='Attractive', bbox_to_anchor=(1.0, 0.75), title_fontsize=12)
    repulsive_label = ax.legend(repulsive_handles,
                                [label.replace('Repulsive ', '') for label in labels if 'Repulsive' in label],
                                loc='center right', title='Repulsive', bbox_to_anchor=(1.0, 0.25), title_fontsize=12)
    ax.add_artist(attractive_label)

    plot_b_vs_amp(df_sim_width_fits, alpha=0.8, ylim=(-0.00046, 0.00053))

    amp_shifts = {0.002: -0.01, 0.006: 0.0, 0.01: +0.01}
    amp_colors = {0.002: 'red', 0.006: 'black', 0.01: 'blue'}
    amp_markers = {0.002: 's', 0.006: 'o', 0.01: '^'}
    plot_z_vs_spread(df_sim_width_fits, amps=list(amp_colors.keys()), amps_colors=amp_colors, amps_markers=amp_markers,
                     amps_x_shifts=amp_shifts, alpha=0.7)

    # plot_b_vs_amp_sig_dep(df_sim_width_fits, alpha=0.8)

    # plot_slope_div_fits_simpars(df_sim_width_fits)
    plot_full_test_from_file()  # Momentum_Conservation_Sim.momentum_conservation_model.py
    # Get current figure and set window title
    fig = plt.gcf()
    fig.canvas.manager.set_window_title('dsig2_vs_total_particles_fit')

    plt.show()


def plot_star_model():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'cf_tprotons_fits.csv'
    df_partitions_fits_name = 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST EV $1fm^3$']))

    # data_sets_plt = ['cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['MUSIC+FIST', 'MUSIC+FIST EV']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_star_model_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
                     'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    v2_cf_in_dir = 'F:/Research/Data_CF/default_resample_epbins1/CF_rapid05_resample_norotate_epbins1_0/'
    v2_cfev_in_dir = 'F:/Research/Data_CFEV/default_resample_epbins1/CFEV_rapid05_resample_norotate_epbins1_0/'
    v2_cfevb342_in_dir = 'F:/Research/Data_CFEVb342/default_resample_epbins1/' \
                         'CFEVb342_rapid05_resample_norotate_epbins1_0/'
    df_name = 'Binomial_Slice_Moments/binom_slice_vars.csv'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_avgs_name = 'dsig_tprotons_avgs_cent8.csv'
    df_partitions_fits_name = 'partitions_fits_cent8.csv'
    df_path = base_path + df_name
    sim_sets = []

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))

    data_sets_plt = ['bes_def', 'ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))

    # data_sets_plt = ['bes_resample_epbins1', 'ampt_new_coal_epbins1']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}
    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}
    v2_cf_vals = {2: read_flow_values(v2_cf_in_dir)}
    v2_cfev_vals = {2: read_flow_values(v2_cfev_in_dir)}
    v2_cfevb342_vals = {2: read_flow_values(v2_cfevb342_in_dir)}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # df = df[df['name'].str.contains('bes')]

    df = df[df['stat'] == stat_plot]
    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])

    stat_binom_vs_protons(df, stat_plot, div_plt, cent_plt, 62, ['raw', 'mix'], 'ampt_new_coal_epbins1',
                          data_sets_labels=data_sets_labels)

    raw_to_mix_stat_err(pd.concat([df_raw, df_mix], ignore_index=True), div_plt, 4, 11, 'bes_def')

    # plt.show()

    dvar_vs_protons(pd.concat([df_raw, df_mix]), div_plt, cent_plt, [62], ['raw', 'mix'], ['ampt_new_coal_epbins1'],
                    plot=True, avg=False)

    dvar_vs_protons(df_diff, div_plt, cent_plt, [62], ['diff'], ['ampt_new_coal_epbins1'], plot=True, avg=False,
                    data_sets_labels=data_sets_labels, data_sets_colors=data_sets_colors)

    dvar_vs_protons(df_diff, div_plt, cent_plt, [62], ['diff'], all_sets_plt, plot=True, avg=False,
                    data_sets_labels=data_sets_labels, data_sets_colors=data_sets_colors)

    dvar_vs_protons(pd.concat([df_raw, df_mix, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)

    dvar_vs_protons_energies(df_diff, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], all_sets_plt,
                             plot=True, avg=False, plot_avg=False, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.00054])

    dsig_avg = dvar_vs_protons_energies(df_diff, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['diff'], all_sets_plt,
                                        plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                                        data_sets_labels=data_sets_labels, y_ranges=[-0.00099, 0.00054])

    plot_protons_avgs_vs_energy(dsig_avg, all_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, title=f'{cent_map[cent_plt]} Centrality, {div_plt}° '
                                                                         f'Partitions, {samples} Samples per Event')

    plt.show()

    dsig_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsig_avgs_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets_plt, plot=False,
                                            avg=True)
        dsig_avgs_div_raw.loc[:, 'name'] = dsig_avgs_div_raw['name'] + '_raw'
        dsig_avgs.append(dsig_avgs_div_raw)

        dsig_avgs_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets_plt, plot=False,
                                            avg=True)
        dsig_avgs_div_mix.loc[:, 'name'] = dsig_avgs_div_mix['name'] + '_mix'
        dsig_avgs.append(dsig_avgs_div_mix)

        dsig_avgs_div_diff = dvar_vs_protons(df_diff, div, cent_plt, energies_fit, ['diff'], all_sets_plt, plot=False,
                                             avg=True)
        dsig_avgs_div_diff.loc[:, 'name'] = dsig_avgs_div_diff['name'] + '_sub'

        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'ampt_new_coal_epbins1_sub',
                                                  'ampt_new_coal_epbins1', v2_ampt_vals, div, cent_plt)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'bes_resample_epbins1_sub',
                                                  'bes_resample_epbins1', v2_star_vals, div, cent_plt)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cf_resample_epbins1_sub',
                                                  'cf_resample_epbins1', v2_cf_vals, div, cent_plt)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cfev_resample_epbins1_sub',
                                                  'cfev_resample_epbins1', v2_cfev_vals, div, cent_plt)
        dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cfevb342_resample_epbins1_sub',
                                                  'cfevb342_resample_epbins1', v2_cfevb342_vals, div, cent_plt)
        dsig_avgs.append(dsig_avgs_div_diff)

    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    if df_tproton_avgs_name is not None:
        dsig_avgs.to_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}', index=False)

    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_raw', '_mix', '_sub']]
        colors = dict(zip(data_sets, ['blue', 'green', 'red']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Raw', ' Mix', ' Sub']]))
        plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels,
                            ylab=r'$\langle\Delta\sigma^2\rangle$')
    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_sub', '']]
        colors = dict(zip(data_sets, ['blue', 'red']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Original', ' v2 Corrected']]))
        plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels)
    plot_dvar_avgs_divs(dsig_avgs, data_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                        data_sets_labels=data_sets_labels)

    df_fits = plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                  data_sets_labels=data_sets_labels, plot_energy_panels=False)
    if df_partitions_fits_name is not None:
        df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_star_var_sys():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
    #                  'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    v2_star_in_dir = 'F:/Research/Data/default/' \
                     'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/'
    sys_dir = 'F:/Research/Data/default_sys/'
    sys_default_dir = 'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/'
    df_name = 'Binomial_Slice_Moments/binom_slice_vars_bes_sys.csv'

    plot = True
    # plot = False
    calc_finals = False
    # calc_finals = True
    threads = 11
    sys_pdf_out_path = f'{base_path}systematic_plots.pdf'
    indiv_pdf_out_path = f'F:/Research/Results/BES_QA_Plots/Systematics/'
    df_def_out_name = 'Bes_with_Sys/binom_slice_vars_bes.csv'
    # df_def_out_name = None
    df_def_dsigma_out_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma.csv'
    # df_def_dsigma_out_name = None
    df_def_dsigma_v2sub_out_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma_v2sub.csv'
    df_def_avgs_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_bes.csv'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    fits_out_base = 'Base_Zero_Fits'
    df_partitions_fits_name = 'Bes_with_Sys/partition_width_fits_bes.csv'
    df_path = base_path + df_name

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 180
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    # divs_all = [60, 120, 180]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [7, 11]
    samples = 72  # For title purposes only

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    # sys_info_dict = {  # Old systematics (in my opinion, better) used in thesis. For all gaussian.
    #     'vz': {'name': 'vz range', 'title': 'vz', 'decimal': None, 'default': None,
    #            'sys_vars': ['low7', 'high-7', 'low-5_vzhigh5'], 'val_unit': ' cm',
    #            'sys_var_order': ['low7', 'low-5_vzhigh5', 'high-7']},
    #     'Efficiency': {'name': 'efficiency', 'title': 'efficiency', 'decimal': 2, 'default': 0,
    #                    'sys_vars': [95.0, 90.0], 'val_unit': '%', 'sys_var_order': [95.0, 90.0, 85.0, 80.0]},
    #     'dca': {'name': 'dca', 'title': 'dca', 'decimal': 1, 'default': 1, 'sys_vars': [0.8, 1.2], 'val_unit': ' cm',
    #             'sys_var_order': [0.5, 0.8, 1.2, 1.5]},
    #     'nsprx': {'name': r'n$\sigma$ proton', 'title': r'n$\sigma$ proton', 'decimal': 1, 'default': 1,
    #               'sys_vars': [0.9, 1.1], 'val_unit': '', 'sys_var_order': [0.75, 0.9, 1.1, 1.25]},
    #     'm2r': {'name': r'$m^2$ range', 'title': 'm2 range', 'decimal': 0, 'default': 0.6, 'sys_vars': [0.4, 0.8],
    #             'val_unit': ' GeV', 'sys_var_order': [0.2, 0.4, 0.8, 1.0]},
    #     'nhfit': {'name': 'nHits fit', 'title': 'nhits fit', 'decimal': 2, 'default': 20, 'sys_vars': [15, 25],
    #               'val_unit': '', 'sys_var_order': [15, 25]},
    #     'sysrefshift': {'name': 'refmult3 shift', 'title': 'ref3 shift', 'decimal': None, 'default': 0,
    #                     'sys_vars': ['-1', '1'], 'val_unit': '', 'sys_var_order': ['-1', '1']},
    #     'dcxyqa': {'name': 'dcaxy qa', 'title': 'dcaxy qa', 'decimal': None, 'default': None,
    #                'sys_vars': ['tight', 'loose'], 'val_unit': '',
    #                'sys_var_order': ['2tight', 'tight', 'loose', '2loose']},
    #     'pileupqa': {'name': 'pile-up qa', 'title': 'pile-up qa', 'decimal': None, 'default': None,
    #                  'sys_vars': ['tight', 'loose'], 'val_unit': '',
    #                  'sys_var_order': ['2tight', 'tight', 'loose', '2loose']},
    #     'mix_rand_': {'name': 'mix rand', 'title': 'mix rand', 'decimal': 1, 'default': 0, 'sys_vars': None,
    #                   'val_unit': '', 'sys_var_order': None},
    #     'all_rand_': {'name': 'all rand', 'title': 'all rand', 'decimal': 1, 'default': 0, 'sys_vars': None,
    #                   'val_unit': '', 'sys_var_order': None},
    # }

    sys_info_dict = {
        'vz': {'name': 'vz range', 'title': 'vz', 'decimal': None, 'default': None,
               'sys_vars': ['low7', 'high-7'], 'val_unit': ' cm',
               'sys_var_order': ['low7', 'low-5_vzhigh5', 'high-7'], 'prior': 'flat_two_side'},
        'Efficiency': {'name': 'efficiency', 'title': 'efficiency', 'decimal': 2, 'default': 0,
                       'sys_vars': [90.0], 'val_unit': '%', 'sys_var_order': [95.0, 90.0, 85.0, 80.0],
                       'prior': 'flat_one_side'},
        'dca': {'name': 'dca', 'title': 'dca', 'decimal': 1, 'default': 1, 'sys_vars': [0.5], 'val_unit': ' cm',
                'sys_var_order': [0.5, 0.8, 1.2, 1.5], 'prior': 'flat_one_side'},
        'nsprx': {'name': r'n$\sigma$ proton', 'title': r'n$\sigma$ proton', 'decimal': 1, 'default': 1,
                  'sys_vars': [0.75], 'val_unit': '', 'sys_var_order': [0.75, 0.9, 1.1, 1.25],
                  'prior': 'flat_one_side'},
        'm2r': {'name': r'$m^2$ range', 'title': 'm2 range', 'decimal': 0, 'default': 0.6, 'sys_vars': [0.4],
                'val_unit': ' GeV', 'sys_var_order': [0.2, 0.4, 0.8, 1.0], 'prior': 'flat_one_side'},
        'nhfit': {'name': 'nHits fit', 'title': 'nhits fit', 'decimal': 2, 'default': 20, 'sys_vars': [25],
                  'val_unit': '', 'sys_var_order': [15, 25], 'prior': 'flat_one_side'},
        'sysrefshift': {'name': 'refmult3 shift', 'title': 'ref3 shift', 'decimal': None, 'default': 0,
                        'sys_vars': ['-1', '1'], 'val_unit': '', 'sys_var_order': ['-1', '1'],
                        'prior': 'flat_two_side'},
        'dcxyqa': {'name': 'dcaxy qa', 'title': 'dcaxy qa', 'decimal': None, 'default': None,
                   'sys_vars': ['tight'], 'val_unit': '', 'sys_var_order': ['2tight', 'tight', 'loose', '2loose'],
                   'prior': 'flat_one_side'},
        'pileupqa': {'name': 'pile-up qa', 'title': 'pile-up qa', 'decimal': None, 'default': None,
                     'sys_vars': ['tight'], 'val_unit': '', 'sys_var_order': ['2tight', 'tight', 'loose', '2loose'],
                     'prior': 'flat_one_side'},
        'mix_rand_': {'name': 'mix rand', 'title': 'mix rand', 'decimal': 1, 'default': 0, 'sys_vars': None,
                      'val_unit': '', 'sys_var_order': None, 'prior': None},
        'all_rand_': {'name': 'all rand', 'title': 'all rand', 'decimal': 1, 'default': 0, 'sys_vars': None,
                      'val_unit': '', 'sys_var_order': None, 'prior': None},
    }

    sys_include_sets = sys_info_dict_to_var_names(sys_info_dict)
    sys_priors = sys_info_dict_to_priors(sys_info_dict)
    sys_include_names = [y for x in sys_include_sets.values() for y in x]

    data_sets_plt = ['bes_def', 'dca08', 'nsprx09', 'm2r4', 'nhfit25']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'green', 'black', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['bes_def', 'dca08', 'nsprx09', 'm2r4', 'nhfit25']))

    df = pd.read_csv(df_path)
    df = df.dropna()
    df['name'] = df['name'].str.replace('bes_sys_', '')
    all_sets = pd.unique(df['name'])
    print(all_sets)
    print(sys_include_names)

    # For all sets where 'ampt' is not in the name, make a copy of all 'bes_def' rows in cent_ref_df for this data_set
    cent_ref_extra_sets = []
    for data_set_name in all_sets:
        print(data_set_name)
        if 'ampt' not in data_set_name and data_set_name != 'bes_def':
            cent_ref_extra_sets.append(cent_ref_df[cent_ref_df['data_set'] == 'bes_def'].assign(data_set=data_set_name))
    cent_ref_extra_sets = pd.concat(cent_ref_extra_sets)
    cent_ref_df = pd.concat([cent_ref_df, cent_ref_extra_sets])
    print(cent_ref_df)
    print(pd.unique(cent_ref_df['data_set']))
    for data_set_name in all_sets:
        print(cent_ref_df[cent_ref_df['data_set'] == data_set_name])

    rand_sets = [set_name for set_name in all_sets if '_rand' in set_name]
    non_rand_sets = [set_name for set_name in all_sets if '_rand' not in set_name]

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}
    v2_sys_vals = {name: {2: read_flow_values(get_set_dir(name, sys_default_dir, sys_dir))} for name in all_sets
                   if name != 'bes_def'}
    v2_sys_vals.update({'bes_def': v2_star_vals})

    df = df[df['stat'] == stat_plot]

    # Get k2 raw, mix, diff systematics
    if df_def_out_name is not None and calc_finals:
        df_def_with_sys = get_sys(df, 'bes_def', sys_include_sets, sys_priors,
                                  group_cols=['divs', 'energy', 'cent', 'data_type', 'total_protons'])
        df_def_with_sys.to_csv(f'{base_path}{df_def_out_name}', index=False)

    # Calculate dsigma with k2 values and get systematics
    df = df[df['stat'] == 'k2']
    df = df.drop('stat', axis=1)
    print('Calc dsigma')
    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])
    df_dsigma_types = pd.concat([df_raw, df_mix, df_diff])
    print('Calc diff nlo error')
    df_dsigma_types = add_diff_nlo_err(df_dsigma_types, group_cols=['energy', 'cent', 'name', 'total_protons'],
                                       exclude_divs=[356, 89])

    if df_def_dsigma_out_name is not None and calc_finals:
        df_def_dsigma = get_sys(df_dsigma_types, 'bes_def', sys_include_sets, sys_priors,
                                group_cols=['divs', 'energy', 'cent', 'data_type', 'total_protons'])
        print(df_def_dsigma)
        df_def_dsigma.to_csv(f'{base_path}{df_def_dsigma_out_name}', index=False)

    # Calculate v2 subtraction for each total_protons value
    if df_def_dsigma_v2sub_out_name is not None and calc_finals:
        print('Calc v2 sub')
        df_dsigma_types['meas'] = df_dsigma_types.apply(lambda row: Measure(row['val'], row['err']), axis=1)
        df_def_dsigma_v2sub = []
        for set_name in sys_include_names + ['bes_def']:
            set_v2sub = subtract_dsigma_flow(df_dsigma_types, set_name, set_name, v2_sys_vals[set_name],
                                             new_only=True, val_col='val', err_col='err', meas_col='meas')
            df_def_dsigma_v2sub.append(set_v2sub)
        df_def_dsigma_v2sub = pd.concat(df_def_dsigma_v2sub)
        df_def_dsigma_v2sub = get_sys(df_def_dsigma_v2sub, 'bes_def', sys_include_sets, sys_priors,
                                      group_cols=['divs', 'energy', 'cent', 'data_type', 'total_protons'])
        df_def_dsigma_v2sub.to_csv(f'{base_path}{df_def_dsigma_v2sub_out_name}', index=False)

    sets_run = all_sets if plot else sys_include_names + ['bes_def']
    dsig_avgs_all, dsig_avgs_diff_v2sub = [], []
    print(sets_run)
    for energy in energies_fit:
        print(f'Energy {energy}GeV')
        jobs = [(df_dsigma_types[df_dsigma_types['name'] == set_i], divs_all, cents, energy, ['raw', 'mix', 'diff'],
                 [set_i], None, False, True) for set_i in sets_run]
        dsig_avgs_div_all = []
        with Pool(threads) as pool:
            for job_i in tqdm.tqdm(pool.istarmap(dvar_vs_protons_cents, jobs), total=len(jobs)):
                dsig_avgs_div_all.append(job_i)
        dsig_avgs_div_all = pd.concat(dsig_avgs_div_all, ignore_index=True)
        dsig_avgs_all.append(dsig_avgs_div_all)
        dsig_avgs_div_diff = dsig_avgs_div_all[dsig_avgs_div_all['data_type'] == 'diff']
        dsig_avgs_div_diff = dsig_avgs_div_diff.drop('data_type', axis=1)
        for data_set in sets_run:
            dsig_avgs_div_diff_set = subtract_dsigma_flow(dsig_avgs_div_diff, data_set,
                                                          data_set, v2_sys_vals[data_set], new_only=True)
            dsig_avgs_diff_v2sub.append(dsig_avgs_div_diff_set)

    if df_def_avgs_out_name is not None and calc_finals:
        dsig_avg_all = pd.concat(dsig_avgs_all, ignore_index=True)
        dsig_avgs_def_sys = get_sys(dsig_avg_all, 'bes_def', sys_include_sets, sys_priors, val_col='avg',
                                    err_col='avg_err', group_cols=['divs', 'energy', 'cent', 'data_type'])
        dsig_avgs_def_sys.to_csv(f'{base_path}{df_def_avgs_out_name}', index=False)

    dsig_avgs_diff_v2sub = pd.concat(dsig_avgs_diff_v2sub, ignore_index=True)
    # if plot:
    #     # dsig_avgs_diff_v2sub = dsig_avgs_diff_v2sub[(dsig_avgs_diff_v2sub['divs'] == 120) &
    #     #                                             (dsig_avgs_diff_v2sub['cent'] == 8)]
    #     plot_sys(dsig_avgs_diff_v2sub, 'bes_def', non_rand_sets, sys_info_dict, val_col='avg', err_col='avg_err',
    #              group_cols=['divs', 'energy', 'cent'], y_label=r'$\langle \Delta \sigma^2 \rangle$',
    #              # pdf_out_path=None)
    #              pdf_out_path=sys_pdf_out_path)
    #     # plot_sys(dsig_avgs_diff_v2sub, 'bes_def', rand_sets, sys_info_dict, val_col='avg', err_col='avg_err',
    #     #          group_cols=['divs', 'energy', 'cent'], plot_bars=False, y_label=r'$\Delta \sigma^2$',
    #     #          pdf_out_path=None)
    #     #          # pdf_out_path=sys_pdf_out_path.replace('.pdf', '_rands.pdf'))
    #     # plt.show()

    if df_def_avgs_v2sub_out_name is not None and calc_finals:
        dsig_avg_diff_v2sub_out = get_sys(dsig_avgs_diff_v2sub, 'bes_def', sys_include_sets, sys_priors,
                                          val_col='avg', err_col='avg_err', group_cols=['divs', 'energy', 'cent'])
        dsig_avg_diff_v2sub_out.to_csv(f'{base_path}{df_def_avgs_v2sub_out_name}', index=False)

    dsig_avgs_v2_sub_div120 = dsig_avgs_diff_v2sub[dsig_avgs_diff_v2sub['divs'] == 120]
    cent_fits_120_refmult = plot_dsig_avg_vs_cent_fit(dsig_avgs_v2_sub_div120, all_sets, fit=True, cent_ref=cent_ref_df,
                                                        ref_type=ref_type, plot=False)
    cent_fits_120_npart = plot_dsig_avg_vs_cent_fit(dsig_avgs_v2_sub_div120, all_sets, fit=True, cent_ref=cent_ref_df,
                                                      ref_type='npart', plot=False)
    param_names = ['a', 'b', 'c']
    for ref_type_i, ref_type_df in zip(['refmult', 'npart'], [cent_fits_120_refmult, cent_fits_120_npart]):
        for param_name in param_names:
            cent_fits_120_ref_type = ref_type_df.copy()
            for other_param_name in param_names:
                if other_param_name == param_name:
                    continue
                cent_fits_120_ref_type.drop(columns=[other_param_name, f'{other_param_name}_err'], inplace=True)
            cent_fits_ref_param_out = get_sys(cent_fits_120_ref_type, 'bes_def', sys_include_sets, sys_priors,
                                          val_col=param_name, err_col=f'{param_name}_err', group_cols=['energy'],
                                          name_col='data_set')
            cent_fits_ref_param_out.to_csv(f'{base_path}/Bes_with_Sys/centrality_fits/{ref_type_i}_fits_120_{param_name}.csv',
                                           index=False)
    input('Enter to continue...')

    # plot_sys(dsig_avgs_diff_v2sub, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='avg', err_col='avg_err',
    #          group_cols=['divs', 'energy', 'cent'], y_label=r'$\langle \Delta \sigma^2 \rangle$',
    #          pdf_out_path=sys_pdf_out_path, indiv_pdf_path=indiv_pdf_out_path)
    plot_sys(dsig_avgs_diff_v2sub, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='avg',
             err_col='avg_err', group_cols=['cent', 'energy', 'divs'], y_label=r'$\langle \Delta \sigma^2 \rangle$',
             pdf_out_path=sys_pdf_out_path, indiv_pdf_path=indiv_pdf_out_path)
    df_120_4_up = dsig_avgs_diff_v2sub[(dsig_avgs_diff_v2sub['cent'] > 4) & (dsig_avgs_diff_v2sub['divs'] == 120)]
    plot_sys_table(df_120_4_up, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='avg',
                   err_col='avg_err', name_col='name', indiv_pdf_path=f'{indiv_pdf_out_path}sys_table_high_cent.pdf')
    df_120_4_down = dsig_avgs_diff_v2sub[(dsig_avgs_diff_v2sub['cent'] <= 4) & (dsig_avgs_diff_v2sub['divs'] == 120)]
    plot_sys_table(df_120_4_down, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='avg',
                   err_col='avg_err', name_col='name', indiv_pdf_path=f'{indiv_pdf_out_path}sys_table_low_cent.pdf')

    df_fits = plot_dvar_avgs_divs(dsig_avgs_diff_v2sub, all_sets, fit=True, plot_energy_panels=False, plot=False)
    if df_partitions_fits_name is not None:
        df_baselines = get_sys(df_fits, 'bes_def', sys_include_sets,  sys_priors, val_col='baseline',
                               err_col='base_err', name_col='data_set', sys_col='base_sys',
                               group_cols=['energy', 'cent'])
        df_zeros = get_sys(df_fits, 'bes_def', sys_include_sets, sys_priors, val_col='zero_mag',
                           err_col='zero_mag_err', name_col='data_set', sys_col='zero_sys',
                           group_cols=['energy', 'cent'])
        df_baselines = df_baselines.drop(columns=['zero_mag', 'zero_mag_err'])
        df_fits_out = pd.concat([df_baselines, df_zeros[['zero_mag', 'zero_mag_err', 'zero_sys']]], axis=1)
        df_fits_out.to_csv(f'{base_path}{df_partitions_fits_name}', index=False)

        print(df_fits_out)
        if plot:
            plot_slope_div_fits(df_fits_out, data_sets_colors, data_sets_labels)
            plot_slope_div_fits_simpars(df_fits_out)
    df_fits = df_fits.rename(columns={'data_set': 'name'})
    plot_sys(df_fits, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='baseline', err_col='base_err',
             group_cols=['cent', 'energy'], name_col='name', indiv_pdf_path=indiv_pdf_out_path,
             pdf_out_path=sys_pdf_out_path.replace('.pdf', '_baseline.pdf'), y_label='Baseline')
    plot_sys(df_fits, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='zero_mag',
             err_col='zero_mag_err', group_cols=['cent', 'energy'], name_col='name',
             pdf_out_path=sys_pdf_out_path.replace('.pdf', '_zeros.pdf'), y_label='Curvature')
    df_fits_4_up = df_fits[df_fits['cent'] > 4]
    plot_sys_table(df_fits_4_up, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='baseline',
                   err_col='base_err', name_col='name',
                   indiv_pdf_path=f'{indiv_pdf_out_path}baseline_sys_table_high_cent.pdf')
    df_fits_4_down = df_fits[df_fits['cent'] <= 4]
    plot_sys_table(df_fits_4_down, 'bes_def', non_rand_sets, sys_info_dict, sys_priors, val_col='baseline',
                   err_col='base_err', name_col='name',
                   indiv_pdf_path=f'{indiv_pdf_out_path}baseline_sys_table_low_cent.pdf')

    plt.show()


def make_models_csv():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    v2_cf_in_dir = 'F:/Research/Data_CF/default_resample_epbins1/CF_rapid05_resample_norotate_epbins1_0/'
    v2_cfev_in_dir = 'F:/Research/Data_CFEV/default_resample_epbins1/CFEV_rapid05_resample_norotate_epbins1_0/'
    v2_cfevb342_in_dir = 'F:/Research/Data_CFEVb342/default_resample_epbins1/' \
                         'CFEVb342_rapid05_resample_norotate_epbins1_0/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_var_epbins1.csv'

    threads = 15
    df_def_out_name = 'Bes_with_Sys/binom_slice_vars_model.csv'
    # df_def_out_name = None
    df_def_dsigma_out_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma.csv'
    # df_def_dsigma_out_name = None
    df_def_dsigma_v2sub_out_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma_v2sub.csv'
    df_def_avgs_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_model.csv'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'
    df_def_avgs_v2sub_raw_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_raw_model.csv'
    df_partitions_fits_name = 'Bes_with_Sys/partition_width_fits_model.csv'
    df_partitions_fits_raw_name = 'Bes_with_Sys/partition_width_fits_raw_model.csv'
    df_partitions_fits_nov2sub_name = 'Bes_with_Sys/partition_width_fits_nov2sub_model.csv'

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    # divs_all = [60, 120, 180]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [7, 11]

    df = pd.read_csv(f'{base_path}{df_name}')
    df = df.dropna()
    df = df[df['name'] != 'bes_resample_epbins1']
    all_sets = pd.unique(df['name'])
    print(all_sets)

    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}
    v2_rp_ampt_vals = {2: read_flow_values(v2_ampt_in_dir, 'v2_rp')}
    v2_cum_ampt_vals = {1: read_flow_values(v2_ampt_in_dir, 'v1_cum_new'),
                        2: read_flow_values(v2_ampt_in_dir, 'v2_cum_new')}
    v2_cf_vals = {2: read_flow_values(v2_cf_in_dir)}
    v2_cfev_vals = {2: read_flow_values(v2_cfev_in_dir)}
    v2_cfevb342_vals = {2: read_flow_values(v2_cfevb342_in_dir)}
    v2_sys_vals = {'ampt_new_coal_epbins1': v2_ampt_vals, 'cf_resample_epbins1': v2_cf_vals,
                   'cfev_resample_epbins1': v2_cfev_vals, 'cfevb342_resample_epbins1': v2_cfevb342_vals,
                   'ampt_new_coal_epbins1_v2rp': v2_rp_ampt_vals, 'ampt_new_coal_epbins1_v2cum': v2_cum_ampt_vals}

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    df = df[df['stat'] == stat_plot]
    df['sys'] = 0

    # Get k2 raw, mix, diff systematics
    df.to_csv(f'{base_path}{df_def_out_name}', index=False)

    # Calculate dsigma with k2 values
    df = df[df['stat'] == 'k2']
    df = df.drop('stat', axis=1)
    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])
    print('Calc dsigma')
    df_dsigma_types = pd.concat([df_raw, df_mix, df_diff])

    # Add v2rp dataset which is just copy of ampt_new_coal_epbins1
    df_ampt_v2rp = df_dsigma_types[df_dsigma_types['name'] == 'ampt_new_coal_epbins1'].copy()
    df_ampt_v2rp['name'] = 'ampt_new_coal_epbins1_v2rp'
    df_ampt_v2cum = df_dsigma_types[df_dsigma_types['name'] == 'ampt_new_coal_epbins1'].copy()
    df_ampt_v2cum['name'] = 'ampt_new_coal_epbins1_v2cum'
    all_sets = np.append(all_sets, ['ampt_new_coal_epbins1_v2rp', 'ampt_new_coal_epbins1_v2cum'])
    df_dsigma_types = pd.concat([df_dsigma_types, df_ampt_v2rp, df_ampt_v2cum], ignore_index=True)

    print('Calc diff nlo error')
    df_dsigma_types = add_diff_nlo_err(df_dsigma_types, group_cols=['energy', 'cent', 'name', 'total_protons'],
                                       exclude_divs=[356, 89])
    df_dsigma_types.to_csv(f'{base_path}{df_def_dsigma_out_name}', index=False)

    # Calculate dsigma with v2 subtracted
    # df_dsigma_types['meas'] = df_dsigma_types.apply(lambda row: Measure(row['val'], row['err']), axis=1)
    # df_def_dsigma_v2sub = subtract_dsigma_flow(df_dsigma_types, 'bes_def', 'bes_def', v2_sys_vals['bes_def'],
    #                                            new_only=True, val_col='val', err_col='err', meas_col='meas')
    # df_def_dsigma_v2sub.to_csv(f'{base_path}{df_def_dsigma_v2sub_out_name}', index=False)

    # Calculate v2 subtraction for each total_protons value
    df_dsigma_types['meas'] = df_dsigma_types.apply(lambda row: Measure(row['val'], row['err']), axis=1)
    df_def_dsigma_v2sub = []
    for set_name in all_sets:
        set_v2sub = subtract_dsigma_flow(df_dsigma_types, set_name, set_name, v2_sys_vals[set_name],
                                         new_only=True, val_col='val', err_col='err', meas_col='meas')
        df_def_dsigma_v2sub.append(set_v2sub)
    df_def_dsigma_v2sub = pd.concat(df_def_dsigma_v2sub)
    df_def_dsigma_v2sub.to_csv(f'{base_path}{df_def_dsigma_v2sub_out_name}', index=False)

    sets_run = all_sets
    dsig_avgs_all, dsig_avgs_diff_v2sub, dsig_avgs_raw_v2sub = [], [], []
    print(sets_run)
    for energy in energies_fit:
        print(f'Energy {energy}GeV')
        jobs = [(df_dsigma_types[df_dsigma_types['name'] == set_i], divs_all, cents, energy, ['raw', 'mix', 'diff'],
                 [set_i], None, False, True) for set_i in sets_run]
        dsig_avgs_div_all = []
        with Pool(threads) as pool:
            for job_i in tqdm.tqdm(pool.istarmap(dvar_vs_protons_cents, jobs), total=len(jobs)):
                dsig_avgs_div_all.append(job_i)
        dsig_avgs_div_all = pd.concat(dsig_avgs_div_all, ignore_index=True)
        dsig_avgs_div_all['sys'] = 0
        dsig_avgs_all.append(dsig_avgs_div_all)
        dsig_avgs_div_diff = dsig_avgs_div_all[dsig_avgs_div_all['data_type'] == 'diff']
        dsig_avgs_div_diff = dsig_avgs_div_diff.drop('data_type', axis=1)
        dsig_avgs_div_raw = dsig_avgs_div_all[dsig_avgs_div_all['data_type'] == 'raw']
        dsig_avgs_div_raw = dsig_avgs_div_raw.drop('data_type', axis=1)
        for data_set in sets_run:
            dsig_avgs_div_diff_set = subtract_dsigma_flow(dsig_avgs_div_diff, data_set,
                                                          data_set, v2_sys_vals[data_set], new_only=True)
            dsig_avgs_diff_v2sub.append(dsig_avgs_div_diff_set)
            dsig_avgs_div_raw_set = subtract_dsigma_flow(dsig_avgs_div_raw, data_set,
                                                         data_set, v2_sys_vals[data_set], new_only=True)
            dsig_avgs_raw_v2sub.append(dsig_avgs_div_raw_set)

    dsig_avg_all = pd.concat(dsig_avgs_all, ignore_index=True)
    dsig_avg_all.to_csv(f'{base_path}{df_def_avgs_out_name}', index=False)

    dsig_avgs_diff_v2sub = pd.concat(dsig_avgs_diff_v2sub, ignore_index=True)
    dsig_avgs_diff_v2sub.to_csv(f'{base_path}{df_def_avgs_v2sub_out_name}', index=False)

    dsig_avgs_raw_v2sub = pd.concat(dsig_avgs_raw_v2sub, ignore_index=True)
    dsig_avgs_raw_v2sub.to_csv(f'{base_path}{df_def_avgs_v2sub_raw_out_name}', index=False)

    df_fits = plot_dvar_avgs_divs(dsig_avgs_diff_v2sub, all_sets, fit=True, plot=True)
    df_fits.to_csv(f'{base_path}{df_partitions_fits_name}', index=False)

    df_fits_raw = plot_dvar_avgs_divs(dsig_avgs_raw_v2sub, all_sets, fit=True, plot=True)
    df_fits_raw.to_csv(f'{base_path}{df_partitions_fits_raw_name}', index=False)

    dsig_avgs_raw = dsig_avg_all[dsig_avg_all['data_type'] == 'raw']
    dsig_avgs_raw = dsig_avgs_raw.drop('data_type', axis=1)

    df_fits_nov2sub = plot_dvar_avgs_divs(dsig_avgs_raw, all_sets, fit=True, plot=True)
    df_fits_nov2sub.to_csv(f'{base_path}{df_partitions_fits_nov2sub_name}', index=False)

    plt.show()


def plot_star_var_rand_sys():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
                     'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    df_name = 'Binomial_Slice_Moments/binom_slice_vars_bes_rand_sys_test.csv'
    # fits_out_base = 'Base_Zero_Fits'
    # df_tproton_avgs_name = 'dsig_tprotons_avgs_cent8.csv'
    # df_partitions_fits_name = 'partitions_fits_cent8.csv'
    df_path = base_path + df_name

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 180
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    # energies_fit = [7, 11, 19, 27, 39, 62]
    energies_fit = [7]
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_def', 'bes_sys_dca08', 'bes_sys_dca12']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'blue', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['dca = 1.0 cm', 'dca = 0.8 cm', 'dca = 1.2 cm']))

    # data_sets_plt = ['bes_def', 'bes_sys_nsprx09', 'bes_sys_nsprx11']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'blue', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['nsigmaproton = 2.0', 'nsigmaproton = 1.8', 'nsigmaproton = 2.2']))
    #
    # data_sets_plt = ['bes_def', 'bes_sys_nsprx09', 'bes_sys_nsprx11']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'blue', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['nsigmaproton = 2.0', 'nsigmaproton = 1.8', 'nsigmaproton = 2.2']))

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}

    df = pd.read_csv(df_path)
    df = df.dropna()
    all_sets = pd.unique(df['name'])
    print(all_sets)
    print(df)

    df = df[df['stat'] == stat_plot]
    df_raw = calc_dsigma(df, 'raw')
    print(df_raw)

    dsig_avg = dvar_vs_protons(df_raw, div_plt, cent_plt, [7], ['raw'], all_sets, plot=True, avg=True,
                               leg=False)

    # print(dsig_avg)
    dsig_avg['set_num'] = dsig_avg['name'].apply(lambda x: int(x.split('_')[-1]))
    dsig_avg['set_group'] = dsig_avg['name'].apply(lambda x: x.replace('_' + x.split('_')[-1], ''))
    # print(dsig_avg)
    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    ax.axhline(np.average(dsig_avg['avg']), color='gray')
    ax.set_title(r'Variation of $\langle\Delta\sigma^2\rangle$ from Randomization Sources')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$')
    ax.set_xlabel('Arbitrary run index')
    fig2, ax2 = plt.subplots()
    # plt.xticks(rotation=20)
    ax2.grid()
    ax2.set_ylabel(r'Standard Deviation of $\langle\Delta\sigma^2\rangle$')
    iterate_set_num = 0
    group_map = {'bes_rand_sys': 'File Order, Resampling, StRefMultCorr', 'bes_strefnoseed_sys': 'StRefMultCorr',
                 'bes_strefseed_sys': 'Resampling'}
    for group in pd.unique(dsig_avg['set_group']):
        group_df = dsig_avg[dsig_avg['set_group'] == group]
        # print(f'{group}\n{group_df}')
        ebar = ax.errorbar(group_df['set_num'] + iterate_set_num, group_df['avg'], yerr=group_df['avg_err'], ls='',
                           marker='o', label=group_map[group])
        ax.plot(group_df['set_num'] + iterate_set_num, np.ones(len(group_df)) * np.average(group_df['avg']),
                color=ebar[0].get_color())
        ax2.scatter(['\n'.join(group_map[group].split(', '))], [np.std(group_df['avg'])], color=ebar[0].get_color(),
                    zorder=10)
        iterate_set_num += max(group_df['set_num']) + 1
    ax.legend()
    ax2.axhline(0, color='gray')
    fig.tight_layout()
    fig2.tight_layout()

    plt.show()

    # dsig_avgs = []
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     print(f'Div {div}')
    #     dsig_avgs_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], data_sets_plt, plot=False,
    #                                         avg=True)
    #     dsig_avgs_div_raw.loc[:, 'name'] = dsig_avgs_div_raw['name'] + '_raw'
    #     dsig_avgs.append(dsig_avgs_div_raw)
    #
    #     dsig_avgs_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], data_sets_plt, plot=False,
    #                                         avg=True)
    #     dsig_avgs_div_mix.loc[:, 'name'] = dsig_avgs_div_mix['name'] + '_mix'
    #     dsig_avgs.append(dsig_avgs_div_mix)
    #
    #     dsig_avgs_div_diff = dvar_vs_protons(df_diff, div, cent_plt, energies_fit, ['diff'], data_sets_plt, plot=False,
    #                                          avg=True)
    #     dsig_avgs_div_diff.loc[:, 'name'] = dsig_avgs_div_diff['name'] + '_sub'
    #
    #     dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'bes_resample_epbins1_sub',
    #                                               'bes_resample_epbins1', v2_star_vals, div, cent_plt)
    #     dsig_avgs.append(dsig_avgs_div_diff)
    #
    # dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    # # if df_tproton_avgs_name is not None:
    # #     dsig_avgs.to_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}', index=False)
    #
    # for data_set in data_sets_plt:
    #     data_sets = [data_set + x for x in ['_raw', '_mix', '_sub']]
    #     colors = dict(zip(data_sets, ['blue', 'green', 'red']))
    #     labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Raw', ' Mix', ' Sub']]))
    #     plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels,
    #                         ylab=r'$\langle\Delta\sigma^2\rangle$')
    # for data_set in data_sets_plt:
    #     data_sets = [data_set + x for x in ['_sub', '']]
    #     colors = dict(zip(data_sets, ['blue', 'red']))
    #     labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Original', ' v2 Corrected']]))
    #     plot_dvar_avgs_divs(dsig_avgs, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels)
    # plot_dvar_avgs_divs(dsig_avgs, data_sets_plt, data_sets_colors=data_sets_colors, fit=True,
    #                     data_sets_labels=data_sets_labels)
    #
    # df_fits = plot_dvar_avgs_divs(dsig_avgs, data_sets_plt, data_sets_colors=data_sets_colors, fit=True,
    #                               data_sets_labels=data_sets_labels, plt_energies=False)
    # # if df_partitions_fits_name is not None:
    # #     df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)
    #
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)
    #
    # plt.show()


def plot_sims():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # df_name = 'binom_slice_stats_cent8_sim.csv'
    df_name = 'binom_slice_stats_cent8_sim_test.csv'
    df_path = base_path + df_name
    sim_sets = []

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    amps = ['002', '006', '01']  # ['002', '006', '01']
    # spreads = ['03', '04', '05', '06', '07', '08', '09', '1', '11', '12']
    spreads = ['1', '12']
    # spreads = ['1', '15']
    # amps = ['0', '002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    #         '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5',
    #         '6', '7', '8', '9', '99']
    # spreads = ['001', '01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25',
    #            '275', '3', '325', '35', '4']  # '375' missing '08' and '09'
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    # stat_plot = 'skewness'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 7]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    if stat_plot == 'skewness':
        exclude_divs.append(180)

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={spread}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    # plt.show()
    proton_fits_div = protons_fits[protons_fits['divs'] == div_plt]
    plot_protons_fits_vs_amp(proton_fits_div, all_sets_plt, data_sets_colors, data_sets_labels)
    # plt.show()
    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    plt.show()
    plot_sigma_fits_interp(sigma_fits)

    plt.show()


def plot_sims_var():
    # plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.figsize"] = (8, 4)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # df_name = 'binom_slice_stats_sim.csv'
    df_name = 'binom_slice_stats_sim_demos.csv'
    df_path = base_path + df_name
    df_dsig_out_path = base_path + 'binom_slice_stats_sim_demos_dsigma.csv'
    df_dsig_avg_out_path = base_path + 'binom_slice_stats_sim_demos_dsigma_avg.csv'
    df_width_fits_out_path = base_path + 'binom_slice_stats_sim_demos_width_fits.csv'
    write_out_dfs = True

    sim_sets = []

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    spreads = ['04', '05', '06', '07', '08', '09', '1', '11', '12']

    # amps = ['006', '01']  # ['002', '006', '01']

    # amps = ['002', '008']
    # spreads = ['04', '08']

    # spreads = ['03', '04', '05', '06', '07', '08', '09', '1', '11', '12']
    # spreads = ['1']
    # amps = ['0', '002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    #         '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5',
    #         '6', '7', '8', '9', '99']
    # spreads = ['001', '01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25',
    #            '275', '3', '325', '35', '4']  # '375' missing '08' and '09'
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'k2'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 7]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={spread}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw, df_mix, df_diff = calc_dsigma(df, data_types=['raw', 'mix', 'diff'])

    dvar_vs_protons(df_raw, div_plt, cent_plt, ['sim'], ['raw'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels, ylabel=r'$\Delta\sigma^2$', data_sets_bands=all_sets_plt,
                    title=f'Gaussian Correlation Model: {div_plt}° Partitions, 72 Samples per Event')
    if write_out_dfs:
        df_combined = pd.concat([df_raw, df_mix, df_diff])
        df_combined.to_csv(df_dsig_out_path, index=False)

    dsig_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsig_avgs_div = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets_plt, plot=False,
                                        avg=True)
        dsig_avgs.append(dsig_avgs_div)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    print(dsig_avgs)
    print(pd.unique(dsig_avgs['amp']))
    print(pd.unique(dsig_avgs['spread']))
    plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels, plot_energy_panels=False,
                        ylab=r'$\langle\Delta\sigma^2\rangle$',
                        title=f'Gaussian Correlation Model: 72 Samples per Event')

    if write_out_dfs:
        dsig_avgs.to_csv(df_dsig_avg_out_path, index=False)

    df_fits = plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                  data_sets_labels=data_sets_labels, plot_energy_panels=False,
                                  ylab=r'$\langle\Delta\sigma^2\rangle$',
                                  title=f'Gaussian Correlation Model: 72 Samples per Event')
    # plt.show()
    proton_fits_div = dsig_avgs[dsig_avgs['divs'] == div_plt]
    print(proton_fits_div)
    plot_dsig_fits_vs_amp(proton_fits_div, all_sets_plt, data_sets_colors, data_sets_labels,
                          title=f'Gaussian Correlation Model: {div_plt}° Partitions, 72 Samples per Event')
    if write_out_dfs:
        df_fits.to_csv(df_width_fits_out_path, index=False)

    # # plt.show()
    # # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    # plt.show()
    # plot_sigma_fits_interp(sigma_fits)

    plt.show()


def get_sim_mapping():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_sim.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'sim_tprotons_fits.csv'
    df_partitions_fits_name = 'sim_partitions_fits.csv'
    threads = 11
    df_path = base_path + df_name
    sim_sets = []

    print(f'{base_path}{fits_out_base}{df_tproton_fits_name}')

    # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    # spreads = ['02', '05', '06', '065', '07']
    # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04']
    # spreads = ['05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
            '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5']
    spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25', '275', '3',
               '325', '35', '375', '4']
    for amp in amps:
        for spread in spreads:
            sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
            # sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    sim_sets = sorted(sim_sets, reverse=True)
    sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={round(spread, 2)}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, False, True, False, None, None)
    jobs = [(df, stat_plot, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    with Pool(threads) as pool:
        for proton_fits_div in tqdm.tqdm(pool.istarmap(stat_vs_protons, jobs), total=len(jobs)):
            # print(proton_fits_div)
            protons_fits.append(proton_fits_div)
    print(protons_fits)
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     print(f'Div {div}')
    #     protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
    #                                        plot=False, fit=True)
    #     protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    protons_fits.to_csv(f'{base_path}{fits_out_base}/{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    interpolations = plot_sigma_fits_interp(sigma_fits)

    with open(f'{base_path}{pickle_map_name}', 'ab') as pickle_file:
        pickle.dump(interpolations, pickle_file)

    # with open(f'{base_path}{pickle_map_name}', 'rb') as pickle_file:
    #     interps_pickle = pickle.load(pickle_file)

    for fig_i in plt.get_fignums():
        plt.figure(fig_i)
        window_title = plt.gcf().canvas.manager.get_window_title()
        plt.savefig(f'{out_dir}{window_title}.png')

    plt.show()


def get_sim_mapping_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_sim.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_avgs_name = 'sim_dsig_tprotons_avgs.csv'
    df_partitions_fits_name = 'sim_partitions_fits.csv'
    threads = 15
    df_path = base_path + df_name
    sim_sets = []

    print(f'{base_path}{fits_out_base}{df_tproton_avgs_name}')

    stat_plot = 'k2'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim', 62]
    data_types_plt = ['diff']
    samples = 72  # For title purposes only

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # # amps = ['002', '004', '006', '008', '01']  # ['002', '006', '01']
    # # spreads = ['02', '05', '06', '065', '07']
    # # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04']
    # # spreads = ['05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    #         '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5']
    # spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25', '275', '3',
    #            '325', '35', '375', '4']
    # # amps = ['002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
    # #         '045', '05']
    # # spreads = ['01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2']
    # for amp in amps:
    #     for spread in spreads:
    #         sim_sets.append(f'sim_aclmul_amp{amp}_spread{spread}')
    #         # sim_sets.append(f'sim_clmul_amp{amp}_spread{spread}')
    # sim_sets = sorted(sim_sets, reverse=True)
    # sim_sets = sim_sets[:int(len(sim_sets) / 2)] + sorted(sim_sets[int(len(sim_sets) / 2):])

    sim_sets = []
    for data_set in pd.unique(df['name']):
        if 'sim_' in data_set:
            sim_sets.append(data_set)

    print(sim_sets)
    # sim_sets = sim_sets[:20]

    data_sets_labels = {}
    for sim_set in sim_sets:
        label = ''
        if '_clmul_' in sim_set:
            label += 'Attractive '
        elif '_aclmul_' in sim_set:
            label += 'Repulsive '
        amp, spread = get_name_amp_spread(sim_set)
        label += f'A={amp} σ={round(spread, 2)}'
        data_sets_labels.update({sim_set: label})

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    all_sets_plt = data_sets_plt + sim_sets[:]

    df_diff = calc_dsigma(df[df['stat'] == stat_plot])

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    dsig_avgs = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, False, True, False, None, None, None)
    jobs = [(df_diff, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    with Pool(threads) as pool:
        for dsig_avgs_div in tqdm.tqdm(pool.istarmap(dvar_vs_protons, jobs), total=len(jobs)):
            # print(dsig_avgs_div)
            dsig_avgs.append(dsig_avgs_div)
    print(dsig_avgs)
    # for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
    #     print(f'Div {div}')
    #     protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
    #                                        plot=False, fit=True)
    #     protons_fits.append(protons_fits_div)
    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    if df_tproton_avgs_name is not None:
        dsig_avgs.to_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}', index=False)
    print(dsig_avgs)
    print(pd.unique(dsig_avgs['amp']))
    print(pd.unique(dsig_avgs['spread']))
    df_fits = plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                  data_sets_labels=data_sets_labels, plot_energy_panels=False)
    if df_partitions_fits_name is not None:
        df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)
    interpolations = plot_sigma_fits_interp(sigma_fits)
    #
    # with open(f'{base_path}{pickle_map_name}', 'ab') as pickle_file:
    #     pickle.dump(interpolations, pickle_file)
    #
    # # with open(f'{base_path}{pickle_map_name}', 'rb') as pickle_file:
    # #     interps_pickle = pickle.load(pickle_file)
    #
    # for fig_i in plt.get_fignums():
    #     plt.figure(fig_i)
    #     window_title = plt.gcf().canvas.manager.get_window_title()
    #     plt.savefig(f'{out_dir}{window_title}.png')

    plt.show()


def get_sim_mapping_pm():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_simpm_test.csv'
    pickle_map_name = 'binom_slice_sim_mapping'
    out_dir = 'F:/Research/Results/Sim_Mapping/'
    df_fits_out_name = 'Base_Zero_Fits/simpm.csv'
    threads = 11
    df_path = base_path + df_name
    sim_sets = []

    aps = ['04']  # ['02', '04']
    ams = ['02']  # ['01', '02']
    sps = ['02', '04']  # ['02', '04']
    sms = ['2']  # ['1', '2']
    for (ap, am, sp, sm) in itertools.product(aps, ams, sps, sms):
        sim_sets.append(f'sim_clmultipm_ampplus{ap}_ampminus{am}_spreadplus{sp}_spreadminus{sm}')

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_plt = [39, 'sim']  # [7, 11, 19, 27, 39, 62, 'sim']  # [7, 11, 19, 27, 39, 62]
    energies_fit = ['sim']
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt, data_sets_colors, data_sets_labels = [], None, None
    data_sets_plt = []
    data_sets_labels = {}
    for sim_set in sim_sets:
        ap, am, sp, sm = get_name_amp_spread_pm(sim_set)
        label = f'A+={ap} σ+={round(sp, 2)} A-={am} σ-={round(sm, 2)}'
        data_sets_labels.update({sim_set: label})

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # stat_vs_protons(df, stat_plot, div_plt, cent_plt, [7, 'sim'], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                 data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = []
    other_pars = (cent_plt, energies_fit, data_types_plt, all_sets_plt, None, True, True, False, None, None)
    jobs = [(df, stat_plot, div, *other_pars) for div in np.setdiff1d(np.unique(df['divs']), exclude_divs)]
    stat_vs_protons(*jobs[0])
    stat_vs_protons(*jobs[1])
    with Pool(threads) as pool:
        for proton_fits_div in tqdm.tqdm(pool.istarmap(stat_vs_protons, jobs), total=len(jobs)):
            protons_fits.append(proton_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    print(protons_fits)
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{df_fits_out_name}', index=False)
    # print(df_fits)
    sigma_fits = plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_star_model_onediv():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cent8_no_sim.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies = [7, 11, 19, 27, 39, 62]
    energy_plt = 39
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cf_resample_def', 'cfev_resample_def',
                     'cfevb342_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'olive']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$',
                                                'MUSIC+FIST EV $3.42fm^3$']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [energy_plt], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    protons_fits = stat_vs_protons(df, stat_plot, div_plt, cent_plt, energies, data_types_plt,
                                   all_sets_plt, plot=False, fit=True)
    print(protons_fits)
    plot_protons_fits_vs_energy(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, title=f'0-5% Centrality, {div_plt}° Partitions, '
                                                                         f'{samples} Samples per Event')

    plt.show()


def plot_vs_cent():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    save_fits = False
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = 'bes_ampt_cents_tprotons_fits.csv'
    df_partitions_fits_name = 'bes_ampt_cents_partitions_fits.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    # stat_plot = 'non-excess kurtosis'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 39
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    if stat_plot == 'skewness':
        divs_all.pop(180)
        exclude_divs.append(180)

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    # data_sets_plt = ['bes_resample_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black']))
    # data_sets_energies_cmaps = dict(zip(data_sets_plt, ['brg']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, 8, [energy_plt], ['raw', 'mix'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    stat_vs_protons_cents(df, stat_plot, [div_plt], cents, energy_plt, data_types_plt, all_sets_plt, plot=True,
                          fit=True, plot_fit=True, data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)

    df_slope_fits = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy_plt, data_types_plt, all_sets_plt,
                                          plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                          data_sets_labels=data_sets_labels)

    plot_protons_fits_divs_cents(df_slope_fits, data_sets_plt, fit=True, plot=True, data_sets_colors=data_sets_colors,
                                 data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)

    df_onediv_fits, df_tp_fits, df_divs_fits = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    for energy in energies_fit:
        df_fit = stat_vs_protons_cents(df, stat_plot, [div_plt], cents, energy, data_types_plt, all_sets_plt,
                                       plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                       data_sets_labels=data_sets_labels)
        df_onediv_fits = pd.concat([df_onediv_fits, df_fit])

        df_tp_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy, data_types_plt, all_sets_plt,
                                          plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                          data_sets_labels=data_sets_labels)
        df_tp_fits = pd.concat([df_tp_fits, df_tp_fit])

        df_fits = plot_protons_fits_divs_cents(df_tp_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                               data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)
        # print(f'{energy}GeV {df_fit}')
        df_divs_fits = pd.concat([df_divs_fits, df_fits])

    if save_fits:
        df_tp_fits.to_csv(f'{base_path}{fits_out_base}/{df_tproton_fits_name}', index=False)
        df_divs_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=None,
                          ref_type=None, data_sets_energies_cmaps=data_sets_energies_cmaps)
    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_vs_cent_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    v2_star_in_dir = 'F:/Research/Data/default_resample_epbins1_calcv2_qaonly_test/' \
                     'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_qaonly_test_0/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    v2_cf_in_dir = 'F:/Research/Data_CF/default_resample_epbins1/CF_rapid05_resample_norotate_epbins1_0/'
    v2_cfev_in_dir = 'F:/Research/Data_CFEV/default_resample_epbins1/CFEV_rapid05_resample_norotate_epbins1_0/'
    v2_cfevb342_in_dir = 'F:/Research/Data_CFEVb342/default_resample_epbins1/' \
                         'CFEVb342_rapid05_resample_norotate_epbins1_0/'
    df_name = 'Binomial_Slice_Moments/binom_slice_vars_bes.csv'
    fits_out_base = 'Base_Zero_Fits'
    # df_tproton_avgs_name = 'dsig_tprotons_avgs.csv'
    df_tproton_avgs_name = None
    # df_partitions_fits_name = 'partitions_fits.csv'
    df_partitions_fits_name = None
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'F:/Research/Results/Azimuth_Analysis/{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt')

    sim_sets = []

    stat_plot = 'k2'
    div_plt = 120
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 39
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [7, 11, 19, 27, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['diff']
    samples = 72  # For title purposes only

    # data_sets_plt = ['bes_def', 'ampt']
    data_sets_plt = ['bes_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    v2_star_vals = {2: read_flow_values(v2_star_in_dir)}
    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}
    v2_cf_vals = {2: read_flow_values(v2_cf_in_dir)}
    v2_cfev_vals = {2: read_flow_values(v2_cfev_in_dir)}
    v2_cfevb342_vals = {2: read_flow_values(v2_cfevb342_in_dir)}

    v2_vals = {'bes_def': v2_star_vals, 'ampt': v2_ampt_vals,
               'cf': v2_cf_vals, 'cfev': v2_cfev_vals, 'cfevb342': v2_cfevb342_vals}

    df = pd.read_csv(df_path)
    df = df.dropna()
    all_sets_plt = pd.unique(df['name'])
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df = df[df['stat'] == stat_plot]
    df_diff = calc_dsigma(df)

    dvar_vs_protons(df_diff, div_plt, 8, [39], ['diff'], data_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons_cents(df_diff, [div_plt], cents, energy_plt, data_types_plt, data_sets_plt, plot=True,
                          avg=True, plot_avg=True, data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels,
                          y_ranges=[-0.0051, 0.0013])

    dsig_avgs = []
    for energy in energies_fit:
        print(f'Energy {energy}GeV')
        dsig_avgs_div_diff = dvar_vs_protons_cents(df_diff, divs_all, cents, energy, ['diff'], all_sets_plt,
                                                   plot=False, avg=True)
        dsig_avgs_div_diff.loc[:, 'name'] = dsig_avgs_div_diff['name'] + '_flow_included'
        for data_set in all_sets_plt:
            dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, f'{data_set}_flow_included',
                                                      data_set, v2_vals[data_set])
        # dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'ampt_new_coal_epbins1_flow_included',
        #                                           'ampt_new_coal_epbins1', v2_ampt_vals)
        # dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'bes_resample_epbins1_flow_included',
        #                                           'bes_resample_epbins1', v2_star_vals)
        # dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cf_resample_epbins1_flow_included',
        #                                           'cf_resample_epbins1', v2_cf_vals)
        # dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cfev_resample_epbins1_flow_included',
        #                                           'cfev_resample_epbins1', v2_cfev_vals)
        # dsig_avgs_div_diff = subtract_dsigma_flow(dsig_avgs_div_diff, 'cfevb342_resample_epbins1_flow_included',
        #                                           'cfevb342_resample_epbins1', v2_cfevb342_vals)
        dsig_avgs.append(dsig_avgs_div_diff)

    dsig_avgs = pd.concat(dsig_avgs, ignore_index=True)
    if df_tproton_avgs_name is not None:
        dsig_avgs.to_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}', index=False)

    plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True, plot=True,
                        data_sets_labels=data_sets_labels, plot_energy_panels=True)

    df_fits = plot_dvar_avgs_divs(dsig_avgs, all_sets_plt, data_sets_colors=data_sets_colors, fit=True, plot=False,
                                  data_sets_labels=data_sets_labels, plot_energy_panels=False)
    if df_partitions_fits_name is not None:
        df_fits.to_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}', index=False)

    for energy in energies_fit:
        dsig_avgs_energy = dsig_avgs[dsig_avgs['energy'] == energy]
        plot_dvar_avgs_divs_cents(dsig_avgs_energy, data_sets_plt, plot=True, data_sets_colors=data_sets_colors,
                                  data_sets_labels=data_sets_labels, exclude_divs=exclude_divs, fit=False)
        plot_dvar_avgs_divs_cents(dsig_avgs_energy, data_sets_plt, plot=True, data_sets_colors=data_sets_colors,
                                  data_sets_labels=data_sets_labels, exclude_divs=exclude_divs, fit=True)

    df_onediv_avgs = dsig_avgs[dsig_avgs['divs'] == div_plt]
    plot_protons_avgs_vs_cent(df_onediv_avgs, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                              data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
                              title=f'{div_plt}° Partitions, {samples} Samples per Event',
                              data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels, cent_ref=cent_ref_df, ref_type=ref_type,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_protons_fits_vs_cent(df_onediv_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                           data_sets_labels=data_sets_labels,
    #                           title=f'{div_plt}° Partitions, {samples} Samples per Event', ls='none',
    #                           data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=None,
                          ref_type=None, data_sets_energies_cmaps=data_sets_energies_cmaps)
    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_vs_cent_var_fits():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_avgs_name = 'dsig_tprotons_avgs_old.csv'
    df_partitions_fits_name = 'partitions_fits_old.csv'

    df_fits = pd.read_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}')
    df_divs_avgs = pd.read_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}')

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'F:/Research/Results/Azimuth_Analysis/{cent_ref_name}')
    ref_type = 'ref'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df_bes_old = cent_ref_df.replace('bes_def', 'bes_resample_epbins1')
    cent_ref_df = pd.concat([cent_ref_df,
                             cent_ref_df_bes_old[cent_ref_df_bes_old['data_set'] == 'bes_resample_epbins1']],
                            ignore_index=True)
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')
    print(cent_ref_df)
    print(pd.unique(cent_ref_df['data_set']))

    data_sets_plt = ['bes_def', 'ampt_new_coal_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['tab10', 'Dark2']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    # data_sets_plt = ['bes_resample_epbins1', 'bes_def']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    # data_sets_energies_cmaps = dict(zip(data_sets_plt, ['tab10', 'Dark2']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR bad', 'STAR fixed']))

    energy = 7
    df_divs_avgs_energy = df_divs_avgs[df_divs_avgs['energy'] == energy]
    plot_dvar_avgs_vs_divs_cents(df_divs_avgs_energy, data_sets_plt, data_sets_colors=data_sets_colors,
                                 data_sets_labels=data_sets_labels, fit=True)
    div = 120
    df_div_avgs_div = df_divs_avgs[df_divs_avgs['divs'] == div]
    df_div_avgs_div = df_div_avgs_div.rename({'name': 'data_set', 'avg': 'baseline', 'avg_err': 'base_err'},
                                             axis='columns')
    # plot_div_fits_vs_cent_62res(df_div_avgs_div, ['bes_def'], data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, title=f'BES1 {div}', fit=True, cent_ref=cent_ref_df,
    #                             ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    print(df_div_avgs_div)
    df_div_avgs_div_energy = df_div_avgs_div[df_div_avgs_div['energy'] == energy]
    df_div_avgs_div_energy['zero_mag'] = 0
    df_div_avgs_div_energy['zero_mag_err'] = 0
    plot_div_fits_vs_cent(df_div_avgs_div_energy, ['bes_def'], data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=f'BES1 {div}° {energy} GeV', cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=None, fit=True)

    # plt.show()

    # plot_div_fits_vs_cent(df_fits, ['bes_def'], data_sets_colors=data_sets_colors,
    #                       data_sets_labels=data_sets_labels, title=f'BES1', fit=False, cent_ref=cent_ref_df,
    #                       ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_div_fits_vs_cent(df_fits, ['bes_def'], data_sets_colors=data_sets_colors,
    #                       data_sets_labels=data_sets_labels, title=f'BES1', fit=True, cent_ref=cent_ref_df,
    #                       ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    # plot_div_fits_vs_cent_62res(df_fits, ['ampt_new_coal_epbins1'], data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, title=f'AMPT', fit=True, cent_ref=cent_ref_df,
    #                             ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plot_div_fits_vs_cent_62res(df_fits, ['bes_def'], data_sets_colors=data_sets_colors,
    #                             data_sets_labels=data_sets_labels, title=f'BES1', fit=True, cent_ref=cent_ref_df,
    #                             ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)
    #
    # plt.show()
    #
    plot_div_fits_vs_cent_62res(df_fits, ['bes_def'], data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, title=f'BES1 ref', fit=True, cent_ref=cent_ref_df,
                                ref_type='ref', data_sets_energies_cmaps=data_sets_energies_cmaps)

    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=f'{energy} GeV', fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps, fit_boundary=0)
    # plt.show()

    energies = [7, 11, 19, 27, 39, 62]
    # energies = [62]
    for energy in energies:
        df_fits_energy = df_fits[df_fits['energy'] == energy]
        plot_div_fits_vs_cent(df_fits_energy, data_sets_plt, data_sets_colors=data_sets_colors,
                              data_sets_labels=data_sets_labels, title=f'{energy} GeV', fit=True, cent_ref=cent_ref_df,
                              ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps, fit_boundary=0)

    plt.show()


def plot_vs_cent_var_fit_tests():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_avgs_name = 'dsig_tprotons_avgs_old.csv'
    df_partitions_fits_name = 'partitions_fits_old.csv'

    df_fits = pd.read_csv(f'{base_path}{fits_out_base}/{df_partitions_fits_name}')
    df_divs_avgs = pd.read_csv(f'{base_path}{fits_out_base}/{df_tproton_avgs_name}')
    # print(pd.unique(df_divs_avgs['name']))

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'F:/Research/Results/Azimuth_Analysis/{cent_ref_name}')
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df_bes_old = cent_ref_df.replace('bes_def', 'bes_resample_epbins1')
    cent_ref_df = pd.concat([cent_ref_df,
                             cent_ref_df_bes_old[cent_ref_df_bes_old['data_set'] == 'bes_resample_epbins1']],
                            ignore_index=True)
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    def invx_test(x, a):
        return a / np.sqrt(x)

    def invx_const_test(x, a, b):
        return a / np.sqrt(x) + b

    def invx_sigmoid_test(x, a, b, c, d):
        return a / np.sqrt(x) + b / (1 + np.exp(-c * (x - d)))

    energy = 7
    div = 180
    data_set = 'bes_def'
    # data_set = 'ampt_new_coal_epbins1'
    df_div_energy_avgs = df_divs_avgs[(df_divs_avgs['name'] == data_set) &
                                      (df_divs_avgs['divs'] == div) & (df_divs_avgs['energy'] == energy)]
    cent_ref_df_energy = cent_ref_df[(cent_ref_df['data_set'] == 'bes_def') &
                                     (cent_ref_df['energy'] == energy)][['cent', 'mean_ref_val', 'mean_ref_sd']]
    df_test = pd.merge(df_div_energy_avgs, cent_ref_df_energy, on='cent')

    p01 = [-0.01]
    popt1, pcov1 = cf(invx_test, df_test['mean_ref_val'], df_test['avg'],
                      sigma=df_test['avg_err'], absolute_sigma=True, p0=p01)
    pfit1 = [Measure(val, err) for val, err in zip(popt1, np.sqrt(np.diag(pcov1)))]
    print(pfit1)

    p02 = [-0.01, 0.001]
    popt2, pcov2 = cf(invx_const_test, df_test['mean_ref_val'], df_test['avg'],
                      sigma=df_test['avg_err'], absolute_sigma=True, p0=p02)
    pfit2 = [Measure(val, err) for val, err in zip(popt2, np.sqrt(np.diag(pcov2)))]
    print(pfit2)

    p03 = [-0.02, 0.0004, 0.1, 20]
    popt3, pcov3 = cf(invx_sigmoid_test, df_test['mean_ref_val'], df_test['avg'],
                      sigma=df_test['avg_err'], absolute_sigma=True, p0=p03)
    pfit3 = [Measure(val, err) for val, err in zip(popt3, np.sqrt(np.diag(pcov3)))]
    print(pfit3)

    x_fit = np.linspace(1, 400, 1000)
    fig, ax = plt.subplots()
    ax.set_xlabel('RefMult')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$')
    ax.set_title(f'{data_set} {energy}GeV {div}° Partition Width')
    ax.grid()
    ax.axhline(0, color='gray')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'], yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'],
                ls='none', marker='o', color='black')
    ax.plot(x_fit, invx_test(x_fit, *popt1), color='gray')
    ax.plot(x_fit, invx_const_test(x_fit, *popt2), color='blue')
    ax.plot(x_fit, invx_sigmoid_test(x_fit, *popt3), color='orange')
    ax.plot(x_fit, invx_sigmoid_test(x_fit, *p03), color='orange', ls='--')
    ax.set_ylim(-0.0052, 0.0005)
    ax.set_xlim(-2, 400)

    fig, ax = plt.subplots()
    ax.set_xlabel('RefMult')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$ Fit Resituals')
    ax.set_title(f'{data_set} {energy}GeV {div}° Partition Width')
    ax.grid()
    ax.axhline(0, color='black')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_test(df_test['mean_ref_val'], *popt1),
                yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='gray',
                label=r'$a/\sqrt{M}$')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_const_test(df_test['mean_ref_val'], *popt2),
                yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='blue',
                label=r'$a/\sqrt{M}+b$')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_test(df_test['mean_ref_val'], popt2[0]),
                yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='salmon',
                label=r'$a/\sqrt{M}+b$ b fit')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_sigmoid_test(df_test['mean_ref_val'], *popt3),
                yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='orange',
                label=r'$a/\sqrt{M}+sigmoid$')
    ax.legend()
    fig.tight_layout()

    # plt.show()

    def invx_n_test(x, a, b, c):
        return a / x ** b + c

    # def invx_const_test(x, a, b):
    #     return a / np.sqrt(x) + b
    #
    # def invx_sigmoid_test(x, a, b, c, d):
    #     return a / np.sqrt(x) + b / (1 + np.exp(-c * (x - d)))

    p01 = [-0.01, 0.5, 0.0001]
    popt1, pcov1 = cf(invx_n_test, df_test['mean_ref_val'], df_test['avg'],
                      sigma=df_test['avg_err'], absolute_sigma=True, p0=p01)
    pfit1 = [Measure(val, err) for val, err in zip(popt1, np.sqrt(np.diag(pcov1)))]
    print(pfit1)

    # p02 = [-0.01, 0.001]
    # popt2, pcov2 = cf(invx_const_test, df_test['mean_ref_val'], df_test['avg'],
    #                   sigma=df_test['avg_err'], absolute_sigma=True, p0=p02)
    # pfit2 = [Measure(val, err) for val, err in zip(popt2, np.sqrt(np.diag(pcov2)))]
    # print(pfit2)
    #
    # p03 = [-0.02, 0.0004, 0.1, 20]
    # popt3, pcov3 = cf(invx_sigmoid_test, df_test['mean_ref_val'], df_test['avg'],
    #                   sigma=df_test['avg_err'], absolute_sigma=True, p0=p03)
    # pfit3 = [Measure(val, err) for val, err in zip(popt3, np.sqrt(np.diag(pcov3)))]
    # print(pfit3)

    fig, ax = plt.subplots()
    ax.set_xlabel('RefMult')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$')
    ax.set_title(f'{data_set} {energy}GeV {div}° Partition Width')
    ax.grid()
    ax.axhline(0, color='gray')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'], yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'],
                ls='none', marker='o', color='black')
    ax.plot(x_fit, invx_n_test(x_fit, *popt1), color='gray')
    # ax.plot(x_fit, invx_const_test(x_fit, *popt2), color='blue')
    # ax.plot(x_fit, invx_sigmoid_test(x_fit, *popt3), color='orange')
    # ax.plot(x_fit, invx_sigmoid_test(x_fit, *p03), color='orange', ls='--')
    ax.set_ylim(-0.0052, 0.0005)
    ax.set_xlim(-2, 400)

    fig, ax = plt.subplots()
    ax.set_xlabel('RefMult')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$ Fit Resituals')
    ax.set_title(f'{data_set} {energy}GeV {div}° Partition Width')
    ax.grid()
    ax.axhline(0, color='black')
    ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_n_test(df_test['mean_ref_val'], *popt1),
                yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='gray',
                label=r'$a/M^n$')
    # ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_const_test(df_test['mean_ref_val'], *popt2),
    #             yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='blue',
    #             label=r'$a/\sqrt{M}+b$')
    # ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_test(df_test['mean_ref_val'], popt2[0]),
    #             yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='salmon',
    #             label=r'$a/\sqrt{M}+b$ b fit')
    # ax.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_sigmoid_test(df_test['mean_ref_val'], *popt3),
    #             yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o', color='orange',
    #             label=r'$a/\sqrt{M}+sigmoid$')
    ax.legend()
    fig.tight_layout()

    # plt.show()

    energies = [7, 11, 19, 27, 39, 62]

    fig, ax = plt.subplots()
    ax.set_xlabel('RefMult')
    ax.set_ylabel(r'$\langle\Delta\sigma^2\rangle$')
    ax.set_title(f'{data_set} {div}° Partition Width')
    ax.grid()
    ax.axhline(0, color='gray')
    ax.set_ylim(-0.0052, 0.0005)
    ax.set_xlim(-2, 230)

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('RefMult')
    ax2.set_ylabel(r'$\langle\Delta\sigma^2\rangle$ Fit Resituals')
    ax2.set_title(f'{data_set} {div}° Partition Width')
    ax2.grid()
    ax2.axhline(0, color='black')

    fig3, ax3 = plt.subplots()
    ax3.set_xlabel('RefMult')
    ax3.set_ylabel(r'$\langle\Delta\sigma^2\rangle$ Fit - $1/\sqrt{M}$ Term')
    ax3.set_title(f'{data_set} {div}° Partition Width')
    ax3.grid()
    ax3.axhline(0, color='black')

    x_fit = np.linspace(1, 400, 1000)

    df_fit = []
    for energy in energies:
        df_div_energy_avgs = df_divs_avgs[(df_divs_avgs['name'] == data_set) &
                                          (df_divs_avgs['divs'] == div) & (df_divs_avgs['energy'] == energy)]
        cent_ref_df_energy = cent_ref_df[(cent_ref_df['data_set'] == 'bes_def') &
                                         (cent_ref_df['energy'] == energy)][['cent', 'mean_ref_val', 'mean_ref_sd']]
        df_test = pd.merge(df_div_energy_avgs, cent_ref_df_energy, on='cent')

        p02 = [-0.01, 0.001]
        popt2, pcov2 = cf(invx_const_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p02)
        pfit2 = [Measure(val, err) for val, err in zip(popt2, np.sqrt(np.diag(pcov2)))]
        print(f'{energy} GeV: ', pfit2)

        p03 = [p02[0] * 1.5, p02[1] * 0.9, 0.05, 25]
        popt3, pcov3 = cf(invx_sigmoid_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p03)
        pfit3 = [Measure(val, err) for val, err in zip(popt3, np.sqrt(np.diag(pcov3)))]
        print(f'{energy} GeV: ', pfit3)

        p04 = [-0.01, 0.5, 0.001]
        popt4, pcov4 = cf(invx_n_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p04)
        pfit4 = [Measure(val, err) for val, err in zip(popt4, np.sqrt(np.diag(pcov4)))]
        print(f'{energy} GeV: ', pfit4)
        df_fit.append({'energy': energy, 'a_val': pfit4[0].val, 'a_err': pfit4[0].err, 'b_val': pfit4[1].val,
                       'b_err': pfit4[1].err, 'c_val': pfit4[2].val, 'c_err': pfit4[2].err})

        ebar = ax.errorbar(df_test['mean_ref_val'], df_test['avg'], yerr=df_test['avg_err'],
                           xerr=df_test['mean_ref_sd'], ls='none', marker='o', label=f'{energy} GeV')
        ax.plot(x_fit, invx_const_test(x_fit, *popt2), color=ebar[0].get_color(), ls=':', alpha=0.4)
        ax.plot(x_fit, invx_sigmoid_test(x_fit, *popt3), color=ebar[0].get_color())

        ax2.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_test(df_test['mean_ref_val'], popt2[0]),
                     yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='^',
                     color=ebar[0].get_color(), label=f'{energy} GeV', alpha=0.4)
        ax2.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_sigmoid_test(df_test['mean_ref_val'], *popt3),
                     yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o',
                     color=ebar[0].get_color(), label=f'{energy} GeV')

        ax3.errorbar(df_test['mean_ref_val'], df_test['avg'] - invx_test(df_test['mean_ref_val'], popt3[0]),
                     yerr=df_test['avg_err'], xerr=df_test['mean_ref_sd'], ls='none', marker='o',
                     color=ebar[0].get_color(), label=f'{energy} GeV', alpha=0.5)
        ax3.plot(x_fit, invx_sigmoid_test(x_fit, *popt3) - invx_test(x_fit, popt3[0]), color=ebar[0].get_color(),
                 alpha=0.5)
    ax.legend()
    ax2.legend()
    ax3.legend()

    data_set = 'ampt_new_coal_epbins1'
    df_fit2 = []
    for energy in energies:
        df_div_energy_avgs = df_divs_avgs[(df_divs_avgs['name'] == data_set) &
                                          (df_divs_avgs['divs'] == div) & (df_divs_avgs['energy'] == energy)]
        cent_ref_df_energy = cent_ref_df[(cent_ref_df['data_set'] == 'bes_def') &
                                         (cent_ref_df['energy'] == energy)][['cent', 'mean_ref_val', 'mean_ref_sd']]
        df_test = pd.merge(df_div_energy_avgs, cent_ref_df_energy, on='cent')

        p02 = [-0.01, 0.001]
        popt2, pcov2 = cf(invx_const_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p02)
        pfit2 = [Measure(val, err) for val, err in zip(popt2, np.sqrt(np.diag(pcov2)))]
        print(f'{energy} GeV: ', pfit2)

        p03 = [p02[0] * 1.5, p02[1] * 0.9, 0.05, 25]
        popt3, pcov3 = cf(invx_sigmoid_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p03)
        pfit3 = [Measure(val, err) for val, err in zip(popt3, np.sqrt(np.diag(pcov3)))]
        print(f'{energy} GeV: ', pfit3)

        p04 = [-0.01, 0.5, 0.001]
        popt4, pcov4 = cf(invx_n_test, df_test['mean_ref_val'], df_test['avg'],
                          sigma=df_test['avg_err'], absolute_sigma=True, p0=p04)
        pfit4 = [Measure(val, err) for val, err in zip(popt4, np.sqrt(np.diag(pcov4)))]
        print(f'{energy} GeV: ', pfit4)
        df_fit2.append({'energy': energy + 1, 'a_val': pfit4[0].val, 'a_err': pfit4[0].err, 'b_val': pfit4[1].val,
                        'b_err': pfit4[1].err, 'c_val': pfit4[2].val, 'c_err': pfit4[2].err})

    df_fit = pd.DataFrame(df_fit)
    df_fit2 = pd.DataFrame(df_fit2)

    plt.rcParams["figure.figsize"] = (4.4, 4)

    fig4, ax4 = plt.subplots()
    ax4.set_xlabel('Energy')
    ax4.set_ylabel('a')
    ax4.set_title(f'{div}° Partition Width')
    ax4.errorbar(df_fit['energy'], df_fit['a_val'], df_fit['a_err'], ls='none', marker='o', label='STAR')
    ax4.errorbar(df_fit2['energy'], df_fit2['a_val'], df_fit2['a_err'], ls='none', marker='o', label='AMPT')
    ax4.grid()
    ax4.text(0.55, 0.9, r'$\frac{a}{M^n}+c$', fontsize=20, ha='center', va='top', transform=ax4.transAxes)
    ax4.legend()
    fig4.tight_layout()

    fig4, ax4 = plt.subplots()
    ax4.set_xlabel('Energy')
    ax4.set_ylabel('n')
    ax4.set_title(f'{div}° Partition Width')
    ax4.errorbar(df_fit['energy'], df_fit['b_val'], df_fit['b_err'], ls='none', marker='o', label='STAR')
    ax4.errorbar(df_fit2['energy'], df_fit2['b_val'], df_fit2['b_err'], ls='none', marker='o', label='AMPT')
    ax4.grid()
    ax4.text(0.55, 0.9, r'$\frac{a}{M^n}+c$', fontsize=20, ha='center', va='top', transform=ax4.transAxes)
    ax4.legend()
    fig4.tight_layout()

    fig4, ax4 = plt.subplots()
    ax4.set_xlabel('Energy')
    ax4.set_ylabel('c')
    ax4.axhline(0, color='black')
    ax4.set_title(f'{div}° Partition Width')
    ax4.errorbar(df_fit['energy'], df_fit['c_val'], df_fit['c_err'], ls='none', marker='o', label='STAR')
    ax4.errorbar(df_fit2['energy'], df_fit2['c_val'], df_fit2['c_err'], ls='none', marker='o', label='AMPT')
    ax4.grid()
    ax4.text(0.55, 0.9, r'$\frac{a}{M^n}+c$', fontsize=20, ha='center', va='top', transform=ax4.transAxes)
    ax4.legend()
    fig4.tight_layout()

    plt.show()


def plot_vs_cent_var_fit_tests2():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'

    dsig_avgs_v2sub = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_name}')
    dsig_avgs_v2sub_model = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_model_name}')
    dsig_avgs_v2sub = pd.concat([dsig_avgs_v2sub, dsig_avgs_v2sub_model])

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    # dsig_avgs_v2sub['data_type'] = 'diff'
    df = dsig_avgs_v2sub[dsig_avgs_v2sub['divs'] == 120]
    data_sets = ['bes_def', 'ampt_new_coal_epbins1']

    def n_fit(n, a, b):
        return a / n + b

    def n_pow_fit(n, a, b, c):
        return a / n ** c + b

    def n_fit_lin(n, a, b, c):
        return a / (n + c * np.sqrt(n)) + b

    def n_fit_lin2(n, a, b, c, d):
        return a / np.array(n) + c / np.power(n, d) + b

    def n_n2_fit(n, a, b, c):
        return a / n + c / np.power(n, 2) + b

    fit_func = n_fit_lin
    p0 = (-0.2, 0, 0)
    min_x_fit = 0
    cents_fit = 3
    fit_resid = False
    n_pars = len(inspect.signature(fit_func).parameters) - 1

    energies = [7, 11, 19, 27, 39, 62]
    energy_fit_pars = {i: {'bes': {'vals': [], 'errs': []}, 'ampt': {'vals': [], 'errs': []}} for i in range(n_pars)}
    for energy in energies:
        df_energy = df[df['energy'] == energy]
        fig, ax = plt.subplots()
        ax.grid()
        ax.axhline(0, color='gray')
        ax.set_title(f'{energy} GeV')
        ax.set_xlabel('Refmult')
        ax.set_ylabel(r'$\langle \Delta \sigma^2 \rangle$')

        for data_set in data_sets:
            df_set = df_energy[df_energy['name'] == data_set]
            cent_energy = cent_ref_df[(cent_ref_df['data_set'] == data_set) & (cent_ref_df['energy'] == energy)]
            x = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_val'].iloc[0] for cent in
                 df_set['cent']]
            x_err = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_sd'].iloc[0] for cent in
                     df_set['cent']]

            ax.errorbar(x, df_set['avg'], yerr=df_set['avg_err'], xerr=x_err, ls='none', marker='o',
                        label=data_set)

            filtered_lists = [(xi, y, ye) for xi, y, ye in zip(x, df_set['avg'], df_set['avg_err'])
                              if xi > min_x_fit]
            x_fit, y_fit, yerr_fit = zip(*sorted(filtered_lists))
            if cents_fit is not None:
                x_fit, y_fit, yerr_fit = x_fit[-cents_fit:], y_fit[-cents_fit:], yerr_fit[-cents_fit:]

            try:
                popt, pcov = cf(fit_func, x_fit, y_fit, sigma=yerr_fit, absolute_sigma=True, p0=p0)
            except RuntimeError:
                popt, pcov = p0, np.diag(np.ones(len(p0)))
            perr = np.sqrt(np.diag(pcov))
            xs = np.linspace(min(x), max(x), 1000)
            ax.plot(xs, fit_func(xs, *popt), color='red')
            if 'ampt' in data_set:
                ampt_popt = popt
                x_ampt = x
                for i in range(n_pars):
                    energy_fit_pars[i]['ampt']['vals'].append(popt[i])
                    energy_fit_pars[i]['ampt']['errs'].append(perr[i])
            if 'bes_def' in data_set:
                bes_popt = popt
                x_bes = x
                for i in range(n_pars):
                    energy_fit_pars[i]['bes']['vals'].append(popt[i])
                    energy_fit_pars[i]['bes']['errs'].append(perr[i])
        ax.legend()
        fig.canvas.manager.set_window_title(f'{energy} GeV Fits')
        fig.tight_layout()

        fig2, ax2 = plt.subplots()
        ax2.grid()
        ax2.axhline(0, color='gray')
        ax2.set_title(f'{energy} GeV AMPT Subtracted')
        ax2.set_xlabel('Refmult')
        ax2.set_ylabel(r'$\langle \Delta \sigma^2 \rangle$')

        df_bes = df_energy[df_energy['name'] == 'bes_def']
        bes_sub_ampt = df_bes['avg'] - fit_func(x_bes, *ampt_popt)
        ax2.errorbar(x_bes, bes_sub_ampt, yerr=df_bes['avg_err'], ls='none', marker='o', label='bes - ampt')

        df_ampt = df_energy[df_energy['name'] == 'ampt_new_coal_epbins1']
        ampt_sub_ampt = df_ampt['avg'] - fit_func(x_ampt, *ampt_popt)
        ax2.errorbar(x_ampt, ampt_sub_ampt, yerr=df_ampt['avg_err'], ls='none', marker='o', label='ampt - ampt')
        if fit_resid:
            popt_ampt_resid, pcov_ampt_resid = cf(n_pow_fit, x_ampt, ampt_sub_ampt, sigma=df_ampt['avg_err'],
                                                  absolute_sigma=True)
            xs = np.linspace(min(x_ampt), max(x_ampt), 1000)
            print(f'AMPT {energy} GeV: {popt_ampt_resid}')
            ax2.plot(xs, n_pow_fit(xs, *popt_ampt_resid))

        bes_sub_bes = df_bes['avg'] - fit_func(x_bes, *bes_popt)
        ax2.errorbar(x_bes, bes_sub_bes, yerr=df_bes['avg_err'], ls='none', marker='o', label='bes - bes')
        if fit_resid:
            popt_bes_resid, pcov_bes_resid = cf(n_pow_fit, x_bes, bes_sub_bes, sigma=df_bes['avg_err'],
                                                absolute_sigma=True)
            xs = np.linspace(min(x_bes), max(x_bes), 1000)
            print(f'BES {energy} GeV: {popt_bes_resid}')
            ax2.plot(xs, n_pow_fit(xs, *popt_bes_resid))

        ax2.legend()
        fig2.canvas.manager.set_window_title(f'{energy} GeV Residuals')
        fig2.tight_layout()

    for par_i, par_letter in enumerate(['a', 'b', 'c', 'd', 'e'][:len(popt)]):
        fit_pars = energy_fit_pars[par_i]
        fig3, ax3 = plt.subplots()
        ax3.grid()
        ax3.axhline(0, color='gray')
        ax3.errorbar(energies, fit_pars['bes']['vals'], yerr=fit_pars['bes']['errs'],
                     marker='o', ls='none', label='bes')
        ax3.errorbar(energies, fit_pars['ampt']['vals'], yerr=fit_pars['ampt']['errs'],
                     marker='o', ls='none', label='ampt')
        ax3.set_xlabel('Energy')
        ax3.set_ylabel(par_letter)
        ax3.legend()
        fig3.canvas.manager.set_window_title(f'Parameter {par_letter} vs Energy')
        fig3.tight_layout()

    plt.show()


def plot_vs_cent_var_fit_tests3():
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'

    dsig_avgs_v2sub = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_name}')
    dsig_avgs_v2sub_model = pd.read_csv(f'{base_path}{df_def_avgs_v2sub_out_model_name}')
    dsig_avgs_v2sub = pd.concat([dsig_avgs_v2sub, dsig_avgs_v2sub_model])

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    # dsig_avgs_v2sub['data_type'] = 'diff'
    df = dsig_avgs_v2sub[dsig_avgs_v2sub['divs'] == 120]
    data_sets = ['bes_def', 'ampt_new_coal_epbins1']

    plt.rcParams.update({'font.size': 16})

    def n_fit(n, a, b):
        return a / n + b

    def n_pow_fit(n, a, b, c):
        return a * (1 + b / np.power(n, c))

    def n_pow_fit2(n, a, b, c):
        return a + b / np.power(n, c)

    def n_fit_lin(n, a, b, c):
        return a / (n + c * np.sqrt(n)) + b

    def n_fit_lin2(n, a, b, c, d):
        return a / np.array(n) + c / np.power(n, d) + b

    def n_n2_fit(n, a, b, c):
        return a / n + c / np.power(n, 2) + b

    fit_func = n_pow_fit
    fit_eq = r'$\langle \Delta \sigma^2 \rangle = a \left( 1 + \frac{b}{M^c} \right)$'
    # fit_eq = r'$\langle \Delta \sigma^2 \rangle = a + \frac{b}{M^c}$'
    p0 = (0.0002, -100, 0.5)
    min_x_fit = 0
    # cents_fit = 3
    cents_fit = None
    fit_resid = False
    n_pars = len(inspect.signature(fit_func).parameters) - 1

    energies = [7, 11, 19, 27, 39, 62]
    energy_fit_pars = {i: {'bes': {'vals': [], 'errs': []}, 'ampt': {'vals': [], 'errs': []}} for i in range(n_pars)}
    for energy in energies:
        df_energy = df[df['energy'] == energy]
        fig, ax = plt.subplots()
        ax.grid()
        ax.axhline(0, color='gray')
        ax.set_title(f'{energy} GeV')
        ax.set_xlabel('Refmult')
        ax.set_ylabel(r'$\langle \Delta \sigma^2 \rangle$')

        for data_set in data_sets:
            print(f'Starting {data_set} {energy} GeV')
            df_set = df_energy[df_energy['name'] == data_set]
            cent_energy = cent_ref_df[(cent_ref_df['data_set'] == data_set) & (cent_ref_df['energy'] == energy)]
            x = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_val'].iloc[0] for cent in
                 df_set['cent']]
            x_err = [cent_energy[cent_energy['cent'] == cent][f'mean_{ref_type}_sd'].iloc[0] for cent in
                     df_set['cent']]

            ax.errorbar(x, df_set['avg'], yerr=df_set['avg_err'], xerr=x_err, ls='none', marker='o',
                        label=data_set)

            filtered_lists = [(xi, y, ye) for xi, y, ye in zip(x, df_set['avg'], df_set['avg_err'])
                              if xi > min_x_fit]
            x_fit, y_fit, yerr_fit = zip(*sorted(filtered_lists))
            if cents_fit is not None:
                x_fit, y_fit, yerr_fit = x_fit[-cents_fit:], y_fit[-cents_fit:], yerr_fit[-cents_fit:]

            def objective(params):
                y_fit_obj = fit_func(x_fit, *params)
                return np.sum(np.square((y_fit_obj - y_fit) / yerr_fit))

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    popt, pcov = cf(fit_func, x_fit, y_fit, sigma=yerr_fit, absolute_sigma=True, p0=p0)
                    result = basinhopping(objective, p0, niter=100, T=1.0, stepsize=0.5)
                    p0_basin = result.x
                    print(f'CF popt: {popt}')
                    print(f'BH popt: {p0_basin}')
                    popt_new, pcov_new = cf(fit_func, x_fit, y_fit, sigma=yerr_fit, absolute_sigma=True, p0=p0_basin)
                    popt = popt_new
                    pcov = pcov_new
            except RuntimeError:
                popt, pcov = p0, np.diag(np.ones(len(p0)))
            perr = np.sqrt(np.diag(pcov))
            xs = np.linspace(min(x), max(x), 1000)
            ax.plot(xs, fit_func(xs, *popt), color='red')
            if 'ampt' in data_set:
                ampt_popt = popt
                x_ampt = x
                for i in range(n_pars):
                    energy_fit_pars[i]['ampt']['vals'].append(popt[i])
                    energy_fit_pars[i]['ampt']['errs'].append(perr[i])
            if 'bes_def' in data_set:
                bes_popt = popt
                x_bes = x
                for i in range(n_pars):
                    energy_fit_pars[i]['bes']['vals'].append(popt[i])
                    energy_fit_pars[i]['bes']['errs'].append(perr[i])
        ax.legend()
        fig.canvas.manager.set_window_title(f'{energy} GeV Fits')
        fig.tight_layout()

        fig2, ax2 = plt.subplots()
        ax2.grid()
        ax2.axhline(0, color='gray')
        ax2.set_title(f'{energy} GeV AMPT Subtracted')
        ax2.set_xlabel('Refmult')
        ax2.set_ylabel(r'$\langle \Delta \sigma^2 \rangle$')

        df_bes = df_energy[df_energy['name'] == 'bes_def']
        bes_sub_ampt = df_bes['avg'] - fit_func(x_bes, *ampt_popt)
        ax2.errorbar(x_bes, bes_sub_ampt, yerr=df_bes['avg_err'], ls='none', marker='o', label='bes - ampt')

        df_ampt = df_energy[df_energy['name'] == 'ampt_new_coal_epbins1']
        ampt_sub_ampt = df_ampt['avg'] - fit_func(x_ampt, *ampt_popt)
        ax2.errorbar(x_ampt, ampt_sub_ampt, yerr=df_ampt['avg_err'], ls='none', marker='o', label='ampt - ampt')
        if fit_resid:
            popt_ampt_resid, pcov_ampt_resid = cf(n_pow_fit, x_ampt, ampt_sub_ampt, sigma=df_ampt['avg_err'],
                                                  absolute_sigma=True)
            xs = np.linspace(min(x_ampt), max(x_ampt), 1000)
            print(f'AMPT {energy} GeV: {popt_ampt_resid}')
            ax2.plot(xs, n_pow_fit(xs, *popt_ampt_resid))

        bes_sub_bes = df_bes['avg'] - fit_func(x_bes, *bes_popt)
        ax2.errorbar(x_bes, bes_sub_bes, yerr=df_bes['avg_err'], ls='none', marker='o', label='bes - bes')
        if fit_resid:
            popt_bes_resid, pcov_bes_resid = cf(n_pow_fit, x_bes, bes_sub_bes, sigma=df_bes['avg_err'],
                                                absolute_sigma=True)
            xs = np.linspace(min(x_bes), max(x_bes), 1000)
            print(f'BES {energy} GeV: {popt_bes_resid}')
            ax2.plot(xs, n_pow_fit(xs, *popt_bes_resid))

        ax2.legend()
        fig2.canvas.manager.set_window_title(f'{energy} GeV Residuals')
        fig2.tight_layout()

    fig_fit_pars, ax_fit_pars = plt.subplots(nrows=len(popt), figsize=(10, 8), sharex='all')
    for par_i, par_letter in enumerate(['a', 'b', 'c', 'd', 'e'][:len(popt)]):
        fit_pars = energy_fit_pars[par_i]
        fig3, ax3 = plt.subplots()
        ax3.grid()
        ax3.axhline(0, color='gray')
        ax3.errorbar(energies, fit_pars['bes']['vals'], yerr=fit_pars['bes']['errs'],
                     marker='o', ls='none', label='bes')
        ax3.errorbar(energies, fit_pars['ampt']['vals'], yerr=fit_pars['ampt']['errs'],
                     marker='o', ls='none', label='ampt')
        ax3.set_xlabel('Energy')
        ax3.set_ylabel(par_letter)
        ax3.legend()
        fig3.canvas.manager.set_window_title(f'Parameter {par_letter} vs Energy')
        fig3.tight_layout()

        ax_fit_pars[par_i].axhline(0, color='gray')
        ax_fit_pars[par_i].errorbar(energies, fit_pars['bes']['vals'], yerr=fit_pars['bes']['errs'],
                                    marker='o', ls='none', label='bes', markersize=10, alpha=0.7)
        ax_fit_pars[par_i].errorbar(energies, fit_pars['ampt']['vals'], yerr=fit_pars['ampt']['errs'],
                                    marker='o', ls='none', label='ampt', markersize=10, alpha=0.7)
        ax_fit_pars[par_i].set_ylabel(par_letter)
        if par_i == 0:
            ax_fit_pars[par_i].legend()
        if par_i == len(popt) - 1:
            ax_fit_pars[par_i].set_xlabel('Energy')
            ax_fit_pars[par_i].annotate(fit_eq, (0.1, 0.2), xycoords='axes fraction', ha='left', va='bottom',
                                        fontsize=30)
        else:
            ax_fit_pars[par_i].tick_params(bottom=False, top=False)

    fig_fit_pars.subplots_adjust(top=0.99, bottom=0.075, left=0.115, right=0.99, hspace=0.05, wspace=0.05)

    plt.show()


def plot_vs_cent_nofit():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 60
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 62
    energies_fit = [7, 11, 19, 27, 39, 62]
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_divs_fits = pd.DataFrame()
    for energy in energies_fit:
        df_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy, data_types_plt, all_sets_plt,
                                       plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                       data_sets_labels=data_sets_labels)
        df_fits = plot_protons_fits_divs_cents(df_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                               data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)
        df_divs_fits = pd.concat([df_divs_fits, df_fits])

    plot_div_fits_vs_cent(df_divs_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=None, fit=False, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_vs_cent_fittest():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_cents.csv'
    df_path = base_path + df_name

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'ref'  # 'refn'

    print(cent_ref_df)

    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 60
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    energy_plt = 62
    # energies_fit = [energy_plt]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red']))
    data_sets_energies_cmaps = dict(zip(data_sets_plt, ['Greys', 'Reds']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_fit = stat_vs_protons_cents(df, stat_plot, divs_all, cents, energy_plt, data_types_plt, all_sets_plt,
                                   plot=False, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                                   data_sets_labels=data_sets_labels)
    df_fits = plot_protons_fits_divs_cents(df_fit, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                           data_sets_labels=data_sets_labels, exclude_divs=exclude_divs)

    plot_div_fits_vs_cent(df_fits, data_sets_plt, data_sets_colors=data_sets_colors,
                          data_sets_labels=data_sets_labels, title=f'{energy_plt}GeV', fit=True, cent_ref=cent_ref_df,
                          ref_type=ref_type, data_sets_energies_cmaps=data_sets_energies_cmaps)

    plt.show()


def plot_all_zero_base():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Base_Zero_Fits/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_names = ['partitions_fits.csv', 'sim_partitions_fits.csv']  #
    df_paths = [base_path + df_name for df_name in df_names]

    data_sets_plt = ['bes_resample_epbins1', 'ampt_new_coal_epbins1', 'cfev_resample_epbins1',
                     'cfevb342_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'CFEVb342']))

    df_fits = pd.DataFrame()
    for df_path in df_paths:
        df = pd.read_csv(df_path)
        df = df.dropna()
        df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['data_set'] else row['energy'], axis=1)
        print(f'{df_path}\n{df.head()}\n\n')
        df_fits = pd.concat([df_fits, df])

    print(pd.unique(df_fits[~df_fits['data_set'].str.contains('sim_')]['data_set']))

    plot_base_zeros(df_fits, data_sets_plt, data_sets_labels, data_sets_colors, plot_sims=False)
    plot_base_zeros(df_fits, data_sets_plt, data_sets_labels, data_sets_colors, plot_sims=True)
    plot_base_zeros(df_fits, data_sets_plt, data_sets_labels, data_sets_colors, plot_sims=True, cent=8)

    plt.show()


def plot_closest_sims():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Base_Zero_Fits/'
    # base_path = '/home/dylan/Research/Results/Azimuth_Analysis/'
    df_data_name = 'bes_ampt_cents_tprotons_fits.csv'
    df_sim_name = 'sim_tprotons_fits.csv'

    stat = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    divs_all = [60, 72, 89, 90, 180, 240, 270, 288, 300]
    cent = 8
    energy = 62
    num_sims_plt = 4
    data_name = 'ampt_new_coal_resample_def'  # 'bes_resample_def'
    data_title = 'AMPT'  # 'BES'

    # data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation']))

    data_sets_plt = ['bes_resample_def', 'ampt_new_coal_resample_def', 'cfev_resample_def', 'sim_aclmul_',
                     'sim_clmultipm_']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'purple', 'blue', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'CFEV', 'Simulation', 'Simulation 2 Gaus']))

    df_data = pd.read_csv(base_path + df_data_name)
    df_data = df_data.dropna()
    df_data = df_data[(df_data['name'] == data_name) & (df_data['energy'] == energy) & (df_data['cent'] == cent)]
    print(df_data)

    df_sim = pd.read_csv(base_path + df_sim_name)
    df_sim = df_sim.dropna()
    # df_sim['energy'] = df_sim.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_diffs = []
    for sim_set in pd.unique(df_sim['name']):
        df_sim_set = df_sim[df_sim['name'] == sim_set]
        diff, diff_weight = 0, 0
        for div in divs_all:
            diff += abs(df_data[df_data['divs'] == div]['slope'].iloc[0] -
                        df_sim_set[df_sim_set['divs'] == div]['slope'].iloc[0])
            weight = 10 if div == 180 else 1
            diff_weight += abs(df_data[df_data['divs'] == div]['slope'].iloc[0] -
                               df_sim_set[df_sim_set['divs'] == div]['slope'].iloc[0]) * weight
        diff_180 = abs(df_data[df_data['divs'] == 180]['slope'].iloc[0] -
                       df_sim_set[df_sim_set['divs'] == 180]['slope'].iloc[0])
        df_diffs.append({'sim_set': sim_set, 'diff': diff, 'diff_weight': diff_weight, 'diff_180': diff_180})
    df_diffs = pd.DataFrame(df_diffs)

    close_sims_diff = list(df_diffs.sort_values(by=['diff'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_diff = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_diff)]])

    close_sims_180 = list(df_diffs.sort_values(by=['diff_180'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_180 = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_180)]])

    close_sims_weight = list(df_diffs.sort_values(by=['diff_weight'], ascending=True)['sim_set'][:num_sims_plt])
    df_plt_weight = pd.concat([df_data, df_sim[df_sim['name'].isin(close_sims_weight)]])

    df_plt_amp1 = pd.concat([df_data, df_sim[(df_sim['amp'] == 0.1) &
                                             (df_sim['spread'] < 3) & (df_sim['spread'] > 0.5)]])
    df_plt_amp1 = df_plt_amp1.sort_values(by=['spread'])
    df_plt_spread1 = pd.concat([df_data, df_sim[(df_sim['spread'] == 1) & (df_sim['amp'] < 0.125)]])
    df_plt_spread1 = df_plt_spread1.sort_values(by=['amp'])

    plot_vs_div_width_comp(df_plt_diff, f'{data_title} {energy} GeV Equal Weight')
    plot_vs_div_width_comp(df_plt_180, f'{data_title} {energy} GeV 180 Only')
    plot_vs_div_width_comp(df_plt_weight, f'{data_title} {energy} GeV 180 Weighted x10')
    plot_vs_div_width_comp(df_plt_amp1, f'{data_title} {energy} GeV A=0.1')
    plot_vs_div_width_comp(df_plt_spread1, f'{data_title} {energy} GeV spread=1')

    plt.show()


def plot_ampt_efficiency():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_ampt_eff.csv'
    save_fits = False
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = 'ampt_eff_tprotons_fits.csv'
    df_partitions_fits_name = 'ampt_eff_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_resample_def', 'ampt_new_coal_resample_eff1', 'ampt_new_coal_resample_eff2',
                     'ampt_new_coal_resample_eff3', 'ampt_new_coal_resample_eff4']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['100%', '90%', '80%', '70%', '60%']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [39], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=False, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=False)

    # plt.show()
    # return

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    if save_fits:
        df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    def_set = 'ampt_new_coal_resample_def'
    df_fits_ratio = []  # Divide baseline/zeros by 100% default to see difference
    for energy in pd.unique(df_fits['energy']):
        df_fits_energy = df_fits[df_fits['energy'] == energy]
        df_fits_energy_def = df_fits_energy[df_fits_energy['data_set'] == def_set]
        base_def = Measure(df_fits_energy_def['baseline'].iloc[0], df_fits_energy_def['base_err'].iloc[0])
        zero_def = Measure(df_fits_energy_def['zero_mag'].iloc[0], df_fits_energy_def['zero_mag_err'].iloc[0])
        for ds_i, data_set in enumerate(pd.unique(df_fits_energy['data_set'])):
            if data_set != def_set:
                df_fits_e_ds = df_fits_energy[df_fits_energy['data_set'] == data_set]
                base_div = Measure(df_fits_e_ds['baseline'].iloc[0], df_fits_e_ds['base_err'].iloc[0]) / base_def
                zero_div = Measure(df_fits_e_ds['zero_mag'].iloc[0], df_fits_e_ds['zero_mag_err'].iloc[0]) / zero_def
                energy_shifted = energy - 1 + ds_i / 2.0
                df_fits_ratio.append({'data_set': data_set, 'energy': energy_shifted, 'baseline': base_div.val,
                                      'base_err': base_div.err, 'zero_mag': zero_div.val, 'zero_mag_err': zero_div.err})
    plot_slope_div_fits(pd.DataFrame(df_fits_ratio), data_sets_colors, data_sets_labels, ref_line=1)

    plt.show()


def plot_ampt_efficiency_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    v2_ampt_in_dir = 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/Ampt_rapid05_resample_norotate_epbins1_0/'
    df_name = 'binom_slice_vars_ampt_eff.csv'
    save_fits = False
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = 'ampt_eff_tprotons_fits.csv'
    df_partitions_fits_name = 'ampt_eff_partitions_fits.csv'
    df_path = base_path + df_name

    threads = 10

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    divs_all = [60, 72, 89, 90, 120, 180, 240, 270, 288, 300]
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cents = [1, 2, 3, 4, 5, 6, 7, 8]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['diff']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_resample_def', 'ampt_new_coal_resample_eff1', 'ampt_new_coal_resample_eff2',
                     'ampt_new_coal_resample_eff3', 'ampt_new_coal_resample_eff4']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green']))
    data_sets_labels = dict(zip(data_sets_plt, ['100%', '90%', '80%', '70%', '60%']))

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    df = pd.read_csv(df_path)
    df = df.dropna()
    df = df[df['stat'] == stat_plot]
    print(pd.unique(df['name']))

    v2_ampt_vals = {2: read_flow_values(v2_ampt_in_dir)}

    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])
    df_diff['meas'] = df_diff.apply(lambda row: Measure(row['val'], row['err']), axis=1)

    df_def_dsigma_v2sub = []
    for set_name in data_sets_plt:
        set_v2sub = subtract_dsigma_flow(df_diff, set_name, set_name, v2_ampt_vals,
                                         new_only=True, val_col='val', err_col='err', meas_col='meas')
        df_def_dsigma_v2sub.append(set_v2sub)
    df_dsigma_v2sub = pd.concat(df_def_dsigma_v2sub)

    dvar_vs_protons_energies(df_dsigma_v2sub, [div_plt], cent_plt, energies_fit, data_types_plt,
                             list(reversed(data_sets_plt)), plot=True, avg=True, plot_avg=True, alpha=0.6,
                             data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels,
                             ylabel=r'$\Delta\sigma^2$', kin_loc=(0.62, 0.4), y_ranges=[-0.00068, 0.00012])

    dsig_avgs_v2sub = []
    for energy in energies_fit:
        print(f'Energy {energy}GeV')
        jobs = [(df_diff[df_diff['name'] == set_i], divs_all, cents, energy, ['diff'],
                 [set_i], None, False, True) for set_i in data_sets_plt]
        dsig_avgs = []
        with Pool(threads) as pool:
            for job_i in tqdm.tqdm(pool.istarmap(dvar_vs_protons_cents, jobs), total=len(jobs)):
                dsig_avgs.append(job_i)
        dsig_avgs = pd.concat(dsig_avgs, ignore_index=True).drop('data_type', axis=1)
        for data_set in data_sets_plt:
            dsig_avgs_set = subtract_dsigma_flow(dsig_avgs, data_set, data_set, v2_ampt_vals, new_only=True)
            dsig_avgs_v2sub.append(dsig_avgs_set)

    dsig_avgs_v2sub = pd.concat(dsig_avgs_v2sub, ignore_index=True)

    dsig_avgs_v2sub_cent8 = dsig_avgs_v2sub[dsig_avgs_v2sub['cent'] == cent_plt]
    df_fits = plot_dvar_avgs_divs(dsig_avgs_v2sub_cent8, data_sets_plt, fit=True, data_sets_colors=data_sets_colors,
                                  data_sets_labels=data_sets_labels, plot_energy_panels=True,
                                  ylab=r'$\langle\Delta\sigma^2\rangle$', plot_indiv=False)

    plt.rcParams["figure.figsize"] = (8, 4)

    dsig_avgs_v2sub_div120_cent8 = dsig_avgs_v2sub[(dsig_avgs_v2sub['cent'] == cent_plt) &
                                                   (dsig_avgs_v2sub['divs'] == div_plt)]
    plot_protons_avgs_vs_energy(dsig_avgs_v2sub_div120_cent8, data_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels, alpha=1,
                                title=f'{cent_map[8]} Centrality, {div_plt}° Partitions, {samples} Samples per Event')

    plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    plt.show()


def plot_flow():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_stats_flow.csv'
    save_fits = False
    v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    v2_fit_out_dir = None
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = 'flow_tprotons_fits.csv'
    df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    # data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res75_v207', 'flow_resample_res5_v207',
    #                  'flow_resample_res3_v207', 'flow_resample_res15_v207']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green']))
    # data_sets_labels = dict(zip(data_sets_plt, ['0.99', '0.75', '0.5', '0.3', '0.15']))

    data_sets_plt = ['flow_resample_res15_v207', 'flow_resample_res15_v206', 'flow_resample_res15_v205',
                     'flow_resample_res15_v204', 'flow_resample_res15_v203', 'flow_resample_res15_v202']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple', 'green', 'olive']))
    data_sets_labels = dict(zip(data_sets_plt, ['v2=0.07', 'v2=0.06', 'v2=0.05', 'v2=0.04', 'v2=0.03', 'v2=0.02']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    df_raw_div = df[df['data_type'] == 'raw']
    # for data_set in data_sets_plt:
    #     df_set = df[df['name'] == data_set]
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], ['raw', 'mix'], all_sets_plt, plot=True, fit=False,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False)

    # flow_vs_v2(df, div_plt, '15', v2_fit_out_dir)

    # plt.show()
    # return

    # for data_set in data_sets_plt:
    #     df_set = df[df['name'] == data_set]
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
    #                     y_ranges={'standard deviation': [0.946, 1.045]})
    #     stat_vs_protons(df_set, stat_plot, div_plt, cent_plt, [62], ['raw'], all_sets_plt, plot=True, fit=True,
    #                     data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=False,
    #                     y_ranges={'standard deviation': [0.946, 1.045]})
    #
    #
    # plt.show()
    # return
    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        if v2_fit_out_dir:
            flow_vs_v2(df, div, '15', v2_fit_out_dir)
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    # return
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_protons_fits_divs_flow(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors,
                                data_sets_labels=data_sets_labels)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                                  data_sets_labels=data_sets_labels)
    # if save_fits:
    #     df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    plt.show()


def plot_flow_k2():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_v2_ck2.csv'
    save_fits = False
    v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    v2_fit_out_dir = None
    fits_out_base = 'Base_Zero_Fits/'
    df_tproton_fits_name = None  # 'flow_tprotons_fits.csv'
    df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_resample_res15_v207', 'flow_resample_res15_v205', 'flow_resample_res15_v202']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue']))
    data_sets_labels = dict(zip(data_sets_plt, ['v2=0.07', 'v2=0.05', 'v2=0.02']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)
    df = df[(df['data_type'] == 'raw') & ((df['stat'] == 'k2') | (df['stat'] == 'c2'))]
    df['val'] = df['val'] / (df['total_protons'] * df['divs'] / 360 * (1 - df['divs'] / 360))
    df['err'] = df['err'] / (df['total_protons'] * df['divs'] / 360 * (1 - df['divs'] / 360))

    stat_vs_protons(df, stat_plot, 120, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                    plot=True, fit=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        # if v2_fit_out_dir:
        #     flow_vs_v2(df, div, '15', v2_fit_out_dir)
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        protons_fits.append(protons_fits_div)
    # return
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if save_fits:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_protons_fits_divs_flow(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=False,
    #                                  data_sets_labels=data_sets_labels)
    # if save_fits:
    #     df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)

    plt.show()


def plot_ampt_v2_closure():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_ampt_v2_closure.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
    v2_in_dir = f'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/' \
                f'Ampt_rapid05_resample_norotate_epbins1_0/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = None  # 'cf_tprotons_fits.csv'
    df_partitions_fits_name = None  # 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    pd.set_option('max_columns', None)

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['ampt_new_coal_rp', 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                     'ampt_new_coal_epbins1_v2rpcor']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['Reaction Plane', 'Flow Included', 'Flow Corrected ep',
                                                'Flow Corrected rp']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}
    v2_vals = read_flow_values(v2_in_dir)  # Dictionary of {energy: {cent: v2}}
    v2_rp_vals = read_flow_values(v2_in_dir, 'v2_rp')  # Dictionary of {energy: {cent: v2}}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    # df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))  # Placeholder, need k2
    df_raw.loc[:, 'val'] = (df_raw['val'] ** 2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    # df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'val'] = (df_mix['val'] ** 2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))
    df_diff = df[(df['data_type'] == 'diff') & (df['stat'] == stat_plot)]
    p, tp = df_diff['divs'] / 360, df_diff['total_protons']
    df_diff.loc[:, 'val'] = df_diff['val'] / (tp * (tp - 1))
    df_diff.loc[:, 'err'] = df_diff['err'] / (tp * (tp - 1))

    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    dvar_vs_protons(df_sub, div_plt, cent_plt, [39], ['sub'], all_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)
    dvar_vs_protons_energies(df_sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['sub'], all_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.0014, 0.0011])

    protons_fits = []
    v2_apparents = {energy: [] for energy in energies_fit}
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets_plt, plot=False,
                                           avg=True)
        if div != 180:
            for energy in energies_fit:
                loc_ep = (protons_fits_div['name'] == 'ampt_new_coal_epbins1') & (protons_fits_div['energy'] == energy)
                loc_rp = (protons_fits_div['name'] == 'ampt_new_coal_rp') & (protons_fits_div['energy'] == energy)
                d_avg = protons_fits_div[loc_ep]['avg_meas'].iloc[0] - protons_fits_div[loc_rp]['avg_meas'].iloc[0]
                v2_apparent = np.sqrt(d_avg / flow_cor_coefs[div])
                v2_apparents[energy].append(v2_apparent)
                # print(f'{div}, {energy},  v2: {v2_apparent}')
        protons_fits_div = ampt_v2_closure_sub_dsigma(protons_fits_div, 'ampt_new_coal_epbins1',
                                                      'ampt_new_coal_epbins1_v2rpcor', v2_rp_vals,
                                                      flow_cor_coefs[div], cent_plt)
        protons_fits_div = ampt_v2_closure_sub_dsigma(protons_fits_div, 'ampt_new_coal_epbins1',
                                                      'ampt_new_coal_epbins1_v2cor',
                                                      v2_vals, flow_cor_coefs[div], cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    # plt.show()

    proton_fits_w_flow = protons_fits[protons_fits['name'] == 'ampt_new_coal_epbins1']
    v2_ep_e_vals, v2_ep_e_errs, v2_rp_e_vals, v2_rp_e_errs, v2_op_e_vals = [], [], [], [], []
    for energy in energies_fit:
        df_e = proton_fits_w_flow[proton_fits_w_flow['energy'] == energy]
        # print(df_e)
        divs, avgs = np.array(df_e['divs']), np.array(df_e['avg_meas'])
        # print(divs, avgs)
        coefs = np.array([flow_cor_coefs[div] for div in divs])
        v2 = v2_vals[energy][cent_plt].val
        v2_rp = v2_rp_vals[energy][cent_plt].val
        fig, ax = plt.subplots()
        v2s, chi_2s = [], []
        x_divs = np.linspace(60, 300, 1000)
        for v2_i in np.linspace(v2 * 0, v2 * 1.25, 50):
            cor_avgs = avgs - v2_divs(divs * np.pi / 180, v2_i)
            # cor_avgs = avgs - coefs * v2_i ** 2
            cor_avg_vals, cor_avg_errs = [x.val for x in cor_avgs], [x.err for x in cor_avgs]
            popt, pcov = cf(quad_180, divs, cor_avg_vals, sigma=cor_avg_errs, absolute_sigma=True)
            chi_2 = np.sum((cor_avg_vals - quad_180(divs, *popt)) ** 2 / cor_avg_errs)
            ax.scatter(divs, cor_avg_vals)
            ax.plot(x_divs, quad_180(x_divs, *popt))
            v2s.append(v2_i)
            chi_2s.append(chi_2)
        chi_2s, v2s = zip(*sorted(zip(chi_2s, v2s)))
        v2_op_e_vals.append(v2s[0])
        v2_ep_e_vals.append(v2)
        v2_ep_e_errs.append(v2_vals[energy][cent_plt].err)
        v2_rp_e_vals.append(v2_rp)
        v2_rp_e_errs.append(v2_rp_vals[energy][cent_plt].err)

        fig2, ax2 = plt.subplots()
        ax2.axhline(0, color='black')
        ax2.scatter(v2s, chi_2s)
        ax2.set_title(f'{energy}GeV')
        ax2.axvline(v2, color='green', ls='--', label='v2_ep')
        ax2.axvline(v2_rp, color='red', ls='--', label='v2_rp')
        ax2.legend()
        plt.show()
    fig3, ax3 = plt.subplots()
    ax3.grid()
    ax3.set_xlabel('Energy (GeV)')
    ax3.set_ylabel('v2')
    ax3.errorbar(energies_fit, v2_ep_e_vals, v2_ep_e_errs, ls='none', marker='o', color='green', label='v2_ep')
    ax3.errorbar(energies_fit, v2_rp_e_vals, v2_rp_e_errs, ls='none', marker='o', color='red', label='v2_rp')
    ax3.scatter(energies_fit, v2_op_e_vals, color='blue', label='v2_op')
    ax3.legend()
    # plt.show()

    # plt.show()
    # return

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    v2_divs_plt = np.setdiff1d(np.unique(df['divs']), exclude_divs + [180])
    energy_colors = {7: 'blue', 11: 'orange', 19: 'green', 27: 'red', 39: 'purple', 62: 'brown'}
    avg_v2aps = []
    for energy in energies_fit:
        fig_e, ax_e = plt.subplots(dpi=144)
        avg_v2ap = np.mean(v2_apparents[energy])
        avg_v2aps.append(avg_v2ap)
        ax_e.grid()
        ax.errorbar(v2_divs_plt, [x.val for x in v2_apparents[energy]], yerr=[x.err for x in v2_apparents[energy]],
                    marker='o', ls='none', alpha=0.6, color=energy_colors[energy], label=f'{energy}GeV')
        ax.axhspan(avg_v2ap.val + avg_v2ap.err, avg_v2ap.val - avg_v2ap.err, color=energy_colors[energy], alpha=0.3)
        ax.axhline(avg_v2ap.val, ls='--', color=energy_colors[energy])
        ax_e.errorbar(v2_divs_plt, [x.val for x in v2_apparents[energy]], yerr=[x.err for x in v2_apparents[energy]],
                      marker='o', ls='none', alpha=0.9, color=energy_colors[energy])
        ax_e.set_title(f'{energy} GeV')
        ax_e.axhspan(avg_v2ap.val + avg_v2ap.err, avg_v2ap.val - avg_v2ap.err, color=energy_colors[energy], alpha=0.5)
        ax_e.axhline(avg_v2ap.val, color=energy_colors[energy], ls='--', alpha=0.9)
        fig_e.canvas.manager.set_window_title(f'Apparent v2 {energy}GeV')
        fig_e.tight_layout()
    ax.legend()
    ax.set_ylabel('Apparent v2')
    fig.tight_layout()
    fig.canvas.manager.set_window_title('Apparent v2')

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    ax.errorbar(energies_fit, [x.val for x in avg_v2aps], [x.err for x in avg_v2aps], ls='none', marker='o',
                label='v2 Apparent')
    v2_plt_vals, v2_errs, v2_rp_plt_vals, v2_rp_errs = [], [], [], []
    for energy in energies_fit:
        v2_plt_vals.append(v2_vals[energy][cent_plt].val)
        v2_errs.append(v2_vals[energy][cent_plt].err)
        v2_rp_plt_vals.append(v2_rp_vals[energy][cent_plt].val)
        v2_rp_errs.append(v2_rp_vals[energy][cent_plt].err)
    ax.errorbar(energies_fit, v2_plt_vals, v2_errs, ls='none', marker='o', label='v2 Event Plane')
    ax.errorbar(energies_fit, v2_rp_plt_vals, v2_rp_errs, ls='none', marker='o', label='v2 Reaction Plane')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('v2')
    ax.set_title(f'V2 vs Energy Cent {cent_plt}')
    ax.legend()
    fig.tight_layout()
    fig.canvas.manager.set_window_title(f'V2 vs Energy Cent {cent_plt}')

    fig, ax = plt.subplots(dpi=144)
    ax.grid()
    v2_rp_div_app = []
    for v2_app, v2_rp_val, v2_rp_err in zip(avg_v2aps, v2_rp_plt_vals, v2_rp_errs):
        v2_rp_div_app.append(v2_app / Measure(v2_rp_val, v2_rp_err))
    ax.errorbar(energies_fit, [x.val for x in v2_rp_div_app], [x.err for x in v2_rp_div_app], ls='none', marker='o')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('v2')
    ax.set_title(f'V2 Reaction Plane / V2 Apparent vs Energy Cent {cent_plt}')
    fig.tight_layout()
    fig.canvas.manager.set_window_title(f'V2 Reaction Plane / V2 Apparent vs Energy Cent {cent_plt}')

    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    plot_dvar_avgs_divs(protons_fits, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels)
    # df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
    #                                  data_sets_labels=data_sets_labels, verbose=False)
    # df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_ampt_v2_closure_var():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_stats_ampt_flow_closure.csv'
    vs_in_dir = f'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/' \
                f'Ampt_rapid05_resample_norotate_epbins1_0/'
    df_path = base_path + df_name

    pd.set_option('max_columns', None)

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 7
    energies_fit = [7, 11, 19, 27, 39, 62]

    data_sets_plt = ['ampt_new_coal_rp', 'ampt_new_coal_epbins1', 'ampt_new_coal_epbins1_v2cor',
                     'ampt_new_coal_epbins1_v2rpcor', 'ampt_new_coal_epbins1_flow_rpcor']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red', 'purple', 'orange']))
    data_sets_labels = dict(zip(data_sets_plt, ['Reaction Plane', 'Flow Included', 'Flow Corrected ep',
                                                'Flow Corrected rp', 'Flow Corrected rp 2-6']))

    all_sets_plt = data_sets_plt

    v_orders = [2]

    # Dictionary of {order: {energy: {cent: v2}}}
    v_ep_vals = {i: read_flow_values(vs_in_dir, f'v{i}') for i in [2]}
    v_rp_vals = {i: read_flow_values(vs_in_dir, f'v{i}_rp') for i in v_orders}

    # v2_vals = read_flow_values(vs_in_dir)  # Dictionary of {energy: {cent: v2}}
    # v2_rp_vals = read_flow_values(vs_in_dir, 'v2_rp')  # Dictionary of {energy: {cent: v2}}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))  # Placeholder, need k2
    df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    # df_raw.loc[:, 'err'] = abs(2 * df_raw['val'] * df_raw['err'] / (tp * (tp - 1)))
    # df_raw.loc[:, 'val'] = (df_raw['val']**2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))
    # df_mix.loc[:, 'err'] = abs(2 * df_mix['val'] * df_mix['err'] / (tp * (tp - 1)))
    # df_mix.loc[:, 'val'] = (df_mix['val']**2 - (tp * p * (1 - p))) / (tp * (tp - 1))
    df_diff = df[(df['data_type'] == 'diff') & (df['stat'] == stat_plot)]
    p, tp = df_diff['divs'] / 360, df_diff['total_protons']
    df_diff.loc[:, 'val'] = df_diff['val'] / (tp * (tp - 1))
    df_diff.loc[:, 'err'] = df_diff['err'] / (tp * (tp - 1))
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'

    dvar_vs_protons(df_sub, div_plt, cent_plt, [39], ['sub'], all_sets_plt, plot=True, avg=True,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub, df_diff], ignore_index=True), div_plt, cent_plt, [39],
                    ['raw', 'mix', 'sub', 'diff'], all_sets_plt, plot=True, avg=True,
                    data_sets_labels=data_sets_labels)
    dvar_vs_protons_energies(df_sub, [120], cent_plt, [7, 11, 19, 27, 39, 62], ['sub'], all_sets_plt,
                             plot=True, avg=True, plot_avg=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, y_ranges=[-0.0014, 0.0011])

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets_plt, plot=False,
                                           avg=True)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_v2rpcor', {2: v_rp_vals[2]},
                                                div, cent_plt)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_v2cor',
                                                v_ep_vals, div, cent_plt)
        protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
                                                'ampt_new_coal_epbins1_flow_rpcor', v_rp_vals,
                                                div, cent_plt)
        # protons_fits_div = subtract_dsigma_flow(protons_fits_div, 'ampt_new_coal_epbins1',
        #                                                 'ampt_new_coal_epbins1_flow_cor',
        #                                                 v_ep_vals, div, cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    plot_dvar_avgs_divs(protons_fits, data_sets_plt, data_sets_colors=data_sets_colors, fit=False,
                        data_sets_labels=data_sets_labels, alpha=0.6)

    plt.show()


def plot_flow_v2_closure():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_flow.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = None  # 'cf_tprotons_fits.csv'
    df_partitions_fits_name = None  # 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 6
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['divide']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res15_v207', 'flow_resample_res15_v207_v2sub']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red']))
    data_sets_labels = dict(zip(data_sets_plt, ['Flow Excluded', 'Flow Included', 'Flow Subtracted']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    stat_vs_protons(df, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
                             plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
                             data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, data_types_plt, all_sets_plt,
                                           plot=False, fit=True)
        print(protons_fits_div)
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div, 'flow_resample_res15_v207',
                                               'flow_resample_res15_v207_v2sub', 0.07, flow_cor_coefs[div])
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_flow_v2_closure_raw():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = '/media/ucla/Research/Results/Azimuth_Analysis/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'Binomial_Slice_Moments/binom_slice_stats_flow.csv'
    flow_cor_dir = 'F:/Research/Results/Flow_Correction/'
    fits_out_base = 'Base_Zero_Fits'
    df_tproton_fits_name = None  # 'cf_tprotons_fits.csv'
    df_partitions_fits_name = None  # 'cf_partitions_fits.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'standard deviation'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 6
    energies_fit = [7, 11, 19, 27, 39, 62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    pd.set_option('max_columns', None)

    data_sets_plt = ['flow_resample_res99_v207', 'flow_resample_res15_v207', 'flow_resample_res15_v207_v2sub']
    data_sets_colors = dict(zip(data_sets_plt, ['green', 'blue', 'red']))
    data_sets_labels = dict(zip(data_sets_plt, ['Flow Excluded', 'Flow Included', 'Flow Subtracted']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    flow_cor_coefs = read_v2_slope_coefs(flow_cor_dir)  # Dictionary of {div: coef}

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[df['data_type'] == 'raw']
    p = df_raw['divs'] / 360
    df_raw.loc[:, 'val'] = df_raw['val'] / np.sqrt((df_raw['total_protons'] * p * (1 - p)))
    df_mix = df[df['data_type'] == 'mix']
    p = df_mix['divs'] / 360
    df_mix.loc[:, 'val'] = df_mix['val'] / np.sqrt((df_mix['total_protons'] * p * (1 - p)))
    print(df)
    print(df_raw)

    stat_vs_protons(df_raw, stat_plot, div_plt, cent_plt, [62], ['raw'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    stat_vs_protons(df_mix, stat_plot, div_plt, cent_plt, [62], ['mix'], all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels, star_prelim=True,
                    y_ranges={'standard deviation': [0.946, 1.045]})
    # stat_vs_protons_energies(df, stat_plot, [120], cent_plt, [7, 11, 19, 27, 39, 62], data_types_plt, all_sets_plt,
    #                          plot=True, fit=True, plot_fit=True, data_sets_colors=data_sets_colors,
    #                          data_sets_labels=data_sets_labels, star_prelim=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = stat_vs_protons(df_raw, stat_plot, div, cent_plt, energies_fit, ['raw'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_mix = stat_vs_protons(df_mix, stat_plot, div, cent_plt, energies_fit, ['mix'],
                                               all_sets_plt, plot=False, fit=True)
        # print(protons_fits_div_raw)
        for index, row in protons_fits_div_mix.iterrows():
            if data_sets_labels[row['name']] == 'Flow Included':
                continue
            loc_index = (protons_fits_div_raw.name == row['name']) & (protons_fits_div_raw.energy == row['energy'])
            raw_sub_meas = (protons_fits_div_raw.loc[loc_index, 'slope_meas'] - row['slope_meas']).iloc[0]
            # print(row, raw_sub_meas)
            protons_fits_div_raw.loc[loc_index, 'slope_meas'] = raw_sub_meas
            protons_fits_div_raw.loc[loc_index, 'slope'] = raw_sub_meas.val
            protons_fits_div_raw.loc[loc_index, 'slope_err'] = raw_sub_meas.err

        # print(protons_fits_div_raw)
        # input()
        # print(protons_fits_div_mix)
        protons_fits_div = ampt_v2_closure_sub(protons_fits_div_raw, 'flow_resample_res15_v207',
                                               'flow_resample_res15_v207_v2sub', 0.07, flow_cor_coefs[div], cent_plt)
        protons_fits.append(protons_fits_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    if df_tproton_fits_name:
        protons_fits.to_csv(f'{base_path}{fits_out_base}{df_tproton_fits_name}', index=False)
    print(protons_fits)
    print(pd.unique(protons_fits['amp']))
    print(pd.unique(protons_fits['spread']))
    df_fits = plot_protons_fits_divs(protons_fits, all_sets_plt, data_sets_colors=data_sets_colors, fit=True,
                                     data_sets_labels=data_sets_labels, verbose=False)
    df_fits.to_csv(f'{base_path}{fits_out_base}{df_partitions_fits_name}', index=False)
    # print(df_fits)
    # plot_slope_div_fits(df_fits, data_sets_colors, data_sets_labels)
    # plot_slope_div_fits_simpars(df_fits)

    plt.show()


def plot_flow_eff_test():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_var_cent8_2source_closure_tests.csv'
    df_path = base_path + df_name
    sim_sets = []

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    data_sets_plt = ['flow_eff_res15_v207', 'flow_res15_v207']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'blue']))
    data_sets_labels = dict(zip(data_sets_plt, ['flow_eff_test', 'flow_test']))

    all_sets_plt = data_sets_plt + sim_sets[:]

    df = pd.read_csv(df_path)
    df = df.dropna()
    print(pd.unique(df['name']))

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    df_raw = df[df['data_type'] == 'raw']
    p = df_raw['divs'] / 360
    df_raw.loc[:, 'val'] = df_raw['val'] / (df_raw['total_protons'] * p * (1 - p))
    df_raw.loc[:, 'err'] = df_raw['err'] / (df_raw['total_protons'] * p * (1 - p))
    df_mix = df[df['data_type'] == 'mix']
    p = df_mix['divs'] / 360
    df_mix.loc[:, 'val'] = df_mix['val'] / (df_mix['total_protons'] * p * (1 - p))
    df_mix.loc[:, 'err'] = df_mix['err'] / (df_mix['total_protons'] * p * (1 - p))
    print(df)
    print(df_raw)

    stat_vs_protons(df_raw, stat_plot, div_plt, cent_plt, [62], data_types_plt, all_sets_plt, plot=True, fit=False,
                    data_sets_colors=data_sets_colors, data_sets_labels=data_sets_labels,
                    y_ranges={'k2': [0.88, 1.09]})

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = stat_vs_protons(df_raw, stat_plot, div, cent_plt, energies_fit, ['raw'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        protons_fits_div_mix = stat_vs_protons(df_mix, stat_plot, div, cent_plt, energies_fit, ['mix'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        protons_fits.append(protons_fits_div_mix)

        sub_meas = protons_fits_div_raw['slope_meas'] - protons_fits_div_mix['slope_meas']
        protons_fits_div_sub = protons_fits_div_raw.copy()
        protons_fits_div_sub.loc[:, 'slope'] = [x.val for x in sub_meas]
        protons_fits_div_sub.loc[:, 'slope_err'] = [x.err for x in sub_meas]
        protons_fits_div_sub.loc[:, 'slope_meas'] = sub_meas
        protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_mix['name'] + '_sub'
        protons_fits.append(protons_fits_div_sub)

        protons_fits_div_div = stat_vs_protons(df, stat_plot, div, cent_plt, energies_fit, ['divide'],
                                               all_sets_plt, plot=False, fit=True)
        protons_fits_div_div.loc[:, 'name'] = protons_fits_div_div['name'] + '_div'
        protons_fits.append(protons_fits_div_div)
    protons_fits = pd.concat(protons_fits, ignore_index=True)
    for data_set in data_sets_plt:
        data_sets = [data_set + x for x in ['_raw', '_mix', '_div', '_sub']]
        colors = dict(zip(data_sets, ['blue', 'green', 'red', 'purple']))
        labels = dict(zip(data_sets, [data_sets_labels[data_set] + x for x in [' Raw', ' Mix', ' Div', ' Sub']]))
        plot_protons_fits_divs(protons_fits, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels)
    data_sets = ['flow_eff_res15_v207_raw', 'flow_eff_res15_v207_div', 'flow_res15_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Flow+Efficiency Raw', 'Flow+Efficiency Raw/Mix', 'Flow Raw']))
    plot_protons_fits_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels)

    plt.show()


def plot_anticl_flow_closure_test():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_flow_anticl_convo_test.csv'
    df_name = 'binom_slice_var_cent8_2source_closure_tests.csv'
    df_name = 'binom_slice_var_cent8_v2_anticl_closure.csv'
    # save_fits = False
    # v2_fit_out_dir = 'F:/Research/Results/Flow_Correction/'
    # v2_fit_out_dir = None
    # fits_out_base = 'Base_Zero_Fits/'
    # df_tproton_fits_name = None  # 'flow_tprotons_fits.csv'
    # df_partitions_fits_name = 'flow_partitions_fits.csv'
    df_path = base_path + df_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    # data_sets_plt = ['flow_resample_res15_v207', 'anticlmulti_resample_s05_a05',
    #                  'anticlflow_resample_res15_v207_s05_a05']
    # data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue']))
    # data_sets_labels = dict(zip(data_sets_plt, ['flow', 'anticlustering', 'both']))

    df = pd.read_csv(df_path)
    df = df.dropna()

    # Delete this block if I append v2 value to end of name in calc script
    # df_no_flow = df[df['name'].str.contains('anticlmulti')]
    # df_flow = df[df['name'].str.contains('anticlflow')]
    # df_flow.loc[:, 'name'] = df_flow['name'] + '_v207'
    # df = pd.concat([df_flow, df_no_flow], ignore_index=True)

    all_sets = pd.unique(df['name'])
    print(all_sets)

    # df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    # p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    # df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    # df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    # df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    # p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    # df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    # df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))

    df_raw, df_mix = calc_dsigma(df[df['stat'] == stat_plot], data_types=['raw', 'mix'])

    dvar_vs_protons(df_raw, div_plt, cent_plt, [62], ['raw'], all_sets, plot=True, avg=True)
    dvar_vs_protons(df_mix, div_plt, cent_plt, [62], ['mix'], all_sets, plot=True, avg=True)
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    dvar_vs_protons(df_sub, div_plt, cent_plt, [62], ['sub'], all_sets, plot=True, avg=True)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub], ignore_index=True), div_plt, cent_plt, [62],
                    ['raw', 'mix', 'sub'], all_sets, plot=True, avg=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets, plot=False,
                                               avg=True)
        # protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        # protons_fits_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets, plot=False,
        #                                        avg=True)
        # protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        # protons_fits.append(protons_fits_div_mix)
        #
        # protons_fits_div_sub = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets, plot=False,
        #                                        avg=True)
        # protons_fits.append(protons_fits_div_sub.copy())
        # protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_sub['name'] + '_sub'
        # protons_fits.append(protons_fits_div_sub)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    # plot_protons_fits_divs_flow(protons_fits, all_sets)

    for flow_set in [set_name for set_name in all_sets if 'anticlflow' in set_name and '_eff_' not in set_name]:
        no_flow_set = flow_set[:flow_set.find('_v2')].replace('anticlflow', 'anticlmulti')

        print(pd.unique(protons_fits['name']))
        print(flow_set)
        anticl_plus_v2 = protons_fits[protons_fits['name'] == flow_set]
        fits = [anticl_plus_v2.copy(), protons_fits[protons_fits['name'] == no_flow_set]]
        divs, avgs = np.array(anticl_plus_v2['divs']), np.array(anticl_plus_v2['avg_meas'])
        print(divs, avgs)
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        v2 = 0.07
        v2s, chi_2s = [], []
        x_divs = np.linspace(60, 300, 1000)
        for v2_i in np.linspace(0, v2 * 1.25, 50):
            cor_avgs = avgs - v2_divs(divs * np.pi / 180, v2_i)
            # print(v2_divs(divs * np.pi / 180, v2_i))
            # print(cor_slopes)
            cor_avg_vals, cor_avg_errs = [x.val for x in cor_avgs], [x.err for x in cor_avgs]
            ax2.errorbar(divs, cor_avg_vals, yerr=cor_avg_errs, ls='none', alpha=0.4, marker='o')
            popt, pcov = cf(quad_180, divs, cor_avg_vals, sigma=cor_avg_errs, absolute_sigma=True)
            chi_2 = np.sum((cor_avg_vals - quad_180(divs, *popt)) ** 2 / cor_avg_errs)
            ax.scatter(divs, cor_avg_vals)
            ax.plot(x_divs, quad_180(x_divs, *popt))
            v2s.append(v2_i)
            chi_2s.append(chi_2)
        fig2, ax2 = plt.subplots()
        ax2.axhline(0, color='black')
        ax2.scatter(v2s, chi_2s)
        chi_2s, v2s = zip(*sorted(zip(chi_2s, v2s)))
        # anticl_plus_v2['name'] = f'corrected_v2{str(v2s[0])[2:]}_'
        anticl_plus_v2.loc[:, 'name'] = 'corrected_v2'
        # new_slope = anticl_plus_v2['slope'] - v2_divs(anticl_plus_v2['divs'] * np.pi / 180, v2s[0])
        # anticl_plus_v2 = anticl_plus_v2.assign(slope=new_slope)
        anticl_plus_v2.loc[:, 'avg'] = anticl_plus_v2['avg'] - v2_divs(anticl_plus_v2['divs'] * np.pi / 180, v2s[0])
        print(anticl_plus_v2)
        print(v2s[0], chi_2s[0])
        fits = pd.concat([*fits, anticl_plus_v2], ignore_index=True)
        # data_sets_colors.update({'corrected_v2': 'olive'})
        plot_dsigma_fits_divs_flow(fits, [flow_set, no_flow_set] + ['corrected_v2'])

    plt.show()


def plot_anticl_flow_closure_test_simple():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_vars_2source_tests.csv'
    df_path = base_path + df_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    df = pd.read_csv(df_path)
    df = df.dropna()

    # Delete this block if I append v2 value to end of name in calc script
    # df_no_flow = df[df['name'].str.contains('anticlmulti')]
    # df_flow = df[df['name'].str.contains('anticlflow')]
    # df_flow.loc[:, 'name'] = df_flow['name'] + '_v207'
    # df = pd.concat([df_flow, df_no_flow], ignore_index=True)

    all_sets = pd.unique(df['name'])
    print(all_sets)

    df_raw, df_mix, df_diff = calc_dsigma(df[df['stat'] == stat_plot], data_types=['raw', 'mix', 'diff'])

    dvar_vs_protons(df_raw, div_plt, cent_plt, [62], ['raw'], all_sets, plot=True, avg=True)
    dvar_vs_protons(df_mix, div_plt, cent_plt, [62], ['mix'], all_sets, plot=True, avg=True)

    dvar_vs_protons(df_diff, div_plt, cent_plt, [62], ['diff'], all_sets, plot=True, avg=True)
    dvar_vs_protons(pd.concat([df_raw, df_mix, df_diff], ignore_index=True), div_plt, cent_plt, [62],
                    ['raw', 'mix', 'diff'], all_sets, plot=True, avg=True)

    dsigma_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsigma_avgs_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets, plot=False, avg=True)
        dsigma_avgs.append(dsigma_avgs_raw)
    dsigma_avgs = pd.concat(dsigma_avgs, ignore_index=True)

    for flow_set in [set_name for set_name in all_sets if 'anticlflow' in set_name and '_eff_' not in set_name]:
        no_flow_set = flow_set[:flow_set.find('_v2')].replace('anticlflow', 'anticlmulti')

        print(flow_set)
        anticl_plus_v2 = dsigma_avgs[dsigma_avgs['name'] == flow_set]
        fits = [anticl_plus_v2.copy(), dsigma_avgs[dsigma_avgs['name'] == no_flow_set]]
        fits = pd.concat([anticl_plus_v2, *fits], ignore_index=True)
        # data_sets_colors.update({'corrected_v2': 'olive'})
        data_sets_labels = None
        if '_s01_a01' in no_flow_set:
            data_sets_labels = {no_flow_set: 'Pure Repulsive Model', flow_set: 'Repulsive Model + Flow'}
        plot_dsigma_fits_divs_flow(fits, [no_flow_set, flow_set], data_sets_labels=data_sets_labels)

    # mu = np.pi
    # sigma = 0.5
    # amp = 0.2
    # func = base_gaus_pdf_wrap
    # func_args = (mu, sigma, amp)
    # func_name = 'Gaussian Cluster'
    # widths = np.linspace(0, 2 * np.pi, 250)
    # dsigmas = []
    # for width in widths:
    #     dsigmas.append(get_partition_variance(func, func_args, width)[0])
    # fig, ax = plt.subplots()
    # ax.plot(np.rad2deg(widths), dsigmas)
    # df_simple_clust = dsigma_avgs[dsigma_avgs['name'] == 'simpleclust']
    # ax.errorbar(df_simple_clust['divs'], df_simple_clust['avg'], yerr=df_simple_clust['avg_err'], ls='none', marker='o',
    #             alpha=0.6)

    simplecl_flow_sets = [x for x in all_sets if 'simpleclust_flow' in x]
    plot_dvar_avgs_divs(dsigma_avgs, simplecl_flow_sets, fit=False, plot_energy_panels=False, plot_indiv=True)
    for flow_set in simplecl_flow_sets:
        df_set = dsigma_avgs[dsigma_avgs['name'] == flow_set]
        divs = np.array(df_set['divs'])
        avgs = np.array(df_set['avg'])
        v2 = 0
        for element in flow_set.split('_'):
            if 'v2' in element:
                v2 = float('0.' + element.strip('v2'))
                break
        y_sub = avgs - v2_divs(divs / 180 * np.pi, v2)
        plt.errorbar(divs, y_sub, yerr=df_set['avg_err'], ls='none',
                     marker='x', label=f'{flow_set} - v2(div)')

    # for flow_set in [set_name for set_name in all_sets if 'simpleclust_flow' in set_name and '_eff_' not in set_name]:
    #     no_flow_set = flow_set[:flow_set.find('_v2')].replace('anticlflow', 'anticlmulti')
    #
    #     print(pd.unique(dsigma_avgs['name']))
    #     print(flow_set)
    #     anticl_plus_v2 = dsigma_avgs[dsigma_avgs['name'] == flow_set]
    #     fits = [anticl_plus_v2.copy(), dsigma_avgs[dsigma_avgs['name'] == no_flow_set]]
    #     fits = pd.concat([anticl_plus_v2, *fits], ignore_index=True)
    #     # data_sets_colors.update({'corrected_v2': 'olive'})
    #     data_sets_labels = None
    #     if '_s01_a01' in no_flow_set:
    #         data_sets_labels = {no_flow_set: 'Pure Repulsive Model', flow_set: 'Repulsive Model + Flow'}
    #     plot_dsigma_fits_divs_flow(fits, [no_flow_set, flow_set], data_sets_labels=data_sets_labels)

    plt.show()


def plot_efficiency_closure_tests():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'C:/Users/Dylan/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    # base_path = 'D:/Transfer/Research/Results/Azimuth_Analysis/'
    df_name = 'binom_slice_var_cent8_efficiency_closure.csv'
    df2_name = 'binom_slice_var_cent8_v2_anticlindep_closure.csv'
    # df_name = 'binom_slice_var_cent8_simpleclust_test.csv'
    df_path = base_path + df_name
    df2_path = base_path + df2_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    samples = 72  # For title purposes only

    df = pd.read_csv(df_path)
    df = pd.concat([df, pd.read_csv(df2_path)], ignore_index=True)
    df = df.dropna()

    # df = df[df['name'].str.contains('simpleclust_eff')]

    all_sets = pd.unique(df['name'])
    print(all_sets)

    df['energy'] = df.apply(lambda row: 'sim' if 'sim_' in row['name'] else row['energy'], axis=1)

    # df_raw = df[(df['data_type'] == 'raw') & (df['stat'] == stat_plot)]
    # p, tp = df_raw['divs'] / 360, df_raw['total_protons']
    # df_raw.loc[:, 'val'] = (df_raw['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    # df_raw.loc[:, 'err'] = df_raw['err'] / (tp * (tp - 1))
    # df_mix = df[(df['data_type'] == 'mix') & (df['stat'] == stat_plot)]
    # p, tp = df_mix['divs'] / 360, df_mix['total_protons']
    # df_mix.loc[:, 'val'] = (df_mix['val'] - (tp * p * (1 - p))) / (tp * (tp - 1))
    # df_mix.loc[:, 'err'] = df_mix['err'] / (tp * (tp - 1))

    df_raw, df_mix = calc_dsigma(df[df['stat'] == stat_plot], data_types=['raw', 'mix'])

    # dvar_vs_protons(df_raw, div_plt, cent_plt, [62], ['raw'], all_sets, plot=True, avg=True)
    # dvar_vs_protons(df_mix, div_plt, cent_plt, [62], ['mix'], all_sets, plot=True, avg=True)
    df_sub = subtract_avgs(df_raw.drop(columns=['data_type']), df_mix.drop(columns=['data_type']),
                           val_col='val', err_col='err')
    df_sub['data_type'] = 'sub'
    # dvar_vs_protons(df_sub, div_plt, cent_plt, [62], ['sub'], all_sets, plot=True, avg=True)
    # dvar_vs_protons(pd.concat([df_raw, df_mix, df_sub], ignore_index=True), div_plt, cent_plt, [62],
    #                 ['raw', 'mix', 'sub'], all_sets, plot=True, avg=True)

    protons_fits = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        protons_fits_div_raw = dvar_vs_protons(df_raw, div, cent_plt, energies_fit, ['raw'], all_sets, plot=False,
                                               avg=True)
        protons_fits_div_raw.loc[:, 'name'] = protons_fits_div_raw['name'] + '_raw'
        protons_fits.append(protons_fits_div_raw)

        protons_fits_div_mix = dvar_vs_protons(df_mix, div, cent_plt, energies_fit, ['mix'], all_sets, plot=False,
                                               avg=True)
        protons_fits_div_mix.loc[:, 'name'] = protons_fits_div_mix['name'] + '_mix'
        protons_fits.append(protons_fits_div_mix)

        protons_fits_div_sub = dvar_vs_protons(df_sub, div, cent_plt, energies_fit, ['sub'], all_sets, plot=False,
                                               avg=True)
        protons_fits.append(protons_fits_div_sub.copy())
        protons_fits_div_sub.loc[:, 'name'] = protons_fits_div_sub['name'] + '_sub'
        protons_fits.append(protons_fits_div_sub)
    protons_fits = pd.concat(protons_fits, ignore_index=True)

    # for data_set in all_sets:
    #     data_sets = [data_set + x for x in ['_raw', '_mix', '_sub']]
    #     colors = dict(zip(data_sets, ['blue', 'green', 'red']))
    #     labels = dict(zip(data_sets, [data_set + x for x in [' Raw', ' Mix', ' Sub']]))
    #     plot_dvar_avgs_divs(protons_fits, data_sets, data_sets_colors=colors, fit=False, data_sets_labels=labels,
    #                         plot_energy_panels=False, title=data_set)

    plt.rcParams["figure.figsize"] = (8, 4)

    data_sets = ['flow_eff_res15_v207_raw', 'flow_eff_res15_v207_sub', 'flow_res15_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Flow + Efficiency Single', 'Flow + Efficiency Single-Mixed', 'Flow Single']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plot_energy_panels=False, alpha=0.6, title='Flow v2=0.07 Efficiency Correction',
                        ylab=r'$\langle\Delta\sigma^2\rangle$')
    xs = np.linspace(0, 360, 1000)
    plt.plot(xs, vn_divs(np.deg2rad(xs), 0.07, 2), color='gray', label='Flow Analytical')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.7))

    data_sets = ['simpleclust_eff_raw', 'simpleclust_eff_sub', 'simpleclust_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Simple Clustering + Efficiency Single',
                                  'Simple Clustering + Efficiency Single-Mixed', 'Simple Clustering Single']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plot_energy_panels=False, alpha=0.6, title='Simple Clustering (s05, a2) Efficiency Correction',
                        ylab=r'$\langle\Delta\sigma^2\rangle$')

    for s, a in [('1', '01'), ('08', '01'), ('05', '01'), ('01', '01')]:
        data_sets = [f'anticlflow_eff_s{s}_a{a}_raw', f'anticlflow_eff_s{s}_a{a}_sub', f'anticlmulti_s{s}_a{a}_raw']
        colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
        labels = dict(zip(data_sets, ['Repulsive Model + Efficiency Single',
                                      'Repulsive Model + Efficiency Single-Mixed', 'Repulsive Model Single']))
        s, a = round(float('0.' + s) * 10, 1), round(float('0.' + a) * 10, 1)
        plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                            plot_energy_panels=False, alpha=0.6, ylab=r'$\langle\Delta\sigma^2\rangle$',
                            title=rf'Repulsive Simulation ($\sigma={s}$, $A={a}$) Efficiency Correction')

    for s, a in [('1', '01'), ('08', '01'), ('05', '01'), ('01', '01')]:
        data_sets = [f'anticlflowindep_eff_s{s}_a{a}_raw', f'anticlflowindep_eff_s{s}_a{a}_sub',
                     f'anticlmulti_s{s}_a{a}_raw']
        colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
        labels = dict(zip(data_sets, ['Repulsive Model (Independent) + Efficiency Single',
                                      'Repulsive Model (Independent) + Efficiency Single-Mixed',
                                      'Repulsive Model Single']))
        s, a = round(float('0.' + s) * 10, 1), round(float('0.' + a) * 10, 1)
        plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                            plot_energy_panels=False, alpha=0.6, ylab=r'$\langle\Delta\sigma^2\rangle$',
                            title=rf'Repulsive Model (Independent) ($\sigma={s}$, $A={a}$) Efficiency Correction')

    data_sets = [f'anticlflowindep_eff_s05_a01_v207_raw', f'anticlflowindep_eff_s05_a01_v207_sub',
                 f'anticlflowindep_s05_a01_v207_raw']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Repulsive Model (Independent) + Efficiency + Flow Single',
                                  'Repulsive Model (Independent) + Efficiency + Flow Single-Mixed',
                                  'Repulsive Model + Flow Single']))
    s, a = round(float('0.' + '05') * 10, 1), round(float('0.' + '01') * 10, 1)
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plot_energy_panels=False, alpha=0.6, ylab=r'$\langle\Delta\sigma^2\rangle$',
                        title=rf'Repulsive Model (Independent) ($\sigma={s}$, $A={a}$) + \
                        Flow v2=0.07 Efficiency Correction')

    data_sets = ['anticlflow_eff_s1_a01_v207_raw', 'anticlflow_eff_s1_a01_v207_sub', '']
    colors = dict(zip(data_sets, ['blue', 'red', 'orange']))
    labels = dict(zip(data_sets, ['Repulsive Model + Flow + Efficiency Single',
                                  'Repulsive Model + Flow + Efficiency Single-Mixed',
                                  'Repulsive Model + Flow Single']))
    plot_dvar_avgs_divs(protons_fits, data_sets, fit=False, data_sets_colors=colors, data_sets_labels=labels,
                        plot_energy_panels=False, alpha=0.6, title='Flow v2=0.07 Efficiency Correction',
                        ylab=r'$\langle\Delta\sigma^2\rangle$')

    plt.show()


def plot_closure_tests():
    plt.rcParams["figure.figsize"] = (8, 5)
    plt.rcParams["figure.dpi"] = 144
    base_path = 'F:/Research/Results/Azimuth_Analysis/Binomial_Slice_Moments/'
    df_name = 'binom_slice_vars_2source_tests.csv'
    df_path = base_path + df_name

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [62]
    data_types_plt = ['raw']
    samples = 72  # For title purposes only

    df = pd.read_csv(df_path)

    all_sets = pd.unique(df['name'])
    print(all_sets)

    df_raw, df_mix, df_diff = calc_dsigma(df[df['stat'] == stat_plot], data_types=['raw', 'mix', 'diff'])
    df_types = pd.concat([df_raw, df_mix, df_diff], ignore_index=True)

    dsigma_avgs = []
    for div in np.setdiff1d(np.unique(df['divs']), exclude_divs):  # All divs except excluded
        print(f'Div {div}')
        dsigma_avgs_types = dvar_vs_protons(df_types, div, cent_plt, energies_fit, ['raw', 'mix', 'diff'], all_sets,
                                            plot=False, avg=True)
        dsigma_avgs.append(dsigma_avgs_types)
    dsigma_avgs = pd.concat(dsigma_avgs, ignore_index=True)

    simple_clust_amps = ['3', '05']
    amp_map = {'05': '0.165', '3': '1.34'}
    for amp in simple_clust_amps:
        amp_val = amp_map[amp]
        simple_clust = dsigma_avgs[(dsigma_avgs['name'] == f'simpleclust_s08_a{amp}')
                                   & (dsigma_avgs['data_type'] == 'raw')]
        simple_clust_eff = dsigma_avgs[(dsigma_avgs['name'] == f'simpleclust_s08_a{amp}_eff')
                                       & (dsigma_avgs['data_type'] == 'raw')]
        simple_clust_eff_mix = dsigma_avgs[(dsigma_avgs['name'] == f'simpleclust_s08_a{amp}_eff')
                                           & (dsigma_avgs['data_type'] == 'mix')]
        simple_clust_eff_cor = dsigma_avgs[(dsigma_avgs['name'] == f'simpleclust_s08_a{amp}_eff') &
                                           (dsigma_avgs['data_type'] == 'diff')]
        plot_closures(simple_clust, simple_clust_eff, simple_clust_eff_cor, alpha=0.7, df_bkg=simple_clust_eff_mix,
                      title=rf'Simple Clustering (A={amp_val}, $\sigma$=0.8) Efficiency Correction')

        v2s = ['1', '07', '05', '03']
        for v2 in v2s:
            v2_val = float(f'0.{v2}')
            flow_clust = dsigma_avgs[(dsigma_avgs['name'] == f'simpleclust_s08_a{amp}_flow_v2{v2}') &
                                     (dsigma_avgs['data_type'] == 'raw')]
            flow_clust_cor, flow_clust_v2 = flow_clust.copy(), flow_clust.copy()
            flow_clust_v2['avg_meas'] = flow_clust['avg_meas'] * 0 + v2_divs(np.deg2rad(flow_clust['divs']), v2_val)
            flow_clust_v2['avg'] = flow_clust['avg'] * 0 + v2_divs(np.deg2rad(flow_clust['divs']), v2_val)
            flow_clust_v2['avg_err'] = flow_clust['avg_err'] * 0
            flow_clust_cor['avg_meas'] = flow_clust['avg_meas'] - v2_divs(np.deg2rad(flow_clust['divs']), v2_val)
            flow_clust_cor['avg'] = flow_clust_cor['avg_meas'].apply(lambda x: x.val)
            flow_clust_cor['avg_err'] = flow_clust_cor['avg_meas'].apply(lambda x: x.err)
            plot_closures(simple_clust, flow_clust, flow_clust_cor, alpha=0.7, df_bkg=flow_clust_v2,
                          title=rf'Simple Clustering (A={amp_val}, $\sigma$=0.8) Flow (v2={v2_val:.2f}) Correction')

        # Flow and efficiency with simple clustering
        v2_val = float(f'0.07')
        simple_clust = dsigma_avgs[(dsigma_avgs['name'] == 'simpleclust_s08_a05') & (dsigma_avgs['data_type'] == 'raw')]
        flow_eff_clust_raw = dsigma_avgs[(dsigma_avgs['name'] == 'simpleclust_s08_a05_flow_v207_eff')
                                         & (dsigma_avgs['data_type'] == 'raw')]
        simple_clust_flow_eff_mix = dsigma_avgs[(dsigma_avgs['name'] == 'simpleclust_s08_a05_eff')
                                                & (dsigma_avgs['data_type'] == 'mix')]
        simple_clust_flow_eff_cor = dsigma_avgs[(dsigma_avgs['name'] == 'simpleclust_s08_a05_flow_v207_eff') &
                                                (dsigma_avgs['data_type'] == 'diff')]
        flow_clust_cor, flow_clust_v2 = simple_clust_flow_eff_cor.copy(), simple_clust_flow_eff_cor.copy()
        flow_clust_v2['avg_meas'] = (simple_clust_flow_eff_cor['avg_meas'] * 0 +
                                     v2_divs(np.deg2rad(simple_clust_flow_eff_cor['divs']), v2_val))
        flow_clust_v2['avg'] = (simple_clust_flow_eff_cor['avg'] * 0 +
                                v2_divs(np.deg2rad(simple_clust_flow_eff_cor['divs']), v2_val))
        flow_clust_v2['avg_err'] = simple_clust_flow_eff_cor['avg_err'] * 0
        flow_clust_cor['avg_meas'] = (simple_clust_flow_eff_cor['avg_meas'] -
                                      v2_divs(np.deg2rad(simple_clust_flow_eff_cor['divs']), v2_val))
        flow_clust_cor['avg'] = flow_clust_cor['avg_meas'].apply(lambda x: x.val)
        flow_clust_cor['avg_err'] = flow_clust_cor['avg_meas'].apply(lambda x: x.err)
        title = rf'Simple Clustering (A={amp_val}, $\sigma$=0.8) Flow (v2={v2_val:.2f}) + Efficiency Correction'
        plot_closures(simple_clust, flow_eff_clust_raw, flow_clust_cor, title=title,
                      df_bkg=[simple_clust_flow_eff_mix, flow_clust_v2])

    # Anti-clustering
    amp = '01'
    amp_map = {'01': '-0.01'}
    anti_clust_sigmas = ['01', '05', '08', '1']
    sigma_map = {'01': '0.1', '05': '0.5', '08': '0.8', '1': '1'}
    v2 = '07'
    v2_val = float(f'0.{v2}')
    for sigma in anti_clust_sigmas:
        amp_val = amp_map[amp]
        sigma_val = sigma_map[sigma]
        anti_clust = dsigma_avgs[(dsigma_avgs['name'] == f'anticlmulti_s{sigma}_a{amp}')
                                 & (dsigma_avgs['data_type'] == 'raw')]
        anti_clust_eff = dsigma_avgs[(dsigma_avgs['name'] == f'anticlflow_eff_s{sigma}_a{amp}')
                                     & (dsigma_avgs['data_type'] == 'raw')]
        anti_clust_eff_mix = dsigma_avgs[(dsigma_avgs['name'] == f'anticlflow_eff_s{sigma}_a{amp}')
                                         & (dsigma_avgs['data_type'] == 'mix')]
        anti_clust_eff_cor = dsigma_avgs[(dsigma_avgs['name'] == f'anticlflow_eff_s{sigma}_a{amp}') &
                                         (dsigma_avgs['data_type'] == 'diff')]
        plot_closures(anti_clust, anti_clust_eff, anti_clust_eff_cor, alpha=0.7, df_bkg=anti_clust_eff_mix,
                      title=rf'Anti-Clustering (A={amp_val}, $\sigma$={sigma_val}) Efficiency Correction')

        flow_anti_clust = dsigma_avgs[(dsigma_avgs['name'] == f'anticlflow_s{sigma}_a{amp}_v2{v2}') &
                                      (dsigma_avgs['data_type'] == 'raw')]
        flow_anti_clust_cor, flow_anti_clust_v2 = flow_anti_clust.copy(), flow_anti_clust.copy()
        flow_anti_clust_v2['avg_meas'] = (flow_anti_clust['avg_meas'] * 0 +
                                          v2_divs(np.deg2rad(flow_anti_clust['divs']), v2_val))
        flow_anti_clust_v2['avg'] = flow_anti_clust['avg'] * 0 + v2_divs(np.deg2rad(flow_anti_clust['divs']), v2_val)
        flow_anti_clust_v2['avg_err'] = flow_anti_clust['avg_err'] * 0
        flow_anti_clust_cor['avg_meas'] = (flow_anti_clust['avg_meas'] -
                                           v2_divs(np.deg2rad(flow_anti_clust['divs']), v2_val))
        flow_anti_clust_cor['avg'] = flow_anti_clust_cor['avg_meas'].apply(lambda x: x.val)
        flow_anti_clust_cor['avg_err'] = flow_anti_clust_cor['avg_meas'].apply(lambda x: x.err)
        plot_closures(anti_clust, flow_anti_clust, flow_anti_clust_cor, alpha=0.7, df_bkg=flow_anti_clust_v2,
                      title=rf'Anti-Clustering (A={amp_val}, $\sigma$={sigma_val}) Flow (v2={v2_val:.2f}) Correction')

    plt.show()


def plot_raw_mix_sys_comp():
    plt.rcParams["figure.figsize"] = (6.66, 5)
    plt.rcParams["figure.dpi"] = 144

    presentation_mode = False
    if presentation_mode:
        # plt.rcParams['axes.labelsize'] = 14  # Adjust the value as needed
        # plt.rcParams['axes.titlesize'] = 16  # Adjust the value as needed
        # plt.rcParams['legend.fontsize'] = 14  # Adjust the value as needed
        plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.2
    # else:
    #     plt.rcParams['font.size'] = plt.rcParams['font.size'] * 1.1

    base_path = 'F:/Research/Results/Azimuth_Analysis/'
    # base_path = 'C:/Users/Dyn04/OneDrive - personalmicrosoftsoftware.ucla.edu/OneDrive - UCLA IT Services/Research/UCLA/Results/Azimuth_Analysis/'
    # base_path = 'C:/Users/Dyn04/Research/'
    df_name = 'Binomial_Slice_Moments/binom_slice_vars_bes_sys.csv'
    df_model_name = 'Bes_with_Sys/binom_slice_vars_model.csv'
    df_dsigma_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma.csv'
    df_dsigma_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma.csv'
    df_dsigma_v2sub_name = 'Bes_with_Sys/binom_slice_vars_bes_dsigma_v2sub.csv'
    df_dsigma_v2sub_model_name = 'Bes_with_Sys/binom_slice_vars_model_dsigma_v2sub.csv'
    df_def_avgs_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_bes.csv'
    df_def_avgs_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_model.csv'
    df_def_avgs_v2sub_out_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_bes.csv'
    df_def_avgs_v2sub_out_model_name = 'Bes_with_Sys/dsig_tprotons_avgs_v2sub_model.csv'
    df_partitions_fits_name = 'Bes_with_Sys/partition_width_fits_bes.csv'
    df_partitions_fits_model_name = 'Bes_with_Sys/partition_width_fits_model.csv'

    save_dir = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Presentations/STAR_Paper/Autogen_Figures/'

    cent_map = {8: '0-5%', 7: '5-10%', 6: '10-20%', 5: '20-30%', 4: '30-40%', 3: '40-50%', 2: '50-60%', 1: '60-70%',
                0: '70-80%', -1: '80-90%'}

    stat_plot = 'k2'  # 'standard deviation', 'skewness', 'non-excess kurtosis'
    div_plt = 120
    exclude_divs = [89, 356]  # [60, 72, 89, 90, 180, 240, 270, 288, 300, 356]
    cent_plt = 8
    energies_fit = [7, 11, 19, 27, 39, 62]
    samples = 72  # For title purposes only

    data_sets_plt = ['bes_def', 'ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    data_sets_colors = dict(zip(data_sets_plt, ['black', 'red', 'blue', 'purple']))
    data_sets_labels = dict(zip(data_sets_plt, ['STAR', 'AMPT', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$']))
    data_sets_markers = dict(zip(data_sets_plt, [dict(zip(['raw', 'mix', 'diff'], [x, x, x]))
                                                 for x in ['o', 's', '^', '*']]))
    data_sets_bands = ['ampt_new_coal_epbins1', 'cf_resample_epbins1', 'cfev_resample_epbins1']
    legend_order = ['STAR', 'MUSIC+FIST', 'MUSIC+FIST EV $1fm^3$', 'AMPT']

    cent_ref_name = 'mean_cent_ref.csv'
    cent_ref_df = pd.read_csv(f'{base_path}{cent_ref_name}')
    ref_type = 'refn'  # 'refn'
    cent_ref_df = cent_ref_df.replace('bes_resample_def', 'bes_def')
    cent_ref_df = cent_ref_df.replace('ampt_new_coal_resample_def', 'ampt_new_coal_epbins1')

    df_path = base_path + df_name
    df = pd.read_csv(df_path)
    df = df.dropna()
    df['name'] = df['name'].str.replace('bes_sys_', '')
    all_sets = pd.unique(df['name'])
    print(all_sets)

    vz_sets = [x for x in all_sets if 'efficiency' in x]

    df = df[df['stat'] == stat_plot]

    # Calculate dsigma with k2 values and get systematics
    df = df[df['stat'] == 'k2']
    df = df.drop('stat', axis=1)

    # Get only datasets to plot in this short script to make nlo calculation faster
    df = df[df['energy'] == 39]
    df = df[df['cent'] == cent_plt]
    df = df[df['divs'] == div_plt]

    print('Calc dsigma')
    df_raw, df_mix, df_diff = calc_dsigma(df, ['raw', 'mix', 'diff'])
    df_dsigma_types = pd.concat([df_raw, df_mix, df_diff])
    print('Calc diff nlo error')
    df_dsigma_types = add_diff_nlo_err(df_dsigma_types, group_cols=['energy', 'cent', 'name', 'total_protons'],
                                       exclude_divs=[356, 89])

    dvar_vs_protons(df_dsigma_types, div_plt, cent_plt, [39], ['raw', 'mix', 'diff'], ['bes_def'],
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00085, 0.00055],
                    data_sets_labels={'bes_def': {'raw': 'Single Event', 'mix': 'Mixed Event',
                                                  'diff': 'Single Event - Mixed Event'}},
                    marker_map={'bes_def': {'raw': 'o', 'mix': 's', 'diff': '^'}})  # Mix subtraction demo

    dvar_vs_protons(df_dsigma_types, div_plt, cent_plt, [39], ['diff'], ['bes_def'] + vz_sets,
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00119, 0.00069], xlim=[4.5, 24.5])

    dvar_vs_protons(df_dsigma_types, div_plt, cent_plt, [39], ['raw'], ['bes_def'] + vz_sets,
                    plot=True, avg=False, alpha=1.0, y_ranges=[-0.00119, 0.00069], xlim=[4.5, 24.5])

    plt.show()


if __name__ == '__main__':
    main()
