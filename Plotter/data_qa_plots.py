#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 05 4:42 PM 2022
Created in PyCharm
Created as QGP_Scripts/data_qa_plots

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import uproot
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # read_trees()
    read_qas_uproot()
    print('donzo')


def read_trees():
    base_path = 'F:/Research/BES1_Trees/'
    energies = [39]
    hist_name = ''
    track_attributes = ['refmult', 'refmult2', 'refmult3', 'btof_multi', 'btof_match', 'vx', 'vy', 'vz',
                        'proton.dca', 'proton.pt', 'proton.nsigma', 'proton.rapid', 'proton.charge',
                        'proton.nhits_fit', 'proton.m2']
    max_rapid = 0.5
    charge = 1
    max_dca = 1.0
    min_beta = 0.0
    min_nhits_fit = 20.0
    min_m2, max_m2 = 0.6, 1.2
    min_pt_notof, max_pt_notof = 0.4, 0.8
    min_pt_tof, max_pt_tof = 0.8, 2.0
    max_p_notof, max_p_tof = 1.0, 3.0
    fig_ref3_nproton, ax = plt.subplots(dpi=144)
    ref3_nproton = np.histogram2d([], [], bins=(np.arange(-0.5, 800.5, 1), np.arange(-0.5, 100.5, 1)))
    for energy in energies:
        max_nsigma = 2 if energy != 27 else 1
        trees_path = f'{base_path}{energy}GeV/'
        for file_name in os.listdir(trees_path):
            with uproot.open(f'{trees_path}{file_name}') as file:
                tree_index = 1
                while True:
                    tree_name = f'tree;{tree_index}'
                    if tree_name in file:
                        print(file.classnames())
                        print(file[tree_name].arrays(track_attributes))
                        input()
                        # tracks = file[tree_name].arrays(track_attributes)
                        # tracks = tracks[(abs(tracks['proton.rapid']) < max_rapid) &
                        #                 (tracks['proton.charge'] == charge) &
                        #                 (abs(tracks['proton.nsigma']) < max_nsigma) &
                        #                 (tracks['proton.dca'] < max_dca) & (tracks['proton.nhits_fit'] > min_nhits_fit)]
                        # tof_tracks = tracks[(tracks['proton.pt'] >= min_pt_tof) &
                        #                     (tracks['proton.pt'] <= max_pt_tof) &
                        #                     (tracks['proton.pt'] <= max_p_tof)]
                        # tof_tracks = tof_tracks[(tof_tracks['proton.beta'] > min_beta) &
                        #                         (tof_tracks['proton.m2'] > min_m2) & (tof_tracks['proton.m2'] < max_m2)]
                        # tpc_tracks = tracks[(tracks['proton.pt'] >= min_pt_notof) &
                        #                     (tracks['proton.pt'] <= max_pt_notof) &
                        #                     (tracks['proton.p'] <= max_p_notof)]
                        #
                    else:
                        break
    print('donzo')


# def read_qas():
#     base_path = '/media/ucla/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
#     energies = [7]
#     hist_name = ''
#     fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
#     for energy in energies:
#         qa_path = f'{base_path}{energy}GeV/QA_{energy}GeV.root'
#         file = ROOT.TFile(qa_path, 'READ')
#         file.Get(f'pre_pt_rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0_{energy};1').Draw()
#         input()


def read_qas_uproot():
    data_set = 'BES'
    # data_set = 'AMPT_New_Coal'
    base_paths = {'BES': 'F:/Research/Data/default/'
                         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/',
                  'AMPT_New_Coal': 'F:/Research/Data_Ampt_New_Coal/default_resample_epbins1/'
                                   'Ampt_rapid05_resample_norotate_epbins1_0/'}
    save_paths = {'BES': 'F:/Research/Results/BES_QA_Plots/',
                  'AMPT_New_Coal': 'F:/Research/Results/AMPT_New_Coal_QA_Plots/'}
    base_path = base_paths[data_set]
    save_path = save_paths[data_set]
    energies = [7, 11, 19, 27, 39, 62]
    plot_pannels = (2, 3)
    pre_post_colors = {'pre': 'blue', 'post': 'red'}
    name_convert = {'rapid': 'rapidity', 'nsigma': 'nsigma_proton', 'eta': 'pseudorapidity', 'run': 'run_number',
                    'ref': 'refmult', 'refn': 'refmult3', 'ep': 'event_plane_angle'}

    th1_names = get_thn_pre_posts(base_path, energies[0], 1)
    # th2_names = get_thn_pre_posts(base_path, energies[0], 2)
    print(th1_names)

    figs, axs = {}, {}
    for name in th1_names:
        fig, ax = plt.subplots(*plot_pannels, sharex=True, sharey=False, dpi=144, figsize=(13.33, 6))
        title = name[:name.find('_')] if name[:4] != 'btof' else name[:name.find('_', 5)]
        if title in name_convert:
            title = name_convert[title]
        # title = name_convert[name[:name.find('_')]]
        for ax_i in ax.flat[-plot_pannels[1]:]:
            ax_i.set_xlabel(title)
        fig.suptitle(f'{data_set} {title}')
        fig.canvas.manager.set_window_title(title)
        # fig.subplots_adjust(hspace=0)
        figs.update({name: fig})
        axs.update({name: ax.flat})

        fig, ax = plt.subplots(dpi=144, figsize=(6, 5))
        ax.set_xlabel(title)
        fig.suptitle(f'{data_set} {title} pre cuts')
        fig.canvas.manager.set_window_title(title)
        # fig.subplots_adjust(hspace=0)
        figs.update({f'{name}_pre_cut': fig})
        axs.update({f'{name}_pre_cut': ax})

        fig, ax = plt.subplots(dpi=144, figsize=(6, 5))
        ax.set_xlabel(title)
        fig.suptitle(f'{data_set} {title} post cuts')
        fig.canvas.manager.set_window_title(title)
        # fig.subplots_adjust(hspace=0)
        figs.update({f'{name}_post_cut': fig})
        axs.update({f'{name}_post_cut': ax})

    for energy_index, energy in enumerate(energies):
        qa_path = f'{base_path}{energy}GeV/QA_{energy}GeV.root'
        with uproot.open(qa_path) as file:
            for name in th1_names:
                for pre_post, color in pre_post_colors.items():
                    hist = file[f'{pre_post}_{name}{energy};1'].to_numpy()
                    vals = hist[0]
                    bin_centers = (hist[1][1:] + hist[1][:-1]) / 2
                    bin_widths = hist[1][1:] - hist[1][:-1]
                    # axs[name][energy_index].bar(x=bin_centers, height=vals, width=bin_widths, ls='steps', color=color)
                    axs[name][energy_index].step(bin_centers, vals, color=color, alpha=0.7, label=f'{pre_post} cuts')
                    if 'phi_' in name and data_set == 'AMPT_New_Coal' and pre_post == 'post':
                        axs[name][energy_index].set_ylim(0, 1.1 * max(vals))
                    # sns.histplot(x=bin_centers, y=vals, binwidth=bin_widths, ax=axs[name][energy_index])
                    # sns.histplot(y=vals, bins=hist[1], ax=axs[name][energy_index])
                    if pre_post == 'pre':
                        density = vals / (np.sum(vals) * (bin_centers[-1] - bin_centers[0]))
                        axs[f'{name}_pre_cut'].step(bin_centers, density, alpha=0.7, label=f'{energy}GeV')
                    elif pre_post == 'post':
                        density = vals / (np.sum(vals) * (bin_centers[-1] - bin_centers[0]))
                        axs[f'{name}_post_cut'].step(bin_centers, density, alpha=0.7, label=f'{energy}GeV')
                if energy_index == 1:
                    axs[name][energy_index].legend()

    for name in th1_names:
        for pre_post in ['pre', 'post']:
            axs[f'{name}_{pre_post}_cut'].axhline(0, color='black')
            axs[f'{name}_{pre_post}_cut'].legend()
            figs[f'{name}_{pre_post}_cut'].tight_layout()
            figs[f'{name}_{pre_post}_cut'].savefig(f'{save_path}{figs[name]._suptitle.get_text()}_{pre_post}_cut.png')
            figs[f'{name}_{pre_post}_cut'].savefig(f'{save_path}{figs[name]._suptitle.get_text()}_{pre_post}_cut.pdf')
        for energy_index, energy in enumerate(energies):
            axs[name][energy_index].axhline(0, color='black')
            axs[name][energy_index].text(0.5, 0.5, f'{energy} GeV', alpha=0.05, horizontalalignment='center',
                                         fontsize=60,
                                         verticalalignment='center', transform=axs[name][energy_index].transAxes)
            figs[name].tight_layout()
            figs[name].savefig(f'{save_path}{figs[name]._suptitle.get_text()}.png')
            figs[name].savefig(f'{save_path}{figs[name]._suptitle.get_text()}.pdf')

    # plt.show()


def get_thn_pre_posts(base_path, energy, thn=1):
    tfn_name_list = []
    qa_path = f'{base_path}{energy}GeV/QA_{energy}GeV.root'
    with uproot.open(qa_path) as file:
        print(file.classnames())
        for key, val in file.classnames().items():
            if val[:3] == f'TH{thn}' and key[:4] == 'pre_':
                tfn_name_list.append(key[4:key.rfind('_') + 1])

    return tfn_name_list


if __name__ == '__main__':
    main()
