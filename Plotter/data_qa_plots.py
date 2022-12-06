#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 05 4:42 PM 2022
Created in PyCharm
Created as QGP_Scripts/data_qa_plots

@author: Dylan Neff, Dylan
"""

import os

import uproot
import ROOT
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # read_trees()
    # read_qas()
    read_qas_uproot()
    print('donzo')


def read_trees():
    base_path = 'F:/Research/BES1_Trees/'
    energies = [7]
    hist_name = ''
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
    for energy in energies:
        trees_path = f'{base_path}{energy}GeV/'
        for file_name in os.listdir(trees_path):
            with uproot.open(f'{trees_path}{file_name}') as file:
                print(file.classnames())


def read_qas():
    base_path = '/media/ucla/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
    energies = [7]
    hist_name = ''
    fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)
    for energy in energies:
        qa_path = f'{base_path}{energy}GeV/QA_{energy}GeV.root'
        file = ROOT.TFile(qa_path, 'READ')
        file.Get(f'pre_pt_rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0_{energy};1').Draw()
        input()


def read_qas_uproot():
    base_path = '/media/ucla/Research/Data/default_resample/rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0/'
    energies = [7, 11, 19, 27, 39, 62]
    plot_pannels = (2, 3)
    pre_post_colors = {'pre': 'blue', 'post': 'red'}
    name_convert = {'rapid': 'rapidity', 'nsigma': 'nsigma_proton', 'eta': 'pseudorapidity', 'run': 'run_number',
                    'ref': 'refmult', 'refn': 'refmult3', 'ep': 'event_plane_angle'}
    save_path = '/media/ucla/Research/Results/BES_QA_Plots/'

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
        fig.suptitle(title)
        fig.canvas.manager.set_window_title(title)
        # fig.subplots_adjust(hspace=0)
        figs.update({name: fig})
        axs.update({name: ax.flat})

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
                    # sns.histplot(x=bin_centers, y=vals, binwidth=bin_widths, ax=axs[name][energy_index])
                    # sns.histplot(y=vals, bins=hist[1], ax=axs[name][energy_index])
                if energy_index == 1:
                    axs[name][energy_index].legend()

    for name in th1_names:
        for energy_index, energy in enumerate(energies):
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
