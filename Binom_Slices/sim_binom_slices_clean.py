#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 01 2:10 PM 2021
Created in PyCharm
Created as QGP_Scripts/sim_binom_slices_clean.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.optimize import curve_fit as cf
from scipy.interpolate import interp1d

from AzimuthBinData import AzimuthBinData
from DistStats import DistStats
from Measure import Measure


def main():
    pars = init_pars()
    get_data(pars)


def init_pars():
    pars = {'datasets': define_datasets(),
            'base_path': 'D:/Research/'  # '/home/dylan/Research/'
            }

    return pars


def define_datasets():
    """
    Define sets to read for each dataset. Each dataset can contain multiple sets which will be combined to give
    systematic uncertainties. Here set parameters to define which sets will be read and combined.
    :return: list of dictionaries containing parameters for defining datasets
    """

    par_names = ['name', 'base_ext', 'exact_keys', 'contain_keys']
    par_vals = [['ampt_def', '_Ampt', ['default'], []],
                ['sim_amp01_spread01', '_Sim', ['anticlmulti', 'amp05', 'spread01'], ['single']]]

    datasets = [dict(zip(par_names, dset)) for dset in par_vals]

    return datasets


def get_data(pars):
    df = []
    for dataset in pars['datasets']:
        print(dataset)
        get_dataset(dataset, pars)


def get_dataset(dataset, pars):
    raw_path = pars['base_path'] + 'Data' + dataset['base_ext'] + '/'
    mix_path = pars['base_path'] + 'Data' + dataset['base_ext'] + '_Mix' + '/'

    for d in os.listdir(raw_path):
        good = True
        if not all(x in d.split('_') for x in dataset['exact_keys']):
            good = False
        if not all(x in d for x in dataset['contain_keys']):
            good = False
        if good:
            print(d)


# def get_dist_stats(base_path, divs, cents, sim_pars, stats):
#     sim_energy = 62
#     ampt_energy = 7
#
#     ampt_set = 'default/Ampt_rapid05_n1ratios_'
#
#     sim_dfs = []
#     ampt_dfs = []
#     for div in divs:
#         for cent in cents:
#             for pars in sim_pars:
#                 print(f'Get sim {pars}, {div} div, {cent} cent')
#                 sim_dist = get_sim_dists(base_path + 'Data_Sim/', range(11), div, cent, sim_energy,
#                                          exact_keys=['anticlmulti', *pars], contain_keys=['single'])
#                 sim_dist_mix = get_sim_dists(base_path + 'Data_Sim_Mix/', range(11), div, cent, sim_energy,
#                                              exact_keys=['anticlmulti', *pars], contain_keys=['single'])
#                 stats_dict = get_stats(sim_dist, sim_dist_mix, stats)
#                 pars_dict = {'div': div, 'cent': cent, 'energy': sim_energy, 'amp': pars[0], 'spread': pars[1],
#                              'name': f'sim_{pars[0]}_{pars[1]}'}
#                 sim_dfs.append(pd.DataFrame([{**entry, **pars_dict} for entry in stats_dict]))
#
#             print(f'Get ampt, {div} div, {cent} cent')
#             ampt_dist = get_ampt_dist(base_path + 'Data_Ampt/' + ampt_set, div, cent, ampt_energy, range(60))
#             ampt_dist_mix = get_ampt_dist(base_path + 'Data_Ampt_Mix/' + ampt_set, div, cent, ampt_energy, range(60))
#             stats_dict = get_stats(ampt_dist, ampt_dist_mix, stats)
#             pars_dict = {'div': div, 'cent': cent, 'energy': ampt_energy, 'amp': 'ampt', 'spread': 'ampt',
#                          'name': 'ampt'}
#             ampt_dfs.append(pd.DataFrame([{**entry, **pars_dict} for entry in stats_dict]))
#
#     df = pd.concat(sim_dfs, ignore_index=True)
#     df = df.append(pd.concat(ampt_dfs, ignore_index=True))
#
#     return df


if __name__ == '__main__':
    main()
