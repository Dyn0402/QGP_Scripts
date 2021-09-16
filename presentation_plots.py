#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 15 4:56 PM 2021
Created in PyCharm
Created as QGP_Scripts/presentation_plots

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
from Analyzer.AzimuthBinData import AzimuthBinData


def main():
    divs = 120
    cent = 8
    energy = 7
    set_group = 'default'
    set_name = 'Ampt_rapid05_n1ratios_'
    set_num = 8
    data_set = '_Ampt'
    base_path = '/home/dylan/Research/Data'
    set_path = f'{set_group}/{set_name}{set_num}/{energy}GeV/ratios_divisions_{divs}_centrality_{cent}_local.txt'
    path = f'{base_path}{data_set}/{set_path}'
    path_mix = f'{base_path}{data_set}_Mix/{set_path}'
    raw = AzimuthBinData(path=path, div=divs).get_dist()
    mix = AzimuthBinData(path=path_mix, div=divs).get_dist()

    plot_ratio(raw, mix)

    print('donzo')


def plot_ratio(raw, mix):
    pass


if __name__ == '__main__':
    main()
