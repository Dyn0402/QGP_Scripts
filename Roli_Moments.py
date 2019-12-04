#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 29 2:28 PM 2019
Created in PyCharm
Created as QGP_Scripts/Roli_Moments_Systematics.py

@author: Dylan Neff, dylan
"""

import matplotlib.pyplot as plt
import numpy as np


def main():
    path = '/home/dylan/local_server/dyn0402/Research/Roli_Mix_Systematics/moments/'
    file_pre_mix = 'mixed11_'
    file_pre_data = 'data11_'
    file_suf = '.txt'
    file_list_data = [path + file_pre_data + str(i) + file_suf for i in range(9)]
    file_list_mix = [path+file_pre_mix+str(i)+file_suf for i in range(9)]

    div_data = {}
    plot_data = [[], [], []]
    for data_path, mix_path in zip(file_list_data, file_list_mix):
        kurtosis = {'data': {}, 'mix': {}}
        for path, kind in zip((data_path, mix_path), ('data', 'mix')):
            file = open(path)
            lines = file.readlines()
            for line in lines:
                line = line.strip().split()
                cent, kurt, kurt_err = int(line[1]), float(line[8]), float(line[9])
                print(cent, kurt, kurt_err)
                kurtosis[kind].update({cent: (kurt, kurt_err)})

        print('here')
        for cent in kurtosis['data']:
            print('here2')
            val = kurtosis['data'][cent][0] / kurtosis['mix'][cent][0]
            err = abs(val) * ((kurtosis['data'][cent][1]/kurtosis['data'][cent][0])**2 +
                              (kurtosis['mix'][cent][1]/kurtosis['mix'][cent][0])**2)**0.5
            plot_data[0].append(cent)
            plot_data[1].append(val)
            plot_data[2].append(err)
            if cent in div_data:
                div_data[cent].append((val, err))
            else:
                div_data.update({cent: [(val, err)]})

    print(plot_data)

    plt.errorbar(plot_data[0], plot_data[1], yerr=plot_data[2], fmt='o', color='black',
             ecolor='lightgray', elinewidth=3, capsize=0)
    plt.grid()
    plt.xlabel('Centrality Bin')
    plt.title('11GeV Data Divided by Mixed Kurtosis (Roli\'s code)')
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
