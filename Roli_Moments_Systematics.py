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
    file_pre = 'mixed11_'
    file_suf = '.txt'
    file_list = [path+file_pre+str(i)+file_suf for i in range(9)]

    cent = []
    kurtosis = []
    kurt_err = []
    for file_path in file_list:
        file = open(file_path)
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            cent.append(int(line[1]))
            kurtosis.append(float(line[8]))
            kurt_err.append(float(line[9]))

    grouped = {}
    for cen, kurt, err in zip(cent, kurtosis, kurt_err):
        if cen in grouped:
            grouped[cen].append((kurt, err))
        else:
            grouped.update({cen: [(kurt, err)]})

    sigmas = {key: np.std(list(zip(*value))[0]) for (key, value) in grouped.items()}

    prop = {key: np.sqrt(np.sum(np.asarray(list(zip(*value))[1])**2))/len(grouped) for (key, value) in grouped.items()}

    cent_data = []
    kurt_data = []
    kurt_err_data = []
    file_list = [path + 'data11_' + str(i) + file_suf for i in range(9)]
    for file_path in file_list:
        file = open(file_path)
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            cent_data.append(int(line[1]))
            kurt_data.append(float(line[8]))
            kurt_err_data.append(float(line[9]))

    # plt.style.use('seaborn-whitegrid')
    plt.errorbar(cent, kurtosis, yerr=kurt_err, fmt='o', color='blue', ecolor='lightskyblue', elinewidth=3,
                 capsize=0, label='Mixed')
    plt.errorbar(cent_data, kurt_data, yerr=kurt_err_data, fmt='o', color='red', ecolor='salmon', elinewidth=3,
                 capsize=0, label='Data')
    plt.grid()
    plt.legend(loc='upper left')
    plt.xlabel('Centrality Bin')
    plt.title('11GeV Kurtosis (Roli\'s code)')
    plt.show()

    sig_x, sig_y = zip(*sigmas.items())
    print(sig_x)
    print(sig_y)

    prop_x, prop_y = zip(*prop.items())
    print(prop_x)
    print(prop_y)

    plt.plot(sig_x, sig_y, 'o')
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
