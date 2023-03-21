#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 19 4:00 PM 2023
Created in PyCharm
Created as QGP_Scripts/flat_tests.py

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector


def main():
    main_path = 'F:/Research/BES1_Trees/39GeV/'
    hist_name = 'original_psi_west_cent_8_runkey_11106;1'
    entries, means, mean_errs, job_nums = [], [], [], []
    all_job_nums, job_events = [], []
    for file_name in os.listdir(main_path):
        if '_qa.root' in file_name:
            with uproot.open(main_path + file_name) as file:
                if hist_name in file.classnames():
                    hist = file[hist_name]
                    entries.append(np.sum(hist.values()))
                    mean = np.sum(hist.values() * hist.axis().centers()) / np.sum(hist.values())
                    sd = np.sqrt(np.sum(hist.values() * hist.axis().centers()**2) / np.sum(hist.values()) - mean**2)
                    if np.sum(hist.values()) > 200:
                        means.append(mean)
                        mean_errs.append(sd / np.sqrt(np.sum(hist.values())))
                        job_nums.append(int(file_name.split('_')[2]))
                    if np.sum(hist.values()) > 4000:
                        print(file_name, mean, sd)
                    all_job_nums.append(int(file_name.split('_')[2]))
                    job_events.append(np.sum(hist.values()))
                    # print(hist.values())
                    # print(hist.axis().centers())
                    # print(np.mean(hist.values() * hist.axis().centers()))
                    # hist.
        # print(file_name)
    fig, ax = plt.subplots()
    ax.hist(means)
    ax.axvline(np.pi / 2, color='red')

    fig2, ax2 = plt.subplots(figsize=(10, 9))
    ax2.grid()
    ax2.axvline(np.pi / 2, ls='--', color='red')
    ax2.errorbar(means, job_nums, xerr=mean_errs, ls='none', marker='o')
    ax2.set_xlabel('Mean')
    ax2.set_ylabel('Job Number')
    fig2.tight_layout()

    fig3, ax3 = plt.subplots()
    ax3.scatter(all_job_nums, job_events)
    ax3.set_xlabel('Job Number')
    ax3.set_ylabel('Number of Events in Job')

    plt.show()



    print('donzo')


if __name__ == '__main__':
    main()
