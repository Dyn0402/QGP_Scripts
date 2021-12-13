#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 12 12:53 PM 2021
Created in PyCharm
Created as QGP_Scripts/resample_time_plot

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def main():
    path = '/home/dylan/Desktop/Old_PC_Tree_Reader_Break_Output.txt'
    df = []
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            if len(line) == 11 and line[2] == 'complete':
                energy = int(line[0].strip('GeV'))
                progress = float(line[1].strip('%'))
                time = float(line[7].strip('s'))
                df.append({'energy': energy, 'progress': progress, 'time': time})
    df = pd.DataFrame(df)
    print(df[df.energy == 7])
    fig, ax = plt.subplots()
    for energy in np.unique(df.energy):
        df_energy = df[df.energy == energy]
        ax.scatter(df_energy.time / 60 / 60, df_energy.progress, label=f'{energy}GeV', marker='o', alpha=0.8)
    ax.legend()
    ax.grid()
    ax.set_xlabel('Run Time (hrs)')
    ax.set_ylabel('Progress (%)')
    fig.tight_layout()
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
