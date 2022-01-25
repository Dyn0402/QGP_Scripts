#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 21 4:07 PM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_dat_single_file.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

import os


def main():
    dats_path = f'D:/Research/AMPT_Trees/ref_fix3_dat_most_central/string_melting/7GeV/'
    pid_select = -2212

    px_list, py_list, pz_list, rapid_list = [], [], [], []
    for dat in os.listdir(dats_path)[:1]:
        with open(dats_path + dat, 'r') as file:
            lines = file.readlines()
            index = 0
            num_lines = len(lines)
            while index < num_lines:
                n_particles = int(lines[index].strip().split()[2])
                index += 1
                event_particle_count = 0
                while event_particle_count < n_particles:
                    if index >= num_lines:
                        break
                    pid, px, py, pz, m, x, y, z, t = lines[index].strip().split()
                    pid = int(pid)
                    if pid == pid_select:
                        px, py, pz, m = [float(x) for x in (px, py, pz, m)]
                        rapid_list.append(rapidity(px, py, pz, m))
                        px_list.append(px)
                        py_list.append(py)
                        pz_list.append(pz)
                    event_particle_count += 1
                    index += 1

    for name, p_list in zip(('px', 'py', 'pz', 'rapidity'), (px_list, py_list, pz_list, rapid_list)):
        fig, ax = plt.subplots()
        ax.hist(p_list, bins=100)
        ax.set_xlabel(name)
    plt.show()
    print('donzo')


def rapidity(px, py, pz, m):
    e = np.sqrt(m**2 + px**2 + py**2 + pz**2)
    return np.arctanh(pz / e)


if __name__ == '__main__':
    main()
