#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 19 12:55 PM 2020
Created in PyCharm
Created as QGP_Scripts/cent_def_defaults.py

@author: Dylan Neff, dylan
"""

import numpy as np


def main():
    path = '/home/dylan/git/Research/QGP_Fluctuations/Tree_Reader/StRoot/StRefMultCorr/Centrality_def_refmult2.txt'
    energies = {7: [], 11: [], 19: [], 27: [], 39: [], 62: []}
    full_string = 'map<int, map<int, unsigned>> refnn {'
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            line = line.strip().split()
            energy = int(float(line[1]))
            if energy in energies.keys():
                ref_list = []
                for col in line[6:22]:
                    ref_list.append(int(col))
                print(f'energy {energy}: {ref_list}')
                if len(energies[energy]) > 0:
                    if not ref_list == energies[energy]:
                        print(f'{energy} different lists')
                else:
                    energies[energy] = ref_list

    for energy in sorted(energies.keys()):
        ref_list_9bin = []
        skip = False
        for ref in energies[energy][:-2]:
            if skip:
                skip = False
            else:
                ref_list_9bin.append(ref)
                skip = True
        ref_list_9bin.append(energies[energy][-2])
        ref_list_9bin.append(energies[energy][-1])
        ref_diff = [i - j for i, j in zip(ref_list_9bin[1:], ref_list_9bin[:-1])]
        ref_mid = [int(ref_list_9bin[i] + j / 2) for i, j in enumerate(ref_diff)]
        ref_mid.insert(0, 1)
        ref_mid.append(int(ref_list_9bin[-1] * 1.2))

        print(f'{energy}: {ref_list_9bin}  len = {len(ref_list_9bin)}')
        out_string = f'{{{energy}, {{'
        for i, j in enumerate(ref_mid):
            out_string += f'{{{i - 1}, {j}}}, '
        out_string = out_string[:-2] + '}}, '

        print(out_string)
        full_string += out_string

    print()
    print(full_string[:-2] + '};')
    print()
    print('donzo')


if __name__ == '__main__':
    main()
