#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 20 10:37 AM 2021
Created in PyCharm
Created as QGP_Scripts/fix_resample_proton_shift

@author: Dylan Neff, dylan
"""

import os


def main():
    """
    Due to bug all total protons were shifted by 1 (except default raw dists). Go through resample files and correct.
    :return:
    """
    # base_path = '/home/dylan/Research/'
    base_path = 'D:/Transfer/Research/'
    raw_data_paths = ['Data', 'Data_Ampt', 'Data_Sim']
    mix_data_paths = ['Data_Mix', 'Data_Ampt_Mix', 'Data_Sim_Mix']
    total_proton_shift = -1
    fix = True
    fix_time = 1640293741  # seconds since epoch for Thu Dec 23 2021. Only fix if not modified after

    n_fix_files = 0
    for data_type, paths in zip(['raw', 'mix'], [raw_data_paths, mix_data_paths]):
        for data_path in paths:
            data_full_path = f'{base_path}{data_path}/'
            for set_dir in os.listdir(data_full_path):
                if '_resample' in set_dir:
                    set_path = f'{data_full_path}{set_dir}/'
                    print(set_path)
                    for subset_dir in os.listdir(set_path):
                        subset_path = f'{set_path}{subset_dir}/'
                        for energy_dir in os.listdir(subset_path):
                            energy_path = f'{subset_path}{energy_dir}/'
                            for file_name in os.listdir(energy_path):
                                if 'ratios_divisions_' in file_name:
                                    file_path = f'{energy_path}{file_name}'
                                    if os.path.getmtime(file_path) < fix_time:
                                        # print(file_path)
                                        n_fix_files += 1
                                        if fix:
                                            fix_file(file_path, data_type, total_proton_shift)

    print(f'Number of fix files: {n_fix_files}')
    print('donzo')


def fix_file(file_path, file_type, total_proton_shift=-1):
    """
    Shift all total protons (except default raw dist) and rewrite file.
    :param file_path: Path of file to fix
    :param file_type: Either 'raw' or 'mix'. Leave default raw distribution alone
    :param total_proton_shift: How much to shift total protons by. -1 for the bug
    :return:
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        replace = False
        if file_type == 'mix':
            replace = True
        for i_line in range(len(lines)):
            if replace:
                line = lines[i_line]
                line_split = line.split('\t')
                if len(line_split) == 2:
                    line_split[0] = int(line_split[0]) + total_proton_shift
                    lines[i_line] = f'{line_split[0]}\t{line_split[1]}'
            else:
                if 'bootstrap #' in lines[i_line]:
                    replace = True
                else:
                    continue

    with open(file_path, 'w') as file:
        file.writelines(lines)


if __name__ == '__main__':
    main()
