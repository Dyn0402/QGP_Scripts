#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 12 6:59 PM 2023
Created in PyCharm
Created as QGP_Scripts/check_files_identical

@author: Dylan Neff, Dylan
"""

import os


def main():
    base_path = 'F:/Research/Data/default_rand_sys_test/'
    dir1 = 'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/7GeV/'
    dir2 = 'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_2/7GeV/'
    # dir2 = 'rapid05_resample_norotate_strefnoseed_nomix_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0/7GeV/'
    verbose = False

    file_names = [fname for fname in os.listdir(f'{base_path}{dir1}') if '.txt' in fname]

    num_diff = 0
    for file_name in file_names:
        if os.path.exists(f'{base_path}{dir2}{file_name}'):
            with open(f'{base_path}{dir1}{file_name}', 'r') as file1:
                lines1 = file1.readlines()
            with open(f'{base_path}{dir2}{file_name}', 'r') as file2:
                lines2 = file2.readlines()

            if verbose:
                print()
                print(file_name)
            if len(lines1) != len(lines2):
                print(f'Unequal number of lines')
                num_diff += 1
            else:
                files_same = True
                for line_num in range(len(lines1)):
                    if lines1[line_num] != lines2[line_num]:
                        print(f'Lines not same in {file_name}: ')
                        print(lines1[line_num])
                        print(lines2[line_num])
                        print()
                        files_same = False
                if files_same:
                    if verbose:
                        print('Files exactly the same!')
                else:
                    num_diff += 1

    print('\nBetween directories:')
    print(dir1)
    print(dir2)
    print(f'{num_diff} of {len(file_names)} files different')

    print('donzo')


if __name__ == '__main__':
    main()
