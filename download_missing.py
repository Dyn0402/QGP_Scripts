#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 12:21 PM 2019
Created in PyCharm
Created as QGP_Scripts/download_missing.py

@author: Dylan Neff, dylan
"""


import os
import subprocess

# scp -r dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/scratch/trees/output/62GeV/ Trees/


def main():
    download_all()

    print('donzo')


def download_single(energy, ref):
    energy_files = {7: 1685, 11: 568, 19: 1397, 27: 2036, 39: 3580, 62: 3747}
    energy = 19
    ref = 3
    expected = energy_files[energy]
    numbers = list(range(0, expected))
    prefix = ''
    path = f'/media/dylan/SSD_Storage/Research/Trees_Ref{ref}/{energy}GeV/'
    files = os.listdir(path)
    for file in files:
        last_uscore = file.rfind('_')
        last_period = file.rfind('.')
        num = int(file[last_uscore + 1:last_period])
        if num in numbers:
            numbers.pop(next(index for index in range(len(numbers)) if numbers[index] == num))
        pre = file[:last_uscore + 1]
        if pre != prefix:
            print(f'New prefix {pre}')
            prefix = pre

    print(f'\n{len(numbers)} files missing:')
    missing_files = []
    for num in numbers:
        missing_files.append(f'{prefix}{num}.root')
        print(missing_files[-1])

    if len(numbers) > 0:
        res = input(f'\nDownload {len(numbers)} missing files? Enter yes to download or anything else to quit: ')
        if res.strip().lower() in ['yes', 'y']:
            for file in missing_files:
                command = f'scp dneff@rftpexp.rhic.bnl.gov:' \
                          f'/gpfs01/star/pwg/dneff/scratch/trees_ref/output/{energy}GeV/{file} ' \
                          f'{path}'
                print(command)
                os.system(command)


def download_all():
    remote_path = 'dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/scratch/'
    local_path = '/media/dylan/SSD_Storage/Research/'
    remote_tree_prefix = 'trees_old_ref'
    local_tree_prefix = 'Trees_Old_Ref'
    refs = [2, 3]
    energies = [7, 11, 19, 27, 39, 62]
    energy_files = {7: 1685, 11: 568, 19: 1397, 27: 2036, 39: 3580, 62: 3747}
    missing_files = {}
    total_missing = 0
    for ref in refs:
        missing_files.update({ref: {}})
        for energy in energies:
            missing_files[ref].update({energy: []})
            expected = energy_files[energy]
            numbers = list(range(0, expected))
            prefix = ''
            path = local_path + local_tree_prefix + f'{ref}/{energy}GeV'
            # path = f'/media/dylan/SSD_Storage/Research/Trees_Ref{ref}/{energy}GeV/'
            files = os.listdir(path)
            for file in files:
                last_uscore = file.rfind('_')
                last_period = file.rfind('.')
                num = int(file[last_uscore + 1:last_period])
                if num in numbers:
                    numbers.pop(next(index for index in range(len(numbers)) if numbers[index] == num))
                pre = file[:last_uscore + 1]
                if pre != prefix:
                    # print(f'New prefix {pre}')
                    prefix = pre

            print(f'Ref{ref} {energy}GeV \n{len(numbers)} files missing:')
            for num in numbers:
                missing_files[ref][energy].append(f'{prefix}{num}.root')
                total_missing += 1
                # print(missing_files[-1])

    if total_missing > 0:
        res = input(f'\nDownload {total_missing} missing files? Enter yes to download or anything else to quit: ')
        if res.strip().lower() in ['yes', 'y']:
            for ref in missing_files:
                for energy in missing_files[ref]:
                    if len(missing_files[ref][energy]) > 0:
                        local = local_path + local_tree_prefix + f'{ref}/{energy}GeV/'
                        # local = f'/media/dylan/SSD_Storage/Research/Trees_Ref{ref}/{energy}GeV/'
                        files = r'\{'
                        for file in missing_files[ref][energy]:
                            # remote = remote_path + remote_tree_prefix + f'{ref}/output/{energy}GeV/{file}'
                            # print(remote, local)
                            # subprocess.Popen(['gnome-terminal', '--', 'scp', remote, local])
                            files += file + ','
                        files = files[:-1] + r'\}'
                        remote = remote_path + remote_tree_prefix + f'{ref}/output/{energy}GeV/{files}'
                        # # command = f'scp dneff@rftpexp.rhic.bnl.gov:' \
                        # #           f'/gpfs01/star/pwg/dneff/scratch/trees_ref{ref}/output/{energy}GeV/{files} ' \
                        # #           f'{local}'
                        # # print(command)
                        # # os.system(command)
                        command = 'scp ' + remote + ' ' + local
                        print(command)
                        # subprocess.Popen(['gnome-terminal', '--', 'scp', remote, local])
                        subprocess.Popen(command, shell=True)

    else:
        print('All files downloaded!')


if __name__ == '__main__':
    main()
