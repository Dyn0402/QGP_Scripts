#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 12:21 PM 2019
Created in PyCharm
Created as QGP_Scripts/download_missing.py

@author: Dylan Neff, dylan
"""


import os

# scp -r dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/scratch/trees/output/62GeV/ Trees/


def main():
    energy = 19
    expected = 1397
    numbers = list(range(0, expected))
    prefix = ''
    path = f'/media/dylan/SSD_Storage/Research/Trees_Ref3/{energy}GeV/'
    files = os.listdir(path)
    for file in files:
        last_uscore = file.rfind('_')
        last_period = file.rfind('.')
        num = int(file[last_uscore+1:last_period])
        if num in numbers:
            numbers.pop(next(index for index in range(len(numbers)) if numbers[index] == num))
        pre = file[:last_uscore+1]
        if pre != prefix:
            print(f'New prefix {pre}')
            prefix = pre

    print(f'\n{len(numbers)} files missing:')
    missing_files = []
    for num in numbers:
        missing_files.append(f'{prefix}{num}.root')
        print(missing_files[-1])

    if len(numbers) > 0:
        res = input('\nDownload missing? Enter yes to download or anything else to quit: ')
        if res.strip().lower() in ['yes', 'y']:
            for file in missing_files:
                command = f'scp dneff@rftpexp.rhic.bnl.gov:' \
                          f'/gpfs01/star/pwg/dneff/scratch/trees/output/{energy}GeV/{file} ' \
                          f'/home/dylan/Research/Trees/{energy}GeV/'
                print(command)
                os.system(command)

    print('donzo')


if __name__ == '__main__':
    main()
