#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 12:21 PM 2019
Created in PyCharm
Created as QGP_Scripts/list_missing.py

@author: Dylan Neff, dylan
"""


import os

# scp -r dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/scratch/trees/output/62GeV/ Trees/


def main():
    expected = 3747
    numbers = list(range(0, expected))
    prefix = ''
    files = os.listdir('/home/dylan/Research/Trees/62GeV/')
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

    print(numbers)

    print('donzo')


if __name__ == '__main__':
    main()
