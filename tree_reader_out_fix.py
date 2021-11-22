#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 04 11:21 AM 2021
Created in PyCharm
Created as QGP_Scripts/tree_reader_out_fix.py

@author: Dylan Neff, dylan
"""


import os
import shutil


def main():
    check_path = '/home/dylan/Research/Data_Mix/default_sys/'
    dirs = os.listdir(check_path)
    for d in dirs:
        if d[-2:] == '00':
            if d[:-1] in dirs:
                print(f'Found correspondent for : {d}')
                for d2 in os.listdir(os.path.join(check_path, d)):
                    shutil.move(os.path.join(check_path, d, d2), os.path.join(check_path, d[:-1]))
            else:
                print(f'Didn\'t find correspondent for : {d}')

    print('donzo')


if __name__ == '__main__':
    main()
