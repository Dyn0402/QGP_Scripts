#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 04 7:13 AM 2019
Created in PyCharm
Created as QGP_Scripts/find_seg_faults.py

@author: Dylan Neff, dylan
"""

import os


def main():
    path = '/home/dylan/Research/Tree_Maker_Logs/log/7GeV/'
    breaks = 0
    files = 0
    for fpath in os.listdir(path):
        # print(path)
        if '.err' in fpath:
            files += 1
            with open(path+fpath, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if '*** Break *** segmentation violation' in line:
                        breaks += 1
                        print(fpath)
                        break
    print(f'Files: {files}  |  Breaks: {breaks}  |  Percentage Broken: {float(breaks)/files*100}%')
    print('donzo')


if __name__ == '__main__':
    main()
