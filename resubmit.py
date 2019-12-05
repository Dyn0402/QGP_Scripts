#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 04 8:31 AM 2019
Created in PyCharm
Created as QGP_Scripts/resubmit.py

@author: Dylan Neff, dylan
"""


import os
import sys


def main():
    energy = sys.argv[1]
    script_path = '/gpfs01/star/pwg/dneff/scratch/trees/script/' + str(energy) + 'GeV/'
    output_path = '/gpfs01/star/pwg/dneff/scratch/trees/output/' + str(energy) + 'GeV/'
    err_path = '/gpfs01/star/pwg/dneff/scratch/trees/log/' + str(energy) + 'GeV/'
    get_break_list(err_path)
    script_list = get_script_list(script_path)
    out_list = get_out_list(output_path)

    print('donzo')


def get_script_list(path):
    pass


def get_out_list(path):
    pass


def get_break_list(path):
    breaks = 0
    files = 0
    break_list = []
    for fpath in os.listdir(path):
        if '.err' in fpath:
            files += 1
            with open(path + fpath, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if '*** Break *** segmentation violation' in line:
                        breaks += 1
                        break_list.append(fpath)
                        break
    print(f'Files: {files}  |  Breaks: {breaks}  |  Percentage Broken: {float(breaks) / files * 100}%')

    return break_list


if __name__ == '__main__':
    main()
