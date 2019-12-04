#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 8:29 PM 2019
Created in PyCharm
Created as QGP_Scripts/sum_events

@author: Dylan Neff, dylan
"""


def main():
    read_file('/home/dylan/Research/rcf_data_paths/7GeV.txt')
    read_file('/home/dylan/Research/rcf_data_paths/7GeV4.txt')

    print('donzo')


def read_file(path):
    events = 0
    files = 0
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            num = int(line.split('::')[-1])
            events += num
            files += 1
    file_name = path.split('/')[-1]
    print(f'File: {file_name}  |  Total events: {events}  |  Total files: {files}')


if __name__ == '__main__':
    main()
