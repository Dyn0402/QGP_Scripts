#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 8:29 PM 2019
Created in PyCharm
Created as QGP_Scripts/sum_events

@author: Dylan Neff, dylan
"""


def main():
    energy = 7
    all_events, all_files = read_file(f'/home/dylan/Research/rcf_data_paths/{energy}GeVall.txt')
    avail_events, avail_files = read_file(f'/home/dylan/Research/rcf_data_paths/{energy}GeVavail.txt')
    print(f'{all_events}\t{all_files}\t{avail_events}\t{avail_files}')
    # comp_files()

    print('donzo')


def read_file(path):
    events = 0
    files = 0
    max_event = 0
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            num = int(line.split('::')[-1])
            events += num
            if num > max_event:
                max_event = num
            files += 1
    file_name = path.split('/')[-1]
    print(f'File: {file_name}  |  Total events: {events}  |  Total files: {files}  |  Max Events: {max_event}')
    return events, files


def comp_files():
    path1 = '/home/dylan/Research/rcf_data_paths/7GeV.txt';
    path2 = '/home/dylan/Research/rcf_data_paths/7GeV4.txt'
    events = 0
    files = 0
    files1 = []
    files2 = []
    with open(path1, 'r') as file:
        lines = file.readlines()
        for line in lines:
            files1.append([line.split('::')[-2], line.split('::')[-1]])
    with open(path2, 'r') as file:
        lines = file.readlines()
        for line in lines:
            files2.append([line.split('::')[-2], line.split('::')[-1]])
    for file in files2:
        if file not in files1:
            print(f'New file from local: {file}')


if __name__ == '__main__':
    main()
