#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 16 8:40 PM 2023
Created in PyCharm
Created as QGP_Scripts/tree_read_output_crash_files.py

@author: Dylan Neff, dylan
"""


def main():
    print('here')
    output_file_path = '/home/dylan/Desktop/crash_output.txt'
    open_files = []
    with open(output_file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            start_end, file_num = get_file_start_end(line)
            if start_end == 'Start':
                if file_num in open_files:
                    print(f'{file_num} already open!')
                else:
                    open_files.append(file_num)
                    print(open_files)
            elif start_end == 'End':
                if file_num in open_files:
                    open_files.remove(file_num)
                else:
                    print(f'{file_num} not open!')
    print('donzo')


def get_file_start_end(line):
    if 'Starting ' in line:
        file_num = int(line.split(' ')[-1].split('_')[-1].replace('.root', ''))
        return 'Start', file_num
    elif 'Ending ' in line:
        file_num = int(line.split(' ')[-1].split('_')[-1].replace('.root', ''))
        return 'End', file_num
    else:
        return None, None


if __name__ == '__main__':
    main()
