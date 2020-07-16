#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 16 2:25 PM 2020
Created in PyCharm
Created as QGP_Scripts/Brian_Counterpoint.py

@author: Dylan Neff, dylan
"""

kB_per_line = 15.0 / 1000
lines_per_file = 1000
line_content = 'Brian is wrong\n'
file_name_base = 'Brian_Wrong_'
file_extension = '.txt'


def main():
    path = '/home/dylan/Desktop/Test/'
    make_files(path, 1)
    print('donzo')


def make_file(base, name):
    with open(base+name, 'w') as file:
        for i in range(lines_per_file):
            file.write(line_content)


def make_files(base, size):
    """
    Make files at base path
    :param base: Base path to write files
    :param size: Total size of written data requested in GB
    :return:
    """
    num_files = int(size * 1e6 / kB_per_line / lines_per_file)
    print(f'{num_files} files to get {size}Gb')
    for file_num in range(num_files):
        make_file(base, file_name_base+str(file_num)+file_extension)


if __name__ == '__main__':
    main()
