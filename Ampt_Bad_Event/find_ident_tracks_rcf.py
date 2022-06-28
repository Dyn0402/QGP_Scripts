#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 27 5:43 PM 2022
Created in PyCharm
Created as QGP_Scripts/find_ident_tracks_rcf

@author: Dylan Neff, Dyn04
"""

import sys
from datetime import datetime

from find_ampt_identical_tracks_uproot import check_file


def main():
    if len(sys.argv) == 2:
        uproot_finder_rcf(sys.argv[1])
    else:
        print('No input root files!')
    print('donzo')


def uproot_finder_rcf(root_path_list):
    start = datetime.now()
    print(f'Start {start}\n')

    out_file_path = 'out.txt'
    tree_name = 'tree'
    write_mode = 'w'
    max_eta = 1
    ignore_pids = [313, 111]
    track_attributes = ['pid', 'px', 'py', 'pz']

    with open(root_path_list, 'r') as file:
        root_paths = file.readlines()

    jobs_start = datetime.now()
    print(f'Jobs Start {jobs_start}\n')
    bad_trees = []
    for root_path in root_paths:
        bad_trees.append(check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids))

    with open(out_file_path, write_mode) as file:
        for bad_tree in bad_trees:
            for bad_event in bad_tree:
                if bad_event is not None:
                    out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
                    file.write(f'{out_line}\n')
        file.write('\n\nFiles checked:\n')
        file.writelines('\n'.join(root_path))

    end = datetime.now()
    print(f'\nEnd {end}')
    print(f'Run time: {end - start}')
    print(f'Job run time: {end - jobs_start}')


if __name__ == '__main__':
    main()
