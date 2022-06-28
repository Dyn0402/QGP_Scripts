#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 27 3:48 PM 2022
Created in PyCharm
Created as QGP_Scripts/job_manager

@author: Dylan Neff, Dyn04
"""

import os


def main():
    """
    Run from /star/u/dneff/Ampt_Bad_Event/sub
    Submit jobs to find bad Ampt files and fix them
    :return:
    """
    top_path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/'
    file_list_path = '/star/u/dneff/Ampt_Bad_Event/sub/list/root_files.txt'
    sub_path = '/star/u/dneff/git/QGP_Scripts/Ampt_Bad_Event/clean_sub.xml'
    init_sub(top_path, file_list_path, sub_path)
    print('donzo')


def init_sub(top_path, file_list_path, sub_path):
    gen_file_list(top_path, file_list_path)
    submit_job(sub_path)


def gen_file_list(top_path, file_list_path):
    files = get_files(top_path)
    with open(file_list_path, 'w') as file_list:
        file_list.writelines('\n'.join(files))


def submit_job(sub_path):
    os.system(f'star-submit {sub_path}')


def get_files(path):
    """
    Get all root paths recursively from path
    :param path: Top level to find root paths
    :return: List of root file paths
    """
    root_paths = []
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            root_paths.append(root_path)

    return root_paths[:10]


if __name__ == '__main__':
    main()
