#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 27 3:48 PM 2022
Created in PyCharm
Created as QGP_Scripts/job_manager

@author: Dylan Neff, Dyn04
"""

import os
import time


def main():
    """
    Run from /star/u/dneff/Ampt_Bad_Event/sub
    Submit jobs to find bad Ampt files and fix them
    :return:
    """
    # top_path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/'
    # file_list_path = '/star/u/dneff/Ampt_Bad_Event/sub/list/root_files.txt'
    # sub_path = '/star/u/dneff/git/QGP_Scripts/Ampt_Bad_Event/clean_sub.xml'
    # output_path = '/star/u/dneff/Ampt_Bad_Event/sub/output/'
    # user = 'dneff'
    # check_interval = 60

    pars = init_pars()

    files = get_files(pars['top_path'])
    submit_job(files, pars['file_list_path'], pars['sub_path'])
    babysit_job(files, pars)
    combine_outputs(pars[''])
    clean_up()
    print('donzo')


def init_pars():
    """
    Define parameters to be passed around and used in dictionary
    :return: Dictionay containing all parameters
    """
    pars = {
        'top_path': '/gpfs01/star/pwg/dneff/data/AMPT/most_central/',
        'file_list_path': '/star/u/dneff/Ampt_Bad_Event/sub/list/root_files.txt',
        'sub_path': '/star/u/dneff/git/QGP_Scripts/Ampt_Bad_Event/clean_sub.xml',
        'output_path': '/star/u/dneff/Ampt_Bad_Event/sub/output/',
        'user': 'dneff',
        'check_interval': 60,  # seconds

        'condor_flag_left': ':',
        'condor_flag_right': 'jobs',
        'out_split_flag': 'Files checked:\n',
    }

    return pars


def submit_job(files, file_list_path, sub_path):
    """
    Submit jobs to check files
    :param files: .root files to be checked
    :param file_list_path: Path to output file list of .root files found
    :param sub_path: Path to submission xml to be star-submit ed
    """
    write_file_list(files, file_list_path)
    submit_job(sub_path)


def babysit_job(files, pars):
    """
    Wait for initial submission to finish running and then check output to see if all files
    have been checked. If so exit, else resubmit job for missing files.
    :param files: .root file paths that need to be checked. Compare to those that were actually checked.
    :param pars: Dictionary of parameters
    :return:
    """
    finished = False
    while not finished:
        while check_jobs_alive(pars['user'], pars['condor_flag_left'], pars['condor_flag_right']):
            print(f'Jobs alive, waiting {pars["check_interval"]}s to check again')
            time.sleep(pars['check_interval'])

        files_checked = check_outputs(pars['output_dir'], pars['out_split_flag'])
        files_remaining = set(files) - set(files_checked)
        if len(files_remaining) > 0:
            print(f'Resubmitting {len(files_remaining)} missing files')
            submit_job(files_remaining, pars['file_list_path'], pars['sub_path'])
            time.sleep(pars['check_interval'] * 4)  # Wait a while to let submission go through before checking condor
        else:
            finished = True


def check_jobs_alive(user='dneff', flag_left=':', flag_right='jobs'):
    job_status = os.popen(f'condor_q {user} | tail -4').read()
    try:
        jobs_alive = int(job_status[job_status.find(flag_left) + len(flag_left):job_status.find(flag_right)].strip())
    except ValueError:
        print('Bad read of jobs alive!')

    return jobs_alive > 0


def check_outputs(output_dir, flag):
    files_checked = []
    for out_file_path in os.listdir(output_dir):
        with open(out_file_path, 'r') as out_file:
            files_checked.append(out_file.read().split(flag)[-1].strip().split('\n'))

    return files_checked


def combine_outputs(output_path, out_combo_path, flag):
    out_combo_lines = []
    for out_path in os.listdir(output_path):
        with open(out_path, 'r') as out_file:
            out_combo_lines.append(out_file.read().split(flag)[0].strip().split('\n'))

    with open(out_combo_path, 'r') as combo_file:
        combo_file.write('\n'.join(out_combo_lines))


def clean_up(output_path, out_combo_path):
    if input('Clean up?\n')[0].lower() == 'y':
        os.system('clean.sh')
        for out_path in os.lisdir(output_path):
            if out_path != out_combo_path:
                os.remove(out_path)


def write_file_list(files, file_list_path):
    """
    Write file list of all .root files in files
    :param files: .root files found to be written to file_list
    :param file_list_path: Path to output file list of .root files found
    """
    with open(file_list_path, 'w') as file_list:
        file_list.writelines('\n'.join(files))


def submit_job(sub_path):
    os.system(f'star-submit {sub_path}')


def get_files(path):
    """
    Get all root paths recursively from path
    :param path: Top level to find root paths
    :return root_paths: List of root file paths
    """
    root_paths = []
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            root_paths.append(root_path)

    return root_paths[:100]


if __name__ == '__main__':
    main()
