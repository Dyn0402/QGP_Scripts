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
from datetime import datetime

from remove_bad_event_rcf import fix_dataset


def main():
    """
    Run from /star/u/dneff/Ampt_Bad_Event/sub
    Submit jobs to find bad Ampt files and fix them
    :return:
    """
    pars = init_pars()

    files = get_files(pars['top_path'])
    submit_jobs(files, pars['file_list_path'], pars['sub_path'])
    babysit_jobs(files, pars)
    combine_outputs(pars['output_path'], pars['output_combo_path'])

    fix_dataset(pars['output_combo_path'], pars['result_path'],
                pars['bad_sufx'], pars['fix_sufx'], pars['min_identical'], True)

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
        'result_path': '/star/u/dneff/Ampt_Bad_Event/sub/result/',
        'output_combo_path': '/star/u/dneff/Ampt_Bad_Event/sub/result/ampt_bad_events.txt',
        'user': 'dneff',
        'check_interval': 60,  # seconds

        'condor_flag_left': ':',
        'condor_flag_right': 'jobs',
        'out_split_flag': 'Files checked:\n',

        'bad_sufx': '_bad',
        'fix_sufx': '_fix',
        'min_identical': 2,
    }

    return pars


def submit_jobs(files, file_list_path, sub_path):
    """
    Submit jobs to check files
    :param files: .root files to be checked
    :param file_list_path: Path to output file list of .root files found
    :param sub_path: Path to submission xml to be star-submit ed
    """
    write_file_list(files, file_list_path)
    os.system(f'star-submit {sub_path}')


def babysit_jobs(files, pars):
    """
    Wait for initial submission to finish running and then check output to see if all files
    have been checked. If so exit, else resubmit job for missing files.
    :param files: .root file paths that need to be checked. Compare to those that were actually checked.
    :param pars: Dictionary of parameters
    :return:
    """
    start = datetime.now()
    print(f'\nJobs submitted, waiting {pars["check_interval"]}s to check them')
    time.sleep(pars['check_interval'])

    finished = False
    while not finished:
        jobs_alive, job_status = check_jobs_alive(pars['user'], pars['condor_flag_left'], pars['condor_flag_right'])
        while jobs_alive > 0:
            jobs_alive, job_status = check_jobs_alive(pars['user'], pars['condor_flag_left'], pars['condor_flag_right'])
            now = datetime.now()
            print(f' {jobs_alive} jobs alive, waiting {pars["check_interval"]}s to check again.')
            print(f'  Run time {now - start}  ' + ', '.join([f'{num} {cat}' for cat, num in job_status.items()]))
            time.sleep(pars['check_interval'])

        files_checked = check_outputs(pars['output_path'], pars['out_split_flag'])
        files_remaining = set(files) - set(files_checked)
        if len(files_remaining) > 0:
            print(f'Resubmitting {len(files_remaining)} missing files')
            submit_jobs(files_remaining, pars['file_list_path'], pars['sub_path'])
            time.sleep(pars['check_interval'] * 4)  # Wait a while to let submission go through before checking condor
        else:
            finished = True


def check_jobs_alive(user='dneff', flag_left=':', flag_right='jobs'):
    try:
        job_status = os.popen(f'condor_q {user} | tail -4').read()
        job_status.split('\n')
        job_status = [x for x in job_status.split('\n') if 'Total for query:' in x][0]
        categories = ['completed', 'removed', 'idle', 'running', 'held', 'suspended']
        job_status = {cat: int(x.strip(cat).strip) for cat in categories for x in job_status if cat in x}
        jobs_alive = sum(job_status.values())
    except ValueError:
        print('Bad condor job read!')
        jobs_alive, job_status = 1, {}

    return jobs_alive, job_status


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


def clean_up():
    if input('Clean up?\n')[0].lower() == 'y':
        os.system('clean.sh')


def write_file_list(files, file_list_path):
    """
    Write file list of all .root files in files
    :param files: .root files found to be written to file_list
    :param file_list_path: Path to output file list of .root files found
    """
    with open(file_list_path, 'w') as file_list:
        file_list.writelines('\n'.join(files))


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

    return root_paths[:10]


if __name__ == '__main__':
    main()
