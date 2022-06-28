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
    # submit_jobs(files, pars['file_list_path'], pars['sub_path'])
    # babysit_jobs(files, pars)
    combine_outputs(pars['output_path'], pars['output_combo_path'], pars['out_split_flag'], files, pars['list_path'])

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
        'list_path': '/star/u/dneff/Ampt_Bad_Event/sub/list/',
        'output_combo_path': '/star/u/dneff/Ampt_Bad_Event/sub/result/ampt_bad_events.txt',
        'user': 'dneff',
        'check_interval': 10,  # seconds

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
        while True:
            jobs_alive, job_status = check_jobs_alive(pars['user'])
            now = datetime.now()
            print(f' {now - start} run time, waiting {pars["check_interval"]}s to check again.')
            print(f'  {jobs_alive} jobs alive:  ' + ', '.join([f'{num} {cat}' for cat, num in job_status.items()]))
            if jobs_alive <= 0:
                break
            time.sleep(pars['check_interval'])

        files_checked = check_outputs(pars['output_path'], pars['out_split_flag'])
        files_checked = convert_files(files_checked, files, pars['list_path'])
        files_remaining = list(set(files) - set(files_checked))
        print('files_checked:\n', files_checked, '\nfiles:\n', files, '\nfiles_remaining:\n', files_remaining)
        if len(files_remaining) > 0:
            print(f'\n\nResubmitting {len(files_remaining)} missing files\n')
            submit_jobs(files_remaining, pars['file_list_path'], pars['sub_path'])
            print(f'\nJobs resubmitted, waiting {pars["check_interval"]}s to check them')
            time.sleep(pars['check_interval'])
        else:
            finished = True


def check_jobs_alive(user='dneff'):
    try:
        job_status = os.popen(f'condor_q {user} | tail -4').read()
        job_status = [x for x in job_status.split('\n') if 'Total for query:' in x][0]
        categories = ['completed', 'removed', 'idle', 'running', 'held', 'suspended']
        alive_categories = ['idle', 'running', 'held', 'suspended']
        job_status = [x.split(';')[-1].strip() for x in job_status.split(',')]
        job_status = {cat: int(x.split()[0].strip()) for x in job_status for cat in categories if cat in x}
        jobs_alive = sum([num for cat, num in job_status.items() if cat in alive_categories])
    except ValueError:
        print('Bad condor job read!')
        jobs_alive, job_status = 1, {}

    return jobs_alive, job_status


def check_outputs(output_dir, flag):
    files_checked = []
    for out_file_path in os.listdir(output_dir):
        with open(output_dir + out_file_path, 'r') as out_file:
            files_checked.extend(out_file.read().split(flag)[-1].strip().split('\n'))

    return files_checked


def combine_outputs(output_path, out_combo_path, flag, real_files, list_path):
    out_combo_lines = []
    for out_path in os.listdir(output_path):
        with open(output_path + out_path, 'r') as out_file:
            out_read = out_file.read()
        bad_events = out_read.split(flag)[0].strip().split('\n')
        for event in bad_events:
            temp_path = [x.strip('path: ') for x in event.split('\t') if 'path: ' in x]
            if len(temp_path) == 1:
                real_path = convert_files(temp_path, real_files, list_path)[0]
                out_combo_lines.append(event.replace(temp_path, real_path))

            # temp_files = [x for x in event.split('\t') if 'path: ' in x for event in bad_events]
            # temp_files = [x for x in temp_files]
            # temp_files = 1
            # out_combo_lines.extend(convert_files(temp_files, real_files, list_path))

    with open(out_combo_path, 'w') as combo_file:
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

    return root_paths[:20]


def convert_files(temp_files, real_files, list_path):
    """
    A bit complicated but should be reliable unless list files are messed with
    :param temp_files:
    :param real_files:
    :param list_path:
    :return:
    """
    list_names = os.listdir(list_path)

    return_files = []
    for temp_file in temp_files:
        file_name = temp_file.split('/')[-1]
        job_name = temp_file.split('/')[-3]
        list_match = [x for x in list_names if job_name in x]

        if len(list_match) == 0:
            print(f'Can\'t find list for {temp_file}!')
        elif len(list_match) > 1:
            print(f'Duplicate lists? : {list_match} for {temp_file}')
        else:
            list_match = list_match[0]

        with open(list_path + list_match, 'r') as list_file:
            list_lines = list_file.readlines()
        real_path = [line for line in list_lines if file_name in line]

        if len(real_path) == 0:
            print(f'Can\'t find file {temp_file} in {list_match}!')
        if len(real_path) > 1:
            print(f'Duplicate files? : {real_path}')
        else:
            real_path = real_path[0]

        return_files.append(real_path)
        # matches = [x for x in real_files if file_name in x]
        # assert(len(matches) == 1)
        # return_files.append(matches[0])

    return return_files


if __name__ == '__main__':
    main()
