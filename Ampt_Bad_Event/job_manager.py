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

from remove_bad_event import fix_dataset


def main():
    """
    Run from /star/u/dneff/Ampt_Bad_Event/sub
    Submit jobs to find bad Ampt files and fix them
    :return:
    """
    pars = init_pars()

    files = get_files(pars['top_path'])
    res = input('Submit initial jobs before check?  ')
    if len(res) > 0 and res[0].lower() == 'y':
        submit_jobs(files, pars['file_list_path'], pars['sub_path'])
    babysit_jobs(files, pars)
    combine_outputs(pars['output_path'], pars['output_combo_path'], pars['out_split_flag'], pars['list_path'])

    res = input('Fix dataset?  ')
    if len(res) > 0 and res[0].lower() == 'y':
        cwd = os.getcwd()
        os.chdir(pars['fix_tree_cpp_dir'])
        fix_dataset(pars['output_combo_path'], pars['bad_repo_path'],
                    pars['bad_sufx'], pars['fix_sufx'], pars['min_identical'], False)
        os.chdir(cwd)

    clean_up()
    print('donzo')


def init_pars():
    """
    Define parameters to be passed around and used in dictionary
    :return: Dictionay containing all parameters
    """
    pars = {
        'top_path': '/gpfs01/star/pwg/dneff/data/AMPT/min_bias/',
        'file_list_path': '/star/u/dneff/Ampt_Bad_Event/sub/list/root_files.txt',
        'sub_path': '/star/u/dneff/git/QGP_Scripts/Ampt_Bad_Event/clean_sub.xml',
        'output_path': '/star/u/dneff/Ampt_Bad_Event/sub/output/',
        'result_path': '/star/u/dneff/Ampt_Bad_Event/sub/result/',
        'list_path': '/star/u/dneff/Ampt_Bad_Event/sub/list/',
        'log_path': '/star/u/dneff/Ampt_Bad_Event/sub/log/',
        'term_log_path': '/star/u/dneff/Ampt_Bad_Event/sub/log/resubed_terms/',
        'output_combo_path': '/star/u/dneff/Ampt_Bad_Event/sub/result/ampt_bad_events.txt',
        'fix_tree_cpp_dir': '/star/u/dneff/git/QGP_Scripts/Ampt_Bad_Event/',
        'bad_repo_path': '/gpfs01/star/pwg/dneff/data/AMPT/bad_events/run_dir/',
        'user': 'dneff',
        'check_interval': 120,  # seconds

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
            print(f' {now} | run time: {now - start}  Waiting {pars["check_interval"]}s to check again.')
            print(f'  {jobs_alive} jobs alive:  ' + ', '.join([f'{num} {cat}' for cat, num in job_status.items()]))
            if jobs_alive <= 0:
                break
            terminated_jobs = check_terminated(pars['log_path'])
            if len(terminated_jobs) > 0:
                terminated_files = get_job_files(terminated_jobs, pars['list_path'])
                print(f'\n\nResubmitting {len(terminated_files)} files from {len(terminated_jobs)} terminated jobs '
                      f'{terminated_jobs}\n')
                submit_jobs(terminated_files, pars['file_list_path'], pars['sub_path'])
                move_log_files(terminated_jobs, pars['log_path'], pars['term_log_path'])
            time.sleep(pars['check_interval'])

        files_checked = check_outputs(pars['output_path'], pars['out_split_flag'])
        files_checked = convert_files(files_checked, pars['list_path'])
        files_remaining = list(set(files) - set(files_checked))
        print('files_checked:  ', len(files_checked), '\nfiles expected:  ', len(files),
              '\nfiles_remaining:  ', len(files_remaining))
        if len(files_remaining) > 0:
            print(f'First 5 files_checked:\n{files_checked[:5]}\nFirst 5 files expected:\n{files[:5]}\n'
                  f'First 5 files_remaining:\n{files_remaining[:5]}')
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
            new_files = out_file.read().split(flag)[-1].strip().split('\n')
            new_files = [x.strip() for x in new_files]
            files_checked.extend(new_files)

    return files_checked


def combine_outputs(output_path, out_combo_path, flag, list_path):
    out_combo_lines = []
    for out_path in os.listdir(output_path):
        with open(output_path + out_path, 'r') as out_file:
            out_read = out_file.read()
        bad_events = out_read.split(flag)[0].strip().split('\n')
        for event in bad_events:
            temp_path = [x.removeprefix('path: ') for x in event.split('\t') if 'path: ' in x]
            if len(temp_path) == 1:
                real_path = convert_files(temp_path, list_path)[0]
                out_combo_lines.append(event.replace(temp_path[0], real_path))

    out_combo_lines = [*set(out_combo_lines)]  # Remove duplicates

    with open(out_combo_path, 'w') as combo_file:
        combo_file.write('\n'.join(out_combo_lines))


def clean_up():
    if input('Clean up?  ')[0].lower() == 'y':
        os.system('/star/u/dneff/Ampt_Bad_Event/sub/clean.sh')


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

    return root_paths


def convert_files(temp_files, list_path):
    """
    A bit complicated but should be reliable unless list files are messed with
    :param temp_files:
    :param list_path:
    :return:
    """
    list_names = os.listdir(list_path)

    return_files = []
    for temp_file in temp_files:
        file_name = temp_file.split('/')[-1].strip()
        job_name = temp_file.split('/')[-3]
        list_match = [x for x in list_names if job_name + '.' in x]  # without '.', job aaa2 grabs aaa21, aaa22, etc

        if len(list_match) == 0:
            print(f'Can\'t find list for {temp_file}!')
        elif len(list_match) > 1:
            print(f'Duplicate lists? : {list_match} for {temp_file}')
        else:
            list_match = list_match[0]

        with open(list_path + list_match, 'r') as list_file:
            list_lines = list_file.readlines()
        real_path = [line.strip() for line in list_lines if file_name in line]

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


def check_terminated(log_path):
    terminated_jobs = []
    for log_name in os.listdir(log_path):
        if log_name[-4:] == '.err':
            with open(log_path + log_name, 'r') as err_file:
                err_lines = err_file.readlines()
                if len(err_lines) > 0 and 'Terminated' in err_lines[0]:
                    terminated_jobs.append(log_name.strip('.err').strip('err_'))

    return terminated_jobs


def move_log_files(jobs, log_path, move_path):
    log_names = os.listdir(log_path)
    for job in jobs:
        for log_name in log_names:
            if job in log_name:
                try:
                    os.rename(log_path + log_name, move_path + log_name)
                except FileNotFoundError:
                    print(f'File not found? {log_path + log_name} to {move_path + log_name}')


def get_job_files(jobs, list_path):
    list_names = os.listdir(list_path)
    files = []
    for job in jobs:
        for list_name in list_names:
            if job + '.' in list_name:  # without '.', job aaa2 grabs aaa21, aaa22, etc
                with open(list_path + list_name, 'r') as list_file:
                    files.extend(list_file.read().strip().split('\n'))

    return files


if __name__ == '__main__':
    main()
