#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 04 8:31 AM 2019
Created in PyCharm
Created as QGP_Scripts/resubmit.py

@author: Dylan Neff, dylan
"""


import os
import sys


def main():
    resub()
    print('donzo')


def local_test():
    script_list = get_script_list('/home/dylan/Research/Tree_Maker_Logs/script/3GeV/')
    out_list = get_out_list('/home/dylan/Research/Tree_Maker_Logs/script/3GeV/')
    failed_jobs = get_failed_jobs(script_list, out_list)
    failed_jobs = get_failed_jobs(script_list, out_list)
    resub_flag = ask_to_resub(failed_jobs)
    if resub_flag:
        resub_jobs('/home/dylan/Research/Tree_Maker_Logs/script/3GeV/', failed_jobs)


def resub():
    try:
        energy = sys.argv[1]
        ref = sys.argv[2]
    except IndexError:
        print('Need to input energy and ref as command line arguments!')
        return
    script_path = '/gpfs01/star/pwg/dneff/scratch/trees_ref' + str(ref) + '/script/' + str(energy) + 'GeV/'
    output_path = '/gpfs01/star/pwg/dneff/scratch/trees_ref' + str(ref) + '/output/' + str(energy) + 'GeV/'
    err_path = '/gpfs01/star/pwg/dneff/scratch/trees_ref' + str(ref) + '/log/' + str(energy) + 'GeV/'
    break_list = get_break_list(err_path)
    script_list = get_script_list(script_path)
    out_list = get_out_list(output_path)
    failed_jobs = get_failed_jobs(script_list, out_list)
    resub_flag = ask_to_resub(failed_jobs, energy)
    if resub_flag:
        resub_jobs(script_path, failed_jobs)


def get_script_list(path):
    """
    Get .csh files in directory and return the jobid plus job number.
    Tested, seems to work.
    """
    script_list = []
    for fpath in os.listdir(path):
        if '.csh' in fpath:
            script_list.append(fpath[5:-4])

    return script_list


def get_out_list(path):
    """
        Get .csh files in directory and return the jobid plus job number.
        Tested, seems to work.
        """
    out_list = []
    for fpath in os.listdir(path):
        if '.root' in fpath:
            out_list.append(fpath[5:-5])

    return out_list


def get_failed_jobs(script_list, out_list):
    """
    Compare script job submission files to output root files to determine which jobs have failed.
    Tested, seems to work.
    """
    failed_jobs = []
    for job in script_list:
        if job not in out_list:
            failed_jobs.append(job)

    return failed_jobs


def ask_to_resub(incomplete_jobs, energy):
    """
    Display failed jobs and ask user if they should be resubmitted.
    Tested, seems to work.
    """
    if len(incomplete_jobs) <= 0:
        print(f'\nAll jobs for {energy}GeV completed!\n')
        return False
    print(f'\n{len(incomplete_jobs)} failed jobs:')
    for job in incomplete_jobs:
        print(job)
    res = input('\nResubmit failed jobs listed above? Enter "yes" to resubmit and anything else to quit: ')
    if res.strip().lower() == 'yes':
        return True
    else:
        return False


def resub_jobs(script_path, failed_jobs):
    """
    Split condor files and resubmit individual failed jobs.
    """
    split_condor(script_path)
    for job in failed_jobs:
        command = 'condor_submit ' + script_path + 'sched' + job + '.condor'
        os.system(command)


def get_break_list(path):
    breaks = 0
    files = 0
    break_list = []
    for fpath in os.listdir(path):
        if '.err' in fpath:
            files += 1
            with open(path + fpath, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if '*** Break *** segmentation violation' in line:
                        breaks += 1
                        break_list.append(fpath)
                        break
    print(f'Files: {files}  |  Breaks: {breaks}  |  Percentage Broken: {float(breaks) / files * 100}%')

    return break_list


def split_condor(path):
    """
    Split any condor files in directory into individual files.
    Tested, seems to work.
    """
    content_lines = 15
    sep_lines = 1
    for fpath in os.listdir(path):
        if '.condor' in fpath:
            with open(path+fpath, 'r') as file:
                try:
                    sched_id, job_num_low, job_num_high = fpath.split('_')
                except ValueError:
                    continue
                job_num_high = int(job_num_high.split('.')[0])
                job_num_low = int(job_num_low)
                job_num = job_num_high
                lines = file.readlines()
                while job_num >= job_num_low:
                    with open(path+sched_id+'_'+str(job_num)+'.condor', 'w') as new_file:
                        line_low = (job_num_high - job_num) * (content_lines + sep_lines)
                        line_high = (job_num_high - job_num + 1) * (content_lines + sep_lines)
                        # print(path+sched_id+'_'+str(job_num)+'.condor')
                        # print()
                        # print(lines[line_low:line_high])
                        # print()
                        new_file.writelines(lines[line_low:line_high])
                        job_num -= 1


if __name__ == '__main__':
    main()
