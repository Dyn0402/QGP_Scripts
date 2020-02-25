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
import re


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
    # break_list = get_break_list(err_path)
    print('Reading err files for status: ')
    status_lists = get_err_status(err_path)
    print('Checking script files: ')
    script_list = get_script_list(script_path)
    print('Checking output Root files: ')
    out_list = get_out_list(output_path)
    print('Comparing Root files to script files to check for failures: ')
    failed_jobs = get_failed_jobs(script_list, out_list)
    print('Cross checking script/Root file comparison with err file status: ')
    cross_check_failed(failed_jobs, status_lists)
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
            fpath = fpath.strip('.csh').strip('sched')
            script_list.append(fpath)

    return script_list


def get_out_list(path):
    """
        Get .csh files in directory and return the jobid plus job number.
        Tested, seems to work.
        """
    out_list = []
    for fpath in os.listdir(path):
        if '.root' in fpath:
            fpath = fpath.strip('.root').strip('auau_')
            out_list.append(fpath)

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


def get_err_status(path):
    files = 0
    breaks = []
    finished = []
    terminated = []
    running = []
    for fpath in os.listdir(path):
        if '.err' in fpath:
            files += 1
            with open(path + fpath, 'r') as file:
                job = re.sub(r'err_\d+GeV', '', fpath).strip('.err')
                lines = file.readlines()
                alive = True
                for line in lines:
                    if '*** Break *** segmentation violation' in line:
                        breaks.append(job)
                        alive = False
                        break
                if alive and ' Terminated ' in line[0]:
                    terminated.append(job)
                    alive = False
                if alive and ' Done ' in line:
                    finished.append(job)
                    alive = False
                if alive:
                    running.append(job)
    print(f'Files: {files}  |  Breaks: {len(breaks)}  |  Percentage Broken: {float(len(breaks)) / files * 100}%')
    print(f'Files: {files}  |  Terminated: {len(terminated)}  |  Percentage Terminated: '
          f'{float(len(terminated)) / files * 100}%')
    print(f'Files: {files}  |  Finished: {len(finished)}  |  Percentage Finished: '
          f'{float(len(finished)) / files * 100}%')
    print(f'Files: {files}  |  Running: {len(running)}  |  Percentage Running: '
          f'{float(len(running)) / files * 100}%')
    if len(breaks) + len(terminated) + len(finished) + len(running) != files:
        print(f'Bad accounting. {len(breaks) + len(terminated) + len(finished) + len(running)} '
              f'err files categorized out of {files}')

    return {'breaks': breaks, 'terminated': terminated, 'finished': finished, 'running': running}


def cross_check_failed(failed_list, status_lists):
    fails_remaining = []
    for fail_path in failed_list:
        if fail_path in status_lists['breaks']:
            status_lists['breaks'].remove(fail_path)
        elif fail_path in status_lists['terminated']:
            status_lists['terminated'].remove(fail_path)
        elif fail_path in status_lists['running']:
            status_lists['running'].remove(fail_path)
        else:
            fails_remaining.append(fail_path)

    if len(fails_remaining) > 0:
        print(f'No Root file for following but err file says it is finished: ')
        for fail in fails_remaining:
            print(fail)
    if len(status_lists['breaks']) > 0:
        print(f'Root file exists for following but err file says it seg faulted: ')
        for fail in status_lists['breaks']:
            print(fail)
    if len(status_lists['terminated']) > 0:
        print(f'Root file exists for following but err file says it was terminated: ')
        for fail in status_lists['terminated']:
            print(fail)
    if len(status_lists['running']) > 0:
        print(f'Root file exists for following but err file says it is still running: ')
        for fail in status_lists['running']:
            print(fail)


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
