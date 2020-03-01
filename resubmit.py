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
    print('Reading err files for status: ')
    status_lists = get_err_status(err_path)
    script_list = get_script_list(script_path)
    out_list = get_out_list(output_path)
    failed_jobs = get_failed_jobs(script_list, out_list)
    cross_check_failed(failed_jobs, status_lists)
    resub_flag, resub_set = ask_to_resub(status_lists['terminated']+status_lists['breaks'], failed_jobs, energy)
    if resub_flag:
        resub_jobs(script_path, resub_set)


def get_script_list(path):
    """
    Get .csh files in directory and return the jobid plus job number.
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
    """
    failed_jobs = []
    for job in script_list:
        if job not in out_list:
            failed_jobs.append(job)

    return failed_jobs


def ask_to_resub(incomplete_jobs, missing_jobs, energy):
    """
    Display failed jobs and ask user if they should be resubmitted.
    """
    resub_set = []
    if len(incomplete_jobs) == 0 and len(missing_jobs) == 0:
        print(f'\nAll jobs for {energy}GeV completed!\n')
        resub_flag = False
    else:
        print(f'\n{len(missing_jobs)} missing jobs:')
        for job in missing_jobs:
            print(job)
        print(f'\n{len(incomplete_jobs)} stopped jobs:')
        for job in incomplete_jobs:
            print(job)
        res = input('\nResubmit stopped jobs listed above? '
                    '\nEnter "yes" to resubmit stopped jobs only, "missing" to resubmit missing jobs only, '
                    '"both" to resubmit both (without duplicates), and anything else to quit: ')
        if res.strip().lower() == 'yes':
            resub_flag = True
            resub_set = incomplete_jobs
        elif res.strip().lower() == 'missing':
            resub_flag = True
            resub_set = missing_jobs
        elif res.strip().lower() == 'both':
            resub_flag = True
            resub_set = incomplete_jobs.copy()
            resub_set.extend([x for x in missing_jobs if x not in incomplete_jobs])
        else:
            resub_flag = False

    return resub_flag, resub_set


def resub_jobs(script_path, failed_jobs):
    """
    Split condor files and resubmit individual failed jobs.
    """
    split_condor(script_path)
    for job in failed_jobs:
        command = 'condor_submit ' + script_path + 'sched' + job + '.condor'
        os.system(command)


def get_err_status(path):
    """
    Looks at each .err file and determines status of job based on contents. Sorts job strings into corresponding lists.
    :param path: Path to directory containing .err files
    :return: Lists of job names sorted corresponding to status determined of error file
    """
    files = 0
    breaks = []
    finished = []
    terminated = []
    running = []
    grep_write_err = []
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
                if len(lines) > 0:
                    first_line = 0
                    while first_line < len(lines) and '/bin/grep: write error' in lines[first_line]:
                        if first_line == 0:
                            grep_write_err.append(job)
                        first_line += 1
                    if first_line < len(lines):
                        if alive and ' Terminated ' in lines[first_line]:
                            terminated.append(job)
                            alive = False
                        if alive and ' Done ' in lines[first_line]:
                            finished.append(job)
                            alive = False
                if alive:
                    running.append(job)

    print('\nBreak Jobs:')
    for job in breaks:
        print(job)
    print('\ngrep write error Jobs:')
    for job in grep_write_err:
        print(job)
    print('\nTerminated Jobs:')
    for job in terminated:
        print(job)
    print('\nRunning Jobs:')
    for job in running:
        print(job)

    print(f'\nFiles: {files}  |  Breaks: {len(breaks)}  |  Percentage Broken: {float(len(breaks)) / files * 100}%')
    print(f'Files: {files}  |  Terminated: {len(terminated)}  |  Percentage Terminated: '
          f'{float(len(terminated)) / files * 100:.4f}%')
    print(f'Files: {files}  |  Finished: {len(finished)}  |  Percentage Finished: '
          f'{float(len(finished)) / files * 100:.4f}%')
    print(f'Files: {files}  |  Running: {len(running)}  |  Percentage Running: '
          f'{float(len(running)) / files * 100:.4f}%')
    print(f'Files: {files}  |  grep write error: {len(grep_write_err)}  |  Percentage grep write error: '
          f'{float(len(grep_write_err)) / files * 100:.4f}%')
    if len(breaks) + len(terminated) + len(finished) + len(running) != files:
        print(f'Bad accounting. {len(breaks) + len(terminated) + len(finished) + len(running)} '
              f'err files categorized out of {files}')

    return {'breaks': breaks, 'terminated': terminated, 'finished': finished, 'running': running,
            'grep_write_err': grep_write_err}


def cross_check_failed(failed_list, status_lists):
    """
    Cross check job statuses determined by .err files with statuses determined by comparing .csh and .root files.
    .csh/.root comparison only determines if the full process completes while .err files can determine more
    (running, terminated, broken, etc)
    :param failed_list: List of missing .root jobs with resepect to expected from .csh files
    :param status_lists: Lists of job status based on .err files
    :return: -
    """
    breaks = status_lists['breaks'].copy()
    terminated = status_lists['terminated'].copy()
    running = status_lists['running'].copy()
    fails_remaining = []
    for fail_path in failed_list:
        if fail_path in breaks:
            breaks.remove(fail_path)
        elif fail_path in terminated:
            terminated.remove(fail_path)
        elif fail_path in running:
            running.remove(fail_path)
        else:
            fails_remaining.append(fail_path)

    if len(fails_remaining) > 0:
        print(f'No Root file for following but err file says it is finished: ')
        for fail in fails_remaining:
            print(fail)
    if len(breaks) > 0:
        print(f'Root file exists for following but err file says it seg faulted: ')
        for fail in breaks:
            print(fail)
    if len(terminated) > 0:
        print(f'Root file exists for following but err file says it was terminated: ')
        for fail in terminated:
            print(fail)
    if len(running) > 0:
        print(f'Root file exists for following but err file says it is still running: ')
        for fail in running:
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
