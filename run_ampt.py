#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 09 12:59 PM 2020
Created in PyCharm
Created as QGP_Scripts/run_ampt.py

@author: Dylan Neff, dylan

Run AMPT event generator and convert output to filtered root tree.
"""

import subprocess as sp
import sys
from datetime import datetime


def main():
    times = [f'Python start: {str(datetime.now())}']
    print(sys.argv)
    run_id = gen_id(sys.argv[1])
    if len(sys.argv) >= 3:
        run(run_id, times, sys.argv[2])
    else:
        run(run_id, times)
    print(f'Run_id: {run_id}')
    for time in times:
        print(time)
    print('donzo')


def gen_id(job_id):
    """
    Generate unique id for ampt run based on date/time and job id.
    Used as random seed for AMPT as well as output root file name.
    :param job_id: Condor job id input from command line
    :return: Unique id
    """
    max_ampt_num = 2100000001  # Maximum random seed tested that works on AMPT
    now = datetime.now()
    seconds_into_day = now.hour * 3600 + now.minute * 60 + now.second
    run_id = str(seconds_into_day).zfill(5)
    run_id += str(job_id).split('_')[-1] + '1'
    while int(run_id) > max_ampt_num:
        run_id = run_id[1:]

    return run_id


def run(run_id, times, tscript_suf=''):
    """
    Run AMPT with random seed. Once finished run makeAmptroot.C to create root file.
    :param run_id: Unique id for run generated with date/time and job number
    :param times: list of times to print out at end
    :param tscript_suf: suffix on root script, either '' for STAR BES1 cuts or 'all' for no cuts
    :return:
    """
    times.append(f'AMPT start: {str(datetime.now())}')
    sp.run(['./ampt'], input=str(run_id).encode('utf-8'))
    times.append(f'AMPT end ROOT start: {str(datetime.now())}')
    sp.run(['root', '-b', '-q', 'makeAmptroot' + tscript_suf + '.C++("' + str(run_id) + '")'])
    times.append(f'ROOT end: {str(datetime.now())}')


if __name__ == '__main__':
    main()
