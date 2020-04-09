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
    run_id = gen_id(sys.argv[1])
    print(f'Run_id: {run_id}')
    run(run_id)
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


def run(run_id):
    """
    Run AMPT with random seed. Once finished run makeAmptroot.C to create root file.
    :param run_id: Unique id for run generated with date/time and job number
    :return:
    """
    sp.run(['./ampt'], input=str(run_id).encode('utf-8'))
    sp.run(['root', '-b', '-q', 'makeAmptroot.C("' + str(run_id) + '")++'])


if __name__ == '__main__':
    main()
