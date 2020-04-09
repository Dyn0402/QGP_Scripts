#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 09 12:59 PM 2020
Created in PyCharm
Created as QGP_Scripts/run_ampt.py

@author: Dylan Neff, dylan
"""

import subprocess
import sys
import os


def main():
    run_id = gen_id(sys.argv[1])
    set_run_dir(run_id)
    print(f'Run_id: {run_id}')
    run(run_id)
    clean_up(run_id)
    print('donzo')


def gen_id(job_id):
    run_id = str(job_id).split('_')[-1] + '1'
    return run_id


def set_run_dir(run_id):
    subprocess.call(['mkdir', str(run_id)])
    os.chdir(str(run_id))
    # subprocess.call(['cd', str(run_id)])
    subprocess.call(['mkdir', 'ana'])
    subprocess.call(['cp', '../input.ampt', 'ana/'])
    subprocess.call(['cp', '../input.ampt', '.'])
    subprocess.call(['cp', '../ampt', '.'])
    subprocess.call(['cp', '../makeAmptroot.C', '.'])


def run(run_id):
    subprocess.call(['./ampt', str(run_id)])
    subprocess.call(['root', '-b', '-q', 'makeAmptroot.C++'])
    subprocess.call(['mv', 'ana/test.root', f'../test_{run_id}.root'])


def clean_up(run_id):
    # subprocess.call(['cd', '..'])
    os.chdir('..')
    subprocess.call(['rm', '-r', str(run_id)])


if __name__ == '__main__':
    main()
