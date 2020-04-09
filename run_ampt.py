#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 09 12:59 PM 2020
Created in PyCharm
Created as QGP_Scripts/run_ampt.py

@author: Dylan Neff, dylan
"""

import subprocess as sp
import sys
import os
from time import sleep


def main():
    run_id = gen_id(sys.argv[1])
    # set_run_dir(run_id)
    print(f'Run_id: {run_id}')
    run(run_id)
    # clean_up(run_id)
    print('donzo')


def gen_id(job_id):
    run_id = str(job_id).split('_')[-1] + '1'
    return run_id


def set_run_dir(run_id):
    sp.run(['mkdir', str(run_id)])
    os.chdir(str(run_id))
    # sp.run(['cd', str(run_id)])
    sp.run(['mkdir', 'ana'])
    sp.run(['cp', '../input.ampt', 'ana/'])
    sp.run(['cp', '../input.ampt', '.'])
    sp.run(['cp', '../ampt', '.'])
    sp.run(['cp', '../makeAmptroot.C', '.'])


def run(run_id):
    # p = sp.Popen(['./ampt'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.STDOUT)
    # sleep(5)
    # p.communicate(input=str(run_id).encode('utf-8'))
    # p.wait()
    sp.run(['./ampt'], input=str(run_id).encode('utf-8'))
    sp.run(['root', '-b', '-q', 'makeAmptroot.C("' + str(run_id) + '")++'])
    # sp.run(['mv', 'ana/test.root', f'../test_{run_id}.root'])


def clean_up(run_id):
    # sp.run(['cd', '..'])
    os.chdir('..')
    sp.run(['rm', '-r', str(run_id)])


if __name__ == '__main__':
    main()
