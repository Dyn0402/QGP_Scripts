#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 08 9:52 PM 2021
Created in PyCharm
Created as QGP_Scripts/CFSampleRunner

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from time import sleep
import psutil


class Pars:
    def __init__(self):
        self.build_path = '/home/dylan/git/Research/CooperFryeSampler/build/'
        self.converter_path = '/home/dylan/git/Research/CooperFryeSampler/'
        self.temp_path = '/home/dylan/Research/CFSample_Trees/temp'
        self.out_dir = '/home/dylan/Research/CFSample_Trees/'
        self.sampler_name = 'CooperFryeSampler'
        self.ram_buffer = 2
        self.ram_energy = {7: 1.5, 200: 2.5}
        self.sample_submit_sleep = 20


def init_pars():
    pars = Pars()
    # pars = {'build_path': '/home/dylan/git/Research/CooperFryeSampler/build/',
    #         'converter_path': '/home/dylan/git/Research/CooperFryeSampler/',
    #         'temp_path': '/home/dylan/Research/CFSample_Trees/temp',
    #         'out_dir': '/home/dylan/Research/CFSample_Trees/',
    #         'sampler_name': 'CooperFryeSampler',
    #         'ram_buffer': 2,  # Gb of RAM to keep free as buffer
    #         'ram_energy': {7: 1.5, 200: 2.5},
    #         'sample_submit_sleep': 20  # seconds to sleep between sample submission check
    #         }
    samplers = [{'energy': 7, 'rand_seed': 5, 'nevents': 100,
                 'hyper_input': f'{pars["converter_path"]}/input/input.AuAu.7.7.C0-5'}]

    return pars, samplers


def main():
    pars, samplers = init_pars()
    if not os.path.isdir(pars.temp_path):
        os.mkdir(pars.temp_path)
    run_samplers(pars, samplers)
    shutil.rmtree(pars.temp_path)
    # generate_input()  # Generate input file from sample plus passed parameters
    # run_sampler()
    # dat_to_root()
    gen_input_file()
    # dat_to_root()
    print('donzo')


def run_samplers(pars, samplers):
    sampler_index = 0
    while len(samplers) > 0:
        free_ram = psutil.virtual_memory().available / 1e9
        if pars.ram_energy[samplers[sampler_index]['energy']] < free_ram - pars.ram_buffer:
            run_sampler(samplers[sampler_index], pars)  # Run sampler at sampler_index if enough ram available

        sleep(pars.sample_submit_sleep)
        sampler_index += 1
        if sampler_index >= len(samplers):
            sampler_index = 0


def run_sampler(sampler, pars):
    file_name = f'{sampler["energy"]}GeV_rand{sampler["rand_seed"]}_nevent{sampler["nevent"]}'
    cd_build = f'cd {pars.build_path}'
    sampler_call = f'./{pars.sampler_name} {sampler["hyper_input"]} {pars.temp_path}' \
                   f''
    convert_call = f'{pars.out_dir}{sampler["energy"]}GeV/'
    full_command = f''


def run_sampler_test():
    path = '/home/dylan/git/Research/CooperFryeSampler/build/'
    path2 = '/home/dylan/git/Research/CooperFryeSampler/'
    in_path = '/home/dylan/Desktop/test.dat'
    out_path = '/home/dylan/Desktop/test.root'
    info = 'echo "test run"'
    chdir_sample = f'cd {path}'
    sampler_call = f'./CooperFryeSampler ../input/input.AuAu.7.7.C0-5 /home/dylan/Desktop/test.txt'
    chdir_convert = f'cd {path2}'
    convert_call = f'/home/dylan/Software/root/bin/root -b "CFSampleRootConvert.cpp(\\"{in_path}\\", \\"{out_path}\\")"'
    full_command = f'{info}; {chdir_sample}; {sampler_call}; sleep 15; {chdir_convert}; {convert_call}; sleep 15'
    # full_command = convert_call
    print(full_command)
    os.system(f'gnome-terminal -- bash -c \'echo "{info}"; cd {path}; {full_command}\'')


def dat_to_root():
    # while not os.path.isfile('/home/dylan/Desktop/test_finished.txt'):
    #     print('No flag, waiting...')
    #     sleep(5)
    # print('Flag found, convert: ')
    info = 'Convert dat to root: '
    path = '/home/dylan/git/Research/CooperFryeSampler/'
    in_path = '/home/dylan/Desktop/test.dat'
    out_path = '/home/dylan/Desktop/test.root'
    command = f'/home/dylan/Software/root/bin/root "CFSampleRootConvert.cpp(\\"{in_path}\\", \\"{out_path}\\")"'
    print(command)
    # command = f'root \'CFSampleRootConvert.cpp\("{in_path}", "{out_path}"\)\''
    os.system(f'gnome-terminal -- bash -c \'echo "{info}"; cd {path}; {command}; sleep 15\'')


def gen_input_file():
    sample_path = '/home/dylan/git/Research/CooperFryeSampler/input/input.AuAu.7.7.C0-5'
    temp_path = '/home/dylan/Desktop/test.input'
    parameters = {'nevents': 500, 'randomseed': 5, 'hypersurface_file': 'hydro/AuAu.7.7/C0-5/surface_eps_0.26.dat',
                  'output_file': '/home/dylan/Desktop/test.txt'}

    with open(sample_path, 'r') as file:
        sample_lines = file.readlines()

    out_lines = []
    for line in sample_lines:
        first_word = line.strip().split()
        if len(first_word) > 0:
            first_word = first_word[0]
        else:
            first_word = line
        if first_word in parameters:
            old_val = line.strip().split()[1]
            end_string = line[line.find(old_val) + len(old_val):]
            out_lines.append(f'{first_word} {parameters[first_word]}{end_string}')
        else:
            out_lines.append(line)

    with open(temp_path, 'w') as file:
        file.writelines(out_lines)

    print(sample_lines)

    print(out_lines)


if __name__ == '__main__':
    main()
