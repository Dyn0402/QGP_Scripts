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
from random import randint


class Pars:
    def __init__(self):
        self.build_path = '/home/dylan/git/Research/CooperFryeSampler/build/'
        self.converter_path = '/home/dylan/git/Research/CooperFryeSampler/'
        self.temp_path = '/home/dylan/Desktop/CooperFryeSampler_temp/'  # '/home/dylan/Research/CFSample_Trees/temp/'
        self.out_dir = '/home/dylan/Research/CFSample_Trees/'
        self.root_path = '/home/dylan/Software/root/bin/root'
        self.in_sample_path = '/home/dylan/git/Research/CooperFryeSampler/input/input.AuAu.7.7.C0-5'
        self.sampler_name = 'CooperFryeSampler'
        self.converter_name = 'CFSampleRootConvert.cpp'
        self.ram_buffer = 1
        self.ram_energy = {7: 6.4, 19: 6.9, 27: 10.5, 39: 11.6, 62: 12.4}  # ram estimate for each energy in GB
        self.time_energy = {7: 17, 19: 32, 27: 36, 39: 43, 62: 52}  # time estimate per event for each energy in ms
        self.sample_submit_sleep = 20
        self.random_seed_range = (1, 2147483647)


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
    meta_samplers = [{'energy': 62, 'nevents': 1001, 'target_events': 2000,
                      'hypersurface_file': f'{pars.build_path}/hydro/AuAu.62.4/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 27, 'randomseed': 2147483649, 'nevents': 1000,
                     #  'hypersurface_file': f'{pars.build_path}hydro/AuAu.27/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 19, 'randomseed': 5, 'nevents': 1000,
                     #  'hypersurface_file': f'{pars.build_path}hydro/AuAu.19.6/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 39, 'randomseed': 5, 'nevents': 1000,
                     #  'hypersurface_file': f'{pars.build_path}hydro/AuAu.39/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 62, 'randomseed': 5, 'nevents': 1000,
                     #  'hypersurface_file': f'{pars.build_path}hydro/AuAu.62.4/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 27, 'randomseed': 5, 'nevents': 1000,
                     #  'hypersurface_file': f'{pars.build_path}hydro/AuAu.27/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 7, 'randomseed': 83, 'nevents': 100000,
                     #  'hypersurface_file': f'{pars.build_path}/hydro/AuAu.7.7/C0-5/surface_eps_0.26.dat'},
                     # {'energy': 7, 'randomseed': 171, 'nevents': 100000,
                     #  'hypersurface_file': f'{pars.build_path}/hydro/AuAu.7.7/C0-5/surface_eps_0.26.dat'},
                     ]

    return pars, meta_samplers


def main():
    pars, meta_samplers = init_pars()
    if os.path.isdir(pars.temp_path):
        shutil.rmtree(pars.temp_path)
    os.mkdir(pars.temp_path)
    samplers = gen_samplers(meta_samplers, pars)
    run_samplers(pars, samplers)
    # shutil.rmtree(pars.temp_path)
    # generate_input()  # Generate input file from sample plus passed parameters
    # run_sampler()
    # dat_to_root()
    # gen_input_file()
    # dat_to_root()
    print('donzo')


def gen_samplers(meta_samplers, pars):
    samplers = []
    for meta_sampler in meta_samplers:
        for i in range(1 + int(meta_sampler['target_events'] / meta_sampler['nevents'])):
            sampler = {key: val for key, val in meta_sampler.items()}
            sampler.update({'randomseed': randint(pars.random_seed_range[0], pars.random_seed_range[1])})
            samplers.append(sampler)

    return samplers


def run_samplers(pars, samplers):
    # Pair expected ram with each sampler and sort from most RAM usage to least so high RAM processes run first
    # print([pars.ram_energy[x['energy']] for x in samplers])
    # print(samplers)
    # print([x for x in zip([pars.ram_energy[x['energy']] for x in samplers], samplers)])
    samplers = sorted(zip([pars.ram_energy[x['energy']] for x in samplers], samplers), key=lambda y: y[0], reverse=True)

    while len(samplers) > 0:
        free_ram = psutil.virtual_memory().available / 1e9  # See how much RAM available on system
        if free_ram >= min(pars.ram_energy.values()) + pars.ram_buffer:
            for sampler_index in range(len(samplers)):  # Search for largest RAM sampler that will fit and run it
                if pars.ram_energy[samplers[sampler_index][1]['energy']] < free_ram - pars.ram_buffer:
                    run_sampler(samplers[sampler_index][1], pars)  # Run sampler at sampler_index if ram available
                    samplers.pop(sampler_index)  # Remove sampler from list once finished
                    break
        sleep(pars.sample_submit_sleep)  # Wait long enough for submitted process to allocate RAM before checking again

        # if pars.ram_energy[samplers[sampler_index]['energy']] < free_ram - pars.ram_buffer:
        #     run_sampler(samplers[sampler_index], pars)  # Run sampler at sampler_index if enough ram available
        #     samplers.pop(sampler_index)  # Remove sampler from list once finished
        #
        # sleep(pars.sample_submit_sleep)
        # sampler_index += 1
        # if sampler_index >= len(samplers):
        #     sampler_index = 0


def run_sampler(sampler, pars):
    file_name = f'{sampler["energy"]}GeV_rand{sampler["randomseed"]}_nevents{sampler["nevents"]}'
    info = f'echo "Running {file_name}"' \
           f'Estimated time {pars.time_energy[sampler["energy"]] * sampler["nevents"] / 1000 / 60} mins'
    temp_dat_path = f'{pars.temp_path}{file_name}.dat'

    sampler.update({'output_file': temp_dat_path, 'temp_input_path': pars.temp_path + file_name + '.input'})
    gen_input_file(sampler, pars)

    out_root_path = f'{pars.out_dir}{sampler["energy"]}GeV/{file_name}.root'
    cd_build = f'cd {pars.build_path}'
    sampler_call = f'./{pars.sampler_name} {sampler["temp_input_path"]} {temp_dat_path}'
    cd_convert = f'cd {pars.converter_path}'
    convert_call = f'{pars.root_path} -b -q "{pars.converter_name}(\\"{temp_dat_path}\\", \\"{out_root_path}\\")"'
    rm_dat = f'rm {temp_dat_path}'
    full_command = f'{info}; {cd_build}; {sampler_call}; sleep 15; {cd_convert}; {convert_call}; {rm_dat}; read x;'
    print(full_command)
    os.system(f'gnome-terminal -- bash -c \'{full_command}\'')


def gen_input_file(sampler, pars):
    with open(pars.in_sample_path, 'r') as file:
        sample_lines = file.readlines()

    out_lines = []
    for line in sample_lines:
        first_word = line.strip().split()
        if len(first_word) > 0:
            first_word = first_word[0]
        else:
            first_word = line
        if first_word in sampler:
            old_val = line.strip().split()[1]
            end_string = line[line.find(old_val) + len(old_val):]
            out_lines.append(f'{first_word} {sampler[first_word]}{end_string}')
        else:
            out_lines.append(line)

    with open(sampler['temp_input_path'], 'w') as file:
        file.writelines(out_lines)

    print(sample_lines)

    print(out_lines)


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


def gen_input_file_test():
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
