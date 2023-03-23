#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 12 7:21 PM 2021
Created in PyCharm
Created as QGP_Scripts/run_cfsampler

@author: Dylan Neff, dylan
"""

import sys
import os
from random import randint


def main():
    # set_run_from_cmd()
    just_get_input_file()
    print('donzo')


def set_run_from_cmd():
    if len(sys.argv) >= 3:
        hypersurface_name = sys.argv[1]
        nevents = sys.argv[2]
        hypersurface_energy = 'AuAu.7.7'  # Default to 7GeV
        input_sample_name = 'input.AuAu.7.7.C0-5.EVHRG'  # Default to original 7GeV template
        save_all_particles = 0  # Only save protons
    if len(sys.argv) >= 4:
        hypersurface_energy = sys.argv[3]
    if len(sys.argv) >= 5:
        input_sample_name = sys.argv[4]
    if len(sys.argv) >= 6:
        save_all_particles = int(sys.argv[5])
    randomseed = randint(1, 2147483647)
    sampler = {'energy': hypersurface_energy + '_', 'randomseed': randomseed, 'nevents': nevents,
               'hypersurface_file': hypersurface_name, 'input_sample_name': input_sample_name,
               'cfsampler_name': 'CooperFryeSampler', 'converter_name': 'CFSampleRootConvert.cpp',
               'root_path': '/star/u/dneff/Software/root/root-6.14.06/obj/bin/root',
               'save_all_particles': save_all_particles}
    sampler.update({'output_file': f'{sampler["energy"]}nevents_{sampler["nevents"]}_rand_{sampler["randomseed"]}'})
    run(sampler)


def run(sampler):
    input_path = gen_input_file(sampler)
    dat_name = sampler['output_file'] + '.dat'
    output_root = sampler['output_file'] + '.root'
    os.system(f'./{sampler["cfsampler_name"]} {input_path} {dat_name}')
    os.system(f'{sampler["root_path"]} -b -q "{sampler["converter_name"]}(\\"{dat_name}\\", \\"{output_root}\\", '
              f'\\"{sampler["save_all_particles"]}\\")"')


def just_get_input_file():
    sampler = {'energy': 'AuAu.7.7' + '_', 'randomseed': 1, 'nevents': '50000',
               'hypersurface_file': 'surface_eps_0.26.dat', 'input_sample_name': 'input.AuAu.7.7.C0-5.EVHRG',
               'cfsampler_name': 'CooperFryeSampler', 'converter_name': 'CFSampleRootConvert.cpp',
               'root_path': '/star/u/dneff/Software/root/root-6.14.06/obj/bin/root'}
    sampler.update({'output_file': f'{sampler["energy"]}nevents_{sampler["nevents"]}_rand_{sampler["randomseed"]}'})
    gen_input_file(sampler)


def gen_input_file(sampler):
    with open(sampler['input_sample_name'], 'r') as file:
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

    temp_input_path = 'input.' + sampler['output_file']
    with open(temp_input_path, 'w') as file:
        file.writelines(out_lines)

    return temp_input_path


if __name__ == '__main__':
    main()
