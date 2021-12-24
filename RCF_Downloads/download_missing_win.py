#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 05 6:35 PM 2021
Created in PyCharm
Created as QGP_Scripts/download_missing_win.py

@author: Dylan Neff, Dyn04
"""

import os
from subprocess import Popen, PIPE
from time import sleep


def main():
    download()
    print('donzo')


def download():
    data_set = 'AMPT_cent_sm'
    data_sets = {'BES1': {'remote_path_suf': 'BES1/', 'remote_tree_pref': 'trees/output',
                          'local_path': 'C:/Users/Dylan/Research/', 'local_tree_pref': 'BES1_Trees'},
                 'AMPT_Run': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                              'local_path': 'C:/Users/Dylan/Research/', 'local_tree_pref': 'AMPT_Trees/min_bias/default'},
                 'AMPT_Run_mcent_sm': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                              'local_path': 'C:/Users/Dylan/Research/', 'local_tree_pref': 'AMPT_Trees/most_central/string_melting'},
                 'AMPT_cent_def': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'most_central/default',
                                   'local_path': 'C:/Users/Dylan/Research/', 'local_tree_pref': 'AMPT_Trees/most_central/default'},
                 'AMPT_cent_sm': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'most_central/string_melting',
                                   'local_path': 'C:/Users/Dylan/Research/',
                                   'local_tree_pref': 'AMPT_Trees/most_central/string_melting'}}

    energies = [11]  # [7, 11, 15, 19, 27, 39, 62]
    bwlimit = 10  # bandwidth limit per energy in MBPS or None
    size_tolerance = 0.001  # percentage tolerance between remote and local sizes, re-download if different
    file_delay = 1  # seconds to delay between file download calls

    remote_path = 'dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/data/'

    remote_path += data_sets[data_set]['remote_path_suf']
    remote_tree_prefix = data_sets[data_set]['remote_tree_pref']
    local_tree_prefix = data_sets[data_set]['local_tree_pref']
    local_path = data_sets[data_set]['local_path']

    missing_files = {}
    total_missing = 0
    for energy in energies:
        missing_files.update({energy: []})
        expected_files = get_expected_list(energy, remote_path, remote_tree_prefix)
        path = local_path + local_tree_prefix + f'/{energy}GeV/'
        bad_size_files = 0
        for file, file_size in expected_files.items():
            if os.path.exists(path + file):
                size_frac = (os.path.getsize(path + file) - file_size) / file_size
                if abs(size_frac) > size_tolerance:
                    bad_size_files += 1
                    missing_files[energy].append(file)
            else:
                missing_files[energy].append(file)
        total_missing += len(missing_files[energy])
        print(f'{energy}GeV missing {len(missing_files[energy])} of {len(expected_files)} files,'
              f' {bad_size_files} of these mismatched size')

    if total_missing > 0:
        res = input(f'\nDownload {total_missing} missing files? Enter yes to download all; energy number to download'
                    f' a single energy; energy number,number of files to download only first n files for energy;'
                    f' or anything else to quit: \n')
        energy_list = []
        if res.strip().lower() in ['yes', 'y']:
            energy_list = missing_files.keys()
        else:
            try:
                energy_in = int(res.strip().lower())
                if energy_in in missing_files:
                    energy_list.append(energy_in)
            except ValueError:
                try:
                    in_list = res.strip().lower().split(',')
                    energy_in = int(in_list[0].strip())
                    num_files = int(in_list[1].strip())
                    if energy_in in missing_files:
                        energy_list.append(energy_in)
                        missing_files[energy_in] = missing_files[energy_in][:num_files]
                except ValueError:
                    pass
        for energy in energy_list:
            local = local_path + local_tree_prefix + f'/{energy}GeV/'
            if len(missing_files[energy]) > 0:
                # files = r'\{'
                for file in missing_files[energy]:
                    # files += file + ','
                    start_download(file, energy, remote_path, remote_tree_prefix, local, bwlimit)
                    sleep(1)

    else:
        print('All files downloaded!')


def get_expected_list(energy, remote_path, remote_tree_prefix):
    stdout, stderr = Popen(['ssh', f'{remote_path.split(":")[0]}', 'ls -l',
                            f'{remote_path.split(":")[1]}{remote_tree_prefix}/{energy}GeV'],
                           stdout=PIPE, stderr=PIPE).communicate()
    files_str = stdout.decode('UTF-8').split('\n')
    files_dict = {}
    for file in files_str:
        file = file.split()
        if len(file) == 9:
            files_dict.update({file[-1]: int(file[4])})

    return files_dict


# def start_download(files, energy, remote_path, remote_tree_prefix, local, bwlimit=None):
#     files = files[:-1] + r'\}'
#     remote = remote_path + remote_tree_prefix + f'/{energy}GeV/{files}'
#     if bwlimit:
#         command = f'rsync -r -v --bwlimit={bwlimit * 1000} --progress -e ssh ' + remote + ' ' + local
#     else:
#         command = 'rsync -r -v --progress -e ssh ' + remote + ' ' + local
#     info = f'{energy}GeV, {files.count(",") + 1} files:'
#     print(f'{info} {command}')
#     os.system(f'gnome-terminal -- /bin/sh -c \'echo "{info}"; {command}\'')


def start_download(file, energy, remote_path, remote_tree_prefix, local, bwlimit=None):
    # files = files[:-1] + r'\}'
    remote = remote_path + remote_tree_prefix + f'/{energy}GeV/{file}'
    command = 'scp ' + remote + ' ' + local
    info = f'{energy}GeV, {file} files:'
    print(f'{info} {command}')
    os.system(f'start cmd /c {command}')


if __name__ == '__main__':
    main()