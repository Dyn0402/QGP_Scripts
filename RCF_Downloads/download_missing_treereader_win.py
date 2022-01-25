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
    data_set = 'Data_Sim'
    data_sets = {'Data_Sim': {'remote_path_suf': 'tree_reader_data/Data_Sim/',
                              'local_path': 'D:/Research/Data_Sim/'},
                 'Data_Sim_Mix': {'remote_path_suf': 'tree_reader_data/Data_Sim_Mix/',
                                  'local_path': 'D:/Research/Data_Sim_Mix/'}
                 }

    size_tolerance = 0.001  # percentage tolerance between remote and local sizes, re-download if different
    file_delay = 0.1  # seconds to delay between file download calls

    remote_path = 'dneff@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/dneff/'
    rftp_path = 'dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/'

    remote_path += data_sets[data_set]['remote_path_suf']
    rftp_path += data_sets[data_set]['remote_path_suf']
    local_path = data_sets[data_set]['local_path']

    missing_files = []
    expected_files = get_expected_list(rftp_path)
    bad_size_files = 0
    for file, file_size in expected_files.items():
        print(local_path + file, os.path.exists(local_path + file))
        if os.path.exists(local_path + file):
            size_frac = (os.path.getsize(local_path + file) - file_size) / file_size
            # if abs(size_frac) > size_tolerance:
            #     bad_size_files += 1
            #     missing_files.append(file)
        else:
            missing_files.append(file)
    total_missing = len(missing_files)
    print(f'Missing {len(missing_files)} of {len(expected_files)} files')
    # f', {bad_size_files} of these mismatched size')

    if total_missing > 0:
        res = input(f'\nDownload {total_missing} missing files? Enter yes to download all; '
                    f' number of files to download only first n files;'
                    f' or anything else to quit: \n')
        if res.strip().lower() in ['yes', 'y']:
            pass
        else:
            try:
                num_files = int(res.strip())
                missing_files = missing_files[:num_files]
            except ValueError:
                print('Downloading nothing. Exiting.')
                return
        for file in missing_files:
            start_download(file, remote_path, local_path)
            sleep(file_delay)

    else:
        print('All files downloaded!')


def get_expected_list(remote_path):
    stdout, stderr = Popen(['ssh', f'{remote_path.split(":")[0]}', 'ls -l',
                            f'{remote_path.split(":")[1]}'],
                           stdout=PIPE, stderr=PIPE).communicate()
    print(stderr)
    files_str = stdout.decode('UTF-8').split('\n')
    files_dict = {}
    for file in files_str:
        file = file.split()
        if len(file) == 9:
            files_dict.update({file[-1]: int(file[4])})

    return files_dict


def start_download(file, remote_path, local):
    remote = remote_path + f'{file}'
    command = f'sftp -r {remote} {local}{file}'
    # command = 'sftp -r ' + remote + ' ' + local
    info = f'{file} files:'
    print(f'{info} {command}')
    os.system(f'start cmd /c {command}')


if __name__ == '__main__':
    main()
