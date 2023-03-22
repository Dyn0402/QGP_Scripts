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
    # data_set = 'AMPT_cent_sm'
    # data_set = 'CF'
    # data_set = 'CF_b342_lin'
    data_set = 'BES1_flow'
    data_sets = {'BES1': {'remote_path_suf': 'BES1/', 'remote_tree_pref': 'trees/output',
                          'local_path': 'C:/Users/Dylan/Research/', 'local_tree_pref': 'BES1_Trees'},
                 'BES1_flow': {'remote_path_suf': 'BES1/', 'remote_tree_pref': 'trees/output',
                               'local_path': 'F:/Research/', 'local_tree_pref': 'BES1_Trees'},
                 'AMPT_Run': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                              'local_path': 'C:/Users/Dylan/Research/',
                              'local_tree_pref': 'AMPT_Trees/min_bias/default'},
                 'AMPT_Run_mb': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                                 'local_path': 'F:/Research/', 'local_tree_pref': 'AMPT_Trees/min_bias/string_melting'},
                 'AMPT_Run_mcent_sm': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                                       'local_path': 'C:/Users/Dylan/Research/',
                                       'local_tree_pref': 'AMPT_Trees/most_central/string_melting'},
                 'AMPT_cent_def': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'most_central/default',
                                   'local_path': 'C:/Users/Dylan/Research/',
                                   'local_tree_pref': 'AMPT_Trees/most_central/default'},
                 'AMPT_cent_sm': {'remote_path_suf': 'AMPT/', 'local_path': 'F:/Research/',
                                  'remote_tree_pref': 'slim_most_central_new_coal/string_melting',
                                  'local_tree_pref': 'AMPT_Trees_New_Coalescence/slim_most_central/string_melting'},
                 'AMPT_cent_sm_lin': {'remote_path_suf': 'AMPT/', 'local_path': '/media/ucla/Research/',
                                      'remote_tree_pref': 'slim_most_central_new_coal/string_melting',
                                      'local_tree_pref': 'AMPT_Trees_New_Coalescence/slim_most_central/string_melting'},
                 'AMPT_mb_sm': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'min_bias/string_melting',
                                'local_path': 'F:/Research/',
                                'local_tree_pref': 'AMPT_Trees/min_bias/string_melting'},
                 'AMPT_gang': {'remote_path_suf': 'AMPT/', 'remote_tree_pref': 'dylan_run/output',
                               'local_path': 'F:/Research/', 'local_tree_pref': 'AMPT_Trees/gang'},
                 'CF': {'remote_path_suf': 'CooperFrye/', 'remote_tree_pref': 'CooperFrye_protons/output',
                        'local_path': 'F:/Research/', 'local_tree_pref': 'Cooper_Frye_EV_Trees'},
                 'CF_b342': {'remote_path_suf': 'CooperFrye/', 'remote_tree_pref': 'CooperFrye_b342_protons/output',
                             'local_path': 'F:/Research/', 'local_tree_pref': 'Cooper_Frye_EVb342_Trees'},
                 'CF_lin': {'remote_path_suf': 'CooperFrye/', 'remote_tree_pref': 'CooperFrye_protons/output',
                            'local_path': '/media/ucla/Research/', 'local_tree_pref': 'Cooper_Frye_EV_Trees'},
                 'CF_b342_lin': {'remote_path_suf': 'CooperFrye/', 'remote_tree_pref': 'CooperFrye_b342_protons/output',
                                 'local_path': '/media/ucla/Research/', 'local_tree_pref': 'Cooper_Frye_EVb342_Trees'},
                 }

    energies = [7, 11, 19, 27, 39, 62]  # , '2-7TeV_PbPb']
    # energies = [27, 62]
    bw_limit = None  # bandwidth limit per energy in Mbps or None
    size_tolerance = 0.001  # percentage tolerance between remote and local sizes, re-download if different
    file_delay = 0.1  # seconds to delay between file download calls

    remote_path = 'dneff@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/dneff/data/'

    remote_path += data_sets[data_set]['remote_path_suf']
    remote_tree_prefix = data_sets[data_set]['remote_tree_pref']
    local_tree_prefix = data_sets[data_set]['local_tree_pref']
    local_path = data_sets[data_set]['local_path']

    missing_files, bad_size_files = {}, {}
    all_missing = {}
    total_missing = 0
    for energy in energies:
        if type(energy) == int:
            energy = f'{energy}GeV'
        missing_files.update({energy: []})
        bad_size_files.update({energy: []})
        expected_files = get_expected_list(energy, remote_path, remote_tree_prefix)
        path = local_path + local_tree_prefix + f'/{energy}/'
        for file, file_size in expected_files.items():
            if os.path.exists(path + file):
                local_size = os.path.getsize(path + file)
                size_frac = (local_size - file_size) / file_size if file_size > 0 else 1 if local_size > 0 else 0
                if abs(size_frac) > size_tolerance:
                    bad_size_files[energy].append(path + file)
                    missing_files[energy].append(file)
            else:
                missing_files[energy].append(file)
        total_missing += len(missing_files[energy])
        all_missing[energy] = len(missing_files[energy]) == len(expected_files)
        print(f'{energy} missing {len(missing_files[energy])} of {len(expected_files)} files,'
              f' {len(bad_size_files[energy])} of these mismatched size')

    num_bad_files = sum([len(bad_files_energy) for bad_files_energy in bad_size_files.values()])
    if num_bad_files > 0:
        res = input(f'Delete {num_bad_files} bad size files?\n')
        if res.strip().lower()[0] in ['yes', 'y']:
            for bad_energy_files in bad_size_files.values():
                for bad_file in bad_energy_files:
                    # os.system(f'DEL {bad_file}')  # If files say they're busy, this should skip recycle
                    os.remove(bad_file)  # Seemed to have issues with files being busy. If so use above line
    if total_missing > 0:
        res = input(f'\nDownload {total_missing} missing files? Enter yes to download all; energy name to download'
                    f' a single energy; energy name,number of files to download only first n files for energy;'
                    f' or anything else to quit: \n')
        energy_list = []
        if res.strip().lower() in ['yes', 'y']:
            energy_list = missing_files.keys()
        else:
            energy_in = res.strip()
            if energy_in in missing_files:
                energy_list.append(energy_in)
            else:
                try:
                    in_list = res.strip().split(',')
                    energy_in = in_list[0].strip()
                    if energy_in in missing_files:
                        energy_list.append(energy_in)
                        num_files = int(in_list[1].strip())
                        missing_files[energy_in] = missing_files[energy_in][:num_files]
                except ValueError:
                    pass
        for energy in energy_list:
            local = local_path + local_tree_prefix + f'/{energy}/'
            if len(missing_files[energy]) > 0:
                if all_missing[energy]:
                    start_download_all(energy, remote_path, remote_tree_prefix, local, bw_limit)
                    sleep(file_delay)
                else:
                    start_download_sftp(missing_files[energy], energy, remote_path, remote_tree_prefix, local, bw_limit)
                    # for file in missing_files[energy]:
                    #     start_download(file, energy, remote_path, remote_tree_prefix, local)
                    #     sleep(file_delay)

    else:
        print('All files downloaded!')


def get_expected_list(energy, remote_path, remote_tree_prefix):
    cmd = f'echo ls -l|sftp {remote_path.split(":")[0]}:{remote_path.split(":")[1]}{remote_tree_prefix}/{energy}'
    stdout, stderr = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()

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


def start_download(file, energy, remote_path, remote_tree_prefix, local, bw_limit=None):
    remote = remote_path + remote_tree_prefix + f'/{energy}/{file}'
    bw_limit_str = '' if bw_limit is None else f'-l {int(bw_limit * 1000)} '
    command = 'sftp ' + bw_limit_str + remote + ' ' + local + file
    info = f'{energy}, {file} files:'
    print(f'{info} {command}')
    os.system(f'start cmd /c {command}')


def start_download_all(energy, remote_path, remote_tree_prefix, local, bw_limit=None):
    remote = remote_path + remote_tree_prefix + f'/{energy}/*'
    bw_limit_str = '' if bw_limit is None else f'-l {int(bw_limit * 1000)} '
    command = 'sftp ' + bw_limit_str + remote + ' ' + local
    info = f'{energy} all files:'
    print(f'{info} {command}')
    os.system(f'start cmd /c {command}')


def start_download_sftp(files, energy, remote_path, remote_tree_prefix, local, bw_limit=None):
    remote_host, remote_path = remote_path.split(':')
    bat_file_name = f'{energy}_{remote_tree_prefix.split("/")[0]}_sftp_file.bat'

    sftp_file_name = f'{energy}_{remote_tree_prefix.split("/")[0]}_sftp_file.txt'
    sftp_gets = [f'get {remote_path}{remote_tree_prefix}/{energy}/{file} {local}{file}\n' for file in files]
    sftp_gets.insert(0, 'progress\n')  # Show progress of downloads
    with open(sftp_file_name, 'w') as temp_txt:
        temp_txt.writelines(sftp_gets)

    bw_limit_str = '' if bw_limit is None else f'-l {int(bw_limit * 1000)}'
    command = f'sftp -b {sftp_file_name} {bw_limit_str} {remote_host}'

    with open(bat_file_name, 'w') as temp_bat:  # Need bat file to delete itself and sftp batch
        temp_bat.write(f'{command}\n')  # Run sftp batch file
        temp_bat.write(f'del {sftp_file_name}\n')  # Delete sftp batch file
        temp_bat.write(f'start /b "" cmd /c del {bat_file_name}&exit /b')  # Have bat file delete itself

    info = f'{energy}, {len(files)} files:'
    print(f'{info} {command}')
    print('  '.join(sftp_gets))

    os.system(f'start cmd /c {bat_file_name}')  # Run batch file in new terminal


if __name__ == '__main__':
    main()
