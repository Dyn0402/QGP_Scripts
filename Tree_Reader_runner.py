#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 24 3:45 PM 2023
Created in PyCharm
Created as QGP_Scripts/Tree_Reader_runner.py

@author: Dylan Neff, dylan
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


def main():
    exe_path = '/home/dylan/git/Research/QGP_Fluctuations/Tree_Reader/Release/Tree_Reader'
    root_lib_path = '/home/dylan/Software/root/lib'

    energies = [7]

    jobs = {
        'BES1_def_nofileshuffle_test': {
            'default_test': {
                'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_nomix_': [0, 0]
            }
        },
        # 'BES1_def_nofileshuffle': {
        #     'default': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_dca08': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca08_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_dca12': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca12_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nsprx09': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx09_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nsprx11': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx11_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nhfit15': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit15_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nhfit25': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit25_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_m2r8': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r8_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_m2r4': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r4_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_dca05': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca05_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_dca15': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca15_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nsprx075': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx075_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_nsprx125': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx125_m2r6_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_m2r2': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r2_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
        # 'BES1_sys_nofileshuffle_m2r10': {
        #     'default_sys': {
        #         'rapid05_resample_norotate_seed_dca1_nsprx1_m2r10_m2s0_nhfit20_epbins1_calcv2_': [0, 0]
        #     }
        # },
    }

    for energy in energies:
        for job_type, job_type_dict in jobs.items():
            for set_group_name, set_group_dict in job_type_dict.items():
                for set_name, job_nums in set_group_dict.items():
                    command = f'export LD_LIBRARY_PATH={root_lib_path}; stty rows 40 columns 150;' \
                              f'{exe_path} {energy} {job_type} {set_group_name} {set_name} {job_nums[0]} {job_nums[1]}'
                    subprocess.Popen(['gnome-terminal', '--geometry=120x20', '--', 'bash', '-c', command])
                    # # subprocess.Popen(['bash', '-i', '-c', 'source ~/.bashrc'])
                    # os.system('gnome-terminal -- bash -c "{}"'.format(command))


if __name__ == '__main__':
    main()
