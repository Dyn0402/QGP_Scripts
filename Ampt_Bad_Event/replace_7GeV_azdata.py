#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 5:24 PM 2021
Created in PyCharm
Created as QGP_Scripts/replace_7GeV_azdata

@author: Dylan Neff, dylan
"""


import os
import shutil


def main():
    action = True
    set_folder = '7GeV'
    origin_base_path = '/media/dylan/DYLAN_NEFF/Transfer/Research/'
    out_base_path = '/home/dylan/Research/'
    dirs = ['Data_Ampt', 'Data_Ampt_Mix']
    for d in dirs:
        d_path = origin_base_path + d + '/'
        for set_group in os.listdir(d_path):
            sg_path = d_path + set_group + '/'
            sg_out_path = f'{out_base_path}{d}/{set_group}/'
            if not os.path.exists(sg_out_path):
                print(f'Make dir: {sg_out_path}')
                if action:
                    os.mkdir(sg_out_path)
            for set_name in os.listdir(sg_path):
                in_path = sg_path + set_name + '/' + set_folder
                out_path = f'{out_base_path}{d}/{set_group}/{set_name}/{set_folder}'
                print(f'{in_path} --> {out_path}')
                if action:
                    if os.path.exists(out_path):
                        shutil.rmtree(out_path)
                    shutil.copytree(in_path, out_path)
    print('donzo')


if __name__ == '__main__':
    main()
