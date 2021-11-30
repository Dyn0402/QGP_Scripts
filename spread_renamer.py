#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 23 11:28 AM 2021
Created in PyCharm
Created as QGP_Scripts/spread_renamer

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil


def main():
    amps = ['_amp01_', '_amp02_', '_amp03_', '_amp04_']
    spread = '_spread01'
    base_path = '/home/dylan/Research/Data_Sim/'
    for path in os.listdir(base_path):
        dir_path = base_path + path
        good = True
        # for amp in amps:
        #     if amp in path and spread in path:
        #         good = True
        #         break
        if good:
            if os.path.isdir(dir_path):
                for path2 in os.listdir(dir_path):
                    dir_path2 = dir_path + '/' + path2
                    if os.path.isdir(dir_path2):
                        dir_path2_new = dir_path + '/' + path2.replace('spread', 'spread0')
                        print(f'{dir_path2} --> \n{dir_path2_new}')
                        shutil.move(dir_path2, dir_path2_new)
            dir_path_new = base_path + path.replace('spread', 'spread0')
            print(f'{dir_path} --> \n{dir_path_new}')
            shutil.move(dir_path, dir_path_new)
    print('donzo')


if __name__ == '__main__':
    main()
