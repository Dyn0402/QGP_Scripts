#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 21 10:35 PM 2022
Created in PyCharm
Created as QGP_Scripts/eps_pdf_convert

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import os


def main():
    fig_path = '/home/dylan/Downloads/figures/'
    for sub_dir in os.listdir(fig_path):
        for file_path in os.listdir(f'{fig_path}{sub_dir}'):
            if '.eps' in file_path:
                print(f'ps2pdf {fig_path}{sub_dir}/{file_path} {fig_path}{sub_dir}/{file_path.split(".")[0]}.pdf')
                os.system(f'ps2pdf {fig_path}{sub_dir}/{file_path} {fig_path}{sub_dir}/{file_path.split(".")[0]}.pdf')
            os.remove(f'{fig_path}{sub_dir}/{file_path}')
    print('donzo')


if __name__ == '__main__':
    main()
