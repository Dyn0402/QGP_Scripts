#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 14 1:46 AM 2023
Created in PyCharm
Created as QGP_Scripts/v2_neg_res_fixer.py

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
from Binom_Slices.analyze_binom_slices import read_flow_values
import uproot


def main():
    def_path = 'F:/Research/Data/default_sys/'

    for dir_name in os.listdir(def_path):
        set_path = f'{def_path}{dir_name}/'
        v2s = read_flow_values(set_path)
        for energy in v2s:
            for cent in v2s[energy]:
                if np.isnan(v2s[energy][cent].val):
                    qa_path = f'{set_path}{energy}GeV/QA_{energy}GeV.root'
                    with uproot.open(qa_path) as root_file:
                        res_prof = root_file[f'resolution_{dir_name}_{energy}_{cent}']
                        print(dir_name)
                        print(energy)
                        print(cent)
                        print(res_prof, '\n')

    print('donzo')


if __name__ == '__main__':
    main()
