#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 09 4:25 PM 2021
Created in PyCharm
Created as QGP_Scripts/SysReader

@author: Dylan Neff, dylan
"""

import pandas as pd


class SysReader:
    def __init__(self, path):
        self.sys_path = path
        self.cols = ['set', 'data_type', 'energy', 'div', 'cent', 'stat', 'val', 'stat_err', 'sys_err']
        vals = self.read()
        self.values = pd.DataFrame(vals, columns=self.cols)

    def read(self):
        vals = []
        with open(self.sys_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split(',')
                if len(line) != 8:
                    print('Bad read: ' + line)
                else:
                    for i in [3, 4, 5]:  # energy, div, cent
                        line[i] = int(line[i])
                    if line[0] == 'medians':  # Either medians or systematics
                        med_val, stat_err = line[-1].split('±')
                        vals.append(line[1:-1] + [float(med_val), float(stat_err)])
                    if line[0] == 'systematics':
                        if vals[-1][:-2] != line[1:-1]:
                            print('Systematic line doesn\'t match previous medians: ' + line)
                        else:
                            vals[-1].append(float(line[-1]))

        return vals
