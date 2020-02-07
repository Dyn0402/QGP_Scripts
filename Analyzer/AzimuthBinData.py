#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 06 7:27 PM 2020
Created in PyCharm
Created as QGP_Scripts/AzimuthBinData.py

@author: Dylan Neff, dylan
"""


class AzimuthBinData:
    def __init__(self, div='', path=''):
        self.path = path
        self.div = div
        self.data = {}
        if path != '':
            self.read_data()

    def read_data(self):
        with open(self.path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split('\t')
                total_protons = int(line[0])
                bin_protons_string = line[1].split(' ')
                if len(bin_protons_string) > 0:
                    self.data[total_protons] = []
                i = 0
                for entry in bin_protons_string:
                    entry = entry.split(':')
                    while int(entry[0]) > i:
                        self.data[total_protons].append(0)
                        i += 1
                    self.data[total_protons].append(int(entry[1]))
                    i += 1

