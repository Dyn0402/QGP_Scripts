#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 06 8:06 PM 2020
Created in PyCharm
Created as QGP_Scripts/ref3_compare.py

@author: Dylan Neff, dylan
"""


def main():
    git_ref3_def = get_git_ref3_def()
    xiao_ref3_def = get_xiao_ref3_def()
    compare_ref3_def(git_ref3_def, xiao_ref3_def)
    print('donzo')


def get_git_ref3_def():
    path = '/home/dylan/git/Research/QGP_Fluctuations/Tree_Reader/StRefMultCorr2/Param.h'
    ref3_start_line = 596
    ref3_set_step = 604 - 596
    ref3_end_line = 770

    ref3_def = {'runIDs': [], 'z_vertex': [], 'bin_edges': [], 'Norm': None, 'z_pars': [], 'w_pars': [], 'l_pars': []}
    ref3_defs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        line_num = ref3_start_line
        while line_num < ref3_end_line:
            new_def = Ref3Def()
            line = lines[line_num].strip().strip(',').strip('"').split(':')
            ids = line[2].split(',')
            new_def.runID['start'] = int(ids[0])
            new_def.runID['end'] = int(ids[1])
            zs = line[-1].split(',')
            new_def.z_vertex['low'] = float(zs[0])
            new_def.z_vertex['high'] = float(zs[1])
            bins = lines[line_num + 1].strip().strip(',').strip('"').split(',')
            for bin_edge in bins:
                new_def.bin_edges.append(int(bin_edge))
            new_def.norm = int(lines[line_num+2].strip().strip(',').strip('"'))
            zpars = lines[line_num + 3].strip().strip(',').strip('"').split(',')
            for zpar in zpars:
                new_def.z_pars.append(float(zpar))
            wpars = lines[line_num + 4].strip().strip(',').strip('"').split(',')
            for wpar in wpars:
                new_def.w_pars.append(float(wpar))
            lpars = lines[line_num + 5].strip().strip(',').strip('"').split(',')
            for lpar in lpars:
                new_def.l_pars.append(float(lpar))

            ref3_defs.append(new_def)
            line_num += ref3_set_step

    return ref3_defs


def get_xiao_ref3_def():
    path = '/home/dylan/Research/Results/4-7-20/Centrality_def.txt'

    ref3_defs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            line = line.strip().split()
            new_ref3 = Ref3Def()
            new_ref3.runID['start'] = int(line[0])
            new_ref3.runID['end'] = int(line[1])
            new_ref3.z_vertex['low'] = float(line[2])
            new_ref3.z_vertex['high'] = float(line[3])
            for i in range(16):
                new_ref3.bin_edges.append(int(line[4+i]))
            new_ref3.norm = int(line[20])
            for i in range(8):
                new_ref3.z_pars.append(float(line[21+i]))
            for i in range(6):
                new_ref3.w_pars.append(float(line[29+i]))
            new_ref3.w_pars.append(0.)
            new_ref3.w_pars.append(0.)
            for i in range(2):
                new_ref3.l_pars.append(float(line[35+i]))

            ref3_defs.append(new_ref3)

    return ref3_defs


def compare_ref3_def(git, xiao):
    print(f'Number of git: {len(git)}  |  Number of xiao: {len(xiao)}')
    for i in range(len(xiao)):
        if git[i] != xiao[i]:
            print(f'i: {i}\n')
            print('git:')
            print(str(git[i]))
            print('xiao: ')
            print(str(xiao[i]))
            input()


class Ref3Def:
    #  Attributes
    def __init__(self):
        self.runID = {'start': None, 'end': None}
        self.z_vertex = {'low': None, 'high': None}
        self.bin_edges = []
        self.norm = None
        self.z_pars = []
        self.w_pars = []
        self.l_pars = []

    def __eq__(self, other):
        if not isinstance(other, Ref3Def):
            return False
        if self.runID != other.runID:
            return False
        if self.z_vertex != other.z_vertex:
            return False
        if self.bin_edges != other.bin_edges:
            return False
        if self.norm != other.norm:
            return False
        if self.z_pars != other.z_pars:
            return False
        if self.w_pars != other.w_pars:
            return False
        if self.l_pars != other.l_pars:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        string = ''
        string += f'ID: {self.runID["start"]}, {self.runID["end"]}\n'
        string += f'z_vertex: {self.z_vertex["low"]}, {self.z_vertex["high"]}\n'
        string += f'norm: {self.norm}\n'
        string += 'bin_edges: '
        for bin_edge in self.bin_edges:
            string += f'{bin_edge}, '
        string = string[:-2] + '\nz_pars: '
        for z_par in self.z_pars:
            string += f'{z_par}, '
        string = string[:-2] + '\nw_pars: '
        for w_par in self.w_pars:
            string += f'{w_par}, '
        string = string[:-2] + '\nl_pars: '
        for l_par in self.l_pars:
            string += f'{l_par}, '
        string = string[:-2] + '\n'

        return string


if __name__ == '__main__':
    main()
