#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 05 12:55 AM 2020
Created in PyCharm
Created as QGP_Scripts/download_original.py

@author: Dylan Neff, dylan
"""


from subprocess import Popen
from time import sleep


def main():
    remote_path = 'dneff@rftpexp.rhic.bnl.gov:/gpfs01/star/pwg/dneff/scratch/'
    local_path = '/media/dylan/SSD_Storage/Research/'
    remote_tree_prefix = 'trees_old_ref'
    local_tree_prefix = 'Trees_Old_Ref'
    refs = [2, 3]
    energies = [11, 19, 27, 39, 62]
    for ref in refs:
        for energy in energies:
            remote = remote_path + remote_tree_prefix + f'{ref}/output/{energy}GeV'
            local = local_path + local_tree_prefix + f'{ref}/'
            print(remote, local)
            Popen(['gnome-terminal', '--', 'scp', '-r', remote, local])
            # sleep(5)
    print('donzo')


if __name__ == '__main__':
    main()
