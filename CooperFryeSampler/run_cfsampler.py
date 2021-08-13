#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 12 7:21 PM 2021
Created in PyCharm
Created as QGP_Scripts/run_cfsampler

@author: Dylan Neff, dylan
"""

import numpy as np
import sys
import subprocess as sp
from random import randint


def main():
    if len(sys.argv) == 3:
        energy = sys.argv[1]
        hypersurface = sys.argv[2]
        nevents = 10000
    elif len(sys.argv) == 4:
        energy = sys.argv[1]
        hypersurface = sys.argv[2]
        nevents = sys.argv[3]
    randomseed = randint(1, 2147483647)
    print('donzo')


def run():
    pass


if __name__ == '__main__':
    main()
