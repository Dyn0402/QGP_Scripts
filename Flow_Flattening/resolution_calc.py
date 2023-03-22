#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 21 7:16 PM 2023
Created in PyCharm
Created as QGP_Scripts/resolution_calc

@author: Dylan Neff, Dylan
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import uproot
import awkward as ak
import vector


def main():
    tree_dir = 'F:/Research/BES1_Trees/7GeV/'

    total_events = 0
    res_sum = 0
    for file_name in os.listdir(tree_dir):
        with uproot.open(tree_dir + file_name) as file:
            vector.register_awkward()
            psi_east = file['tree']['psi_east'].array()
            psi_west = file['tree']['psi_west'].array()
            num_events = ak.count(psi_west)
            res_vals = np.cos(2 * (psi_west - psi_east))

            total_events += num_events
            res_sum += np.sum(res_vals)
            # if num_events > 0:
            #     print(file_name, psi_west, psi_east)
            #     print(num_events, res_vals)
    print(np.sqrt(res_sum / total_events), total_events)
    print('donzo')


if __name__ == '__main__':
    main()
