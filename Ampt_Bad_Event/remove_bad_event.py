#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 4:27 PM 2021
Created in PyCharm
Created as QGP_Scripts/remove_bad_event

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT


def main():
    # Found: AMPT_Trees/min_bias/string_melting/7GeV/data_741821621.root EVENT:1617 - 1618
    path = '/home/dylan/Research/Ampt_Bad_Event/data_741821621_bad.root'
    new_path = '/home/dylan/Research/Ampt_Bad_Event/data_741821621.root'
    bad_event = 1618
    bad_file = ROOT.TFile(path, 'READ')
    new_file = ROOT.TFile(new_path, 'RECREATE')
    bad_tree = bad_file.tree
    new_tree = bad_tree.CloneTree(0)
    for event in bad_tree:
        if event.event != bad_event:
            new_tree.Fill()
    bad_file.Close()
    new_file.Write()
    new_file.Close()
    print('donzo')


if __name__ == '__main__':
    main()
