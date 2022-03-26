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
import os
import shutil


def main():
    bad_file_list_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_gang7GeV.txt'
    bad_tree_repo = '/home/dylan/Research/Ampt_Bad_Event/'
    bad_tree_sufx = '_bad'
    fix_tree_sufx = '_fix'
    min_identical = 2
    bad_trees = get_bad_event_file(bad_file_list_path, min_identical)
    for tree_path, tree in bad_trees.items():
        # tree_path = tree_path.replace('/media/ucla/Research/AMPT_Trees', '/gpfs01/star/pwg/dneff/data/AMPT')
        tree_path = tree_path.replace('D:/Research', '/media/ucla/Research')
        print(tree_path, tree)
        repo_tree_path = move_tree(tree_path, bad_tree_repo, bad_tree_sufx)
        fix_tree_path = fix_tree(tree, repo_tree_path, bad_tree_sufx, fix_tree_sufx)
        replace_tree(fix_tree_path, tree_path)

    print('donzo')


def get_bad_event_file(path, min_identical):
    bad_events = []
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            bad_event = {}
            line = line.strip().split('\t')
            for item in line:
                item = item.strip().split(': ')
                try:
                    bad_event.update({item[0]: int(item[-1])})
                except ValueError:
                    bad_event.update({item[0]: item[-1]})
            if bad_event['num_identical'] >= min_identical:
                bad_events.append(bad_event)

    bad_trees = {}
    for event in bad_events:
        if event['path'] in bad_trees:
            bad_trees[event['path']].append(event)
        else:
            bad_trees.update({event['path']: [event]})

    return bad_trees


def move_tree(tree_path, repo_path, sufx):
    tree_name = tree_path.split('/')[-1]
    tree_name = tree_name.split('.')
    tree_name = tree_name[0] + sufx + '.' + tree_name[-1]
    repo_tree_path = repo_path + tree_name
    shutil.move(tree_path, repo_tree_path)

    return repo_tree_path


def fix_tree(tree, bad_tree_path, bad_sufx, fix_sufx):
    fix_tree_path = bad_tree_path.replace(bad_sufx, fix_sufx)
    bad_file = ROOT.TFile(bad_tree_path, 'READ')
    new_file = ROOT.TFile(fix_tree_path, 'RECREATE')
    bad_tree = bad_file.tree
    new_tree = bad_tree.CloneTree(0)
    for event in bad_tree:
        for bad_event in tree:
            if event.event != bad_event['event_num']:
                new_tree.Fill()
    bad_file.Close()
    new_file.Write()
    new_file.Close()

    return fix_tree_path


def replace_tree(fixed_tree_path, original_path):
    shutil.copy(fixed_tree_path, original_path)


def remove_specific():
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


if __name__ == '__main__':
    main()
