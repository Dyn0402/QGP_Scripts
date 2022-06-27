#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 06 5:20 PM 2021
Created in PyCharm
Created as QGP_Scripts/remove_bad_event_copy

@author: Dylan Neff, dylan
"""

import os
import subprocess as sp
import shutil


def main():
    """Same as remove_bad_event but no actual tree editing so it doesn't have ROOT dependence.
    Used to copy already fixed files into correct directories in RCF."""
    bad_file_list_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    # bad_file_list_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_gang7GeV.txt'
    bad_tree_repo = '/star/u/dneff/Ampt_Bad_Event/'
    bad_tree_sufx = '_bad'
    fix_tree_sufx = '_fix'
    min_identical = 2
    bad_trees = get_bad_event_file(bad_file_list_path, min_identical)
    for tree_path, tree in bad_trees.items():
        tree_path = tree_path.replace('/media/ucla/Research/AMPT_Trees', '/gpfs01/star/pwg/dneff/data/AMPT')
        # tree_path = tree_path.replace('D:/Research/AMPT_Trees', '/gpfs01/star/pwg/dneff/data/AMPT')
        print(tree_path, tree)
        repo_tree_path = move_tree(tree_path, bad_tree_repo, bad_tree_sufx)
        fix_tree_path = fix_tree(tree, repo_tree_path, bad_tree_sufx, fix_tree_sufx)
        # fix_tree_path = repo_tree_path.replace(bad_tree_sufx, fix_tree_sufx)
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


def replace_tree(fixed_tree_path, original_path):
    try:
        shutil.copy(fixed_tree_path, original_path)
    except PermissionError:
        print(f'Not copied: {original_path}')


def fix_tree(tree, bad_tree_path, bad_sufx, fix_sufx, tree_name='tree'):
    fix_tree_path = bad_tree_path.replace(bad_sufx, fix_sufx)
    bad_events = ','.join([str(bad_event['event_num']) for bad_event in tree])
    cmd = f'fix_tree.cpp("{bad_tree_path}", "{fix_tree_path}", "{tree_name}", {{{bad_events}}})'
    sp.run(['root', '-b', '-q', '-l', cmd])

    return fix_tree_path


if __name__ == '__main__':
    main()
