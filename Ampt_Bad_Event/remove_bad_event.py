#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 4:27 PM 2021
Created in PyCharm
Created as QGP_Scripts/remove_bad_event

@author: Dylan Neff, dylan
"""

import subprocess as sp
import shutil


def main():
    # bad_file_list_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_gang7GeV.txt'
    # bad_tree_repo = '/home/dylan/Research/Ampt_Bad_Event/'
    # bad_file_list_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    # bad_tree_repo = '/media/ucla/Research/Ampt_Bad_Event/'
    bad_file_list_path = 'F:/Research/Ampt_Bad_Event/bad_ampt_events_min_bias2.txt'
    bad_tree_repo = 'F:/Research/Ampt_Bad_Event/'
    bad_tree_sufx = '_bad'
    fix_tree_sufx = '_fix'
    min_identical = 2
    fix_dataset(bad_file_list_path, bad_tree_repo, bad_tree_sufx, fix_tree_sufx, min_identical, test=False,
                fix_method='fix_tree.cpp', replace_text=False)

    print('donzo')


def fix_dataset(bad_file_list_path, bad_tree_repo, bad_sufx='_bad', fix_sufx='_fix', min_identical=2, test=True,
                fix_method='fix_tree.cpp', replace_text=False):
    """
    Fix bad AMPT dataset given bad_file_list_path text file with bad files. For each bad file in list move file to
    bad_tree_repo directory. There create a new file with all events except bad ones. Move this fixed file back to
    original directory with original name.
    :param bad_file_list_path: Path to list of all bad AMPT events and relevant details
    :param bad_tree_repo: Temporary directory to store bad trees and fix them. Good and bad copies remain here at end
    :param bad_sufx: Suffix to add to bad tree file to indicate it is tree with bad events
    :param fix_sufx: Suffix to add to fixed tree file to indicate it is tree copy with bad events removed
    :param min_identical: Minimum number of identical track pairs needed to consider an event bad
    :param test: Flag to run test mode where initial file is copied (not moved) and fixed file is not copied back
    :param fix_method: Determines how to fix trees. If 'pyroot' use pyroot (need import), else call fix_tree.cpp
    :param replace_text: If true run whatever text replacement is hard-coded in. Used to adjust paths on the fly.
    :return:
    """
    bad_trees = get_bad_event_file(bad_file_list_path, min_identical)
    num_trees = len(bad_trees)
    for tree_num, (tree_path, tree) in enumerate(bad_trees.items()):
        print(f'\n\n Tree {tree_num + 1}/{num_trees}')
        if replace_text:
            tree_path = tree_path.replace('/gpfs01/star/pwg/dneff/data/AMPT/most_central',
                                          'F:/Research/AMPT_Trees/slim_most_central')
            # tree_path = tree_path.replace('/most_central',
            #                               '/slim_most_central')
            tree_path = tree_path.replace('.root', '_protons.root')
        repo_tree_path = move_tree(tree_path, bad_tree_repo, bad_sufx, test)
        print(f'{tree_path} moved to {repo_tree_path}')
        if fix_method.lower() == 'pyroot':
            fix_tree_path = fix_tree_pyroot(tree, repo_tree_path, bad_sufx, fix_sufx)
        else:
            fix_tree_path = fix_tree(tree, repo_tree_path, bad_sufx, fix_sufx)
        print(f'Tree fixed: {fix_tree_path}')
        if not test:
            replace_tree(fix_tree_path, tree_path)
            print(f'Fixed tree replaced {tree_path}\n')

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
            if 'num_identical' in bad_event and bad_event['num_identical'] >= min_identical:
                bad_events.append(bad_event)

    bad_trees = {}
    for event in bad_events:
        if event['path'] in bad_trees:
            bad_trees[event['path']].append(event)
        else:
            bad_trees.update({event['path']: [event]})

    return bad_trees


def move_tree(tree_path, repo_path, sufx, copy=True):
    tree_name = tree_path.split('/')[-1].split('\\')[-1]
    tree_name = tree_name.split('.')
    tree_name = tree_name[0] + sufx + '.' + tree_name[-1]
    repo_tree_path = repo_path + tree_name
    if copy:
        shutil.copy(tree_path, repo_tree_path)
    else:
        shutil.move(tree_path, repo_tree_path)

    return repo_tree_path


def fix_tree_pyroot(tree, bad_tree_path, bad_sufx, fix_sufx):
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


def fix_tree(tree, bad_tree_path, bad_sufx, fix_sufx, tree_name='tree'):
    fix_tree_path = bad_tree_path.replace(bad_sufx, fix_sufx)
    bad_events = ','.join([str(bad_event['event_num']) for bad_event in tree])
    cmd = f'fix_tree.cpp("{bad_tree_path}", "{fix_tree_path}", "{tree_name}", {{{bad_events}}})'
    sp.run(['root', '-b', '-q', '-l', cmd])

    return fix_tree_path


def replace_tree(fixed_tree_path, original_path):
    try:
        shutil.copy(fixed_tree_path, original_path)
    except PermissionError:
        print(f'Not copied: {original_path}')


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
