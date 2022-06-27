#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 01 2:17 PM 2021
Created in PyCharm
Created as QGP_Scripts/find_bad_ampt_event

@author: Dylan Neff, dylan
"""

import numpy as np
import uproot
import awkward as ak
import vector
import sys
import os
from datetime import datetime
from multiprocessing import Pool

import tqdm
try:
    import istarmap  # Needed for tqdm
except ModuleNotFoundError:
    import sys
    print(f'Python path: {sys.path}')
    sys.path.append('/star/u/dneff/git/QGP_Scripts/Analyzer')
    import istarmap


def main():
    # file_writer()
    if len(sys.argv) == 2:
        uproot_finder_rcf()
    else:
        uproot_finder()
    print('donzo')


def file_writer():
    out_path = '/home/dylan/Research/Ampt_Bad_Event/root_paths.txt'
    # in_path = '/home/dylan/Research/AMPT_Trees/'
    in_path = '/media/ucla/Research/AMPT_Trees/'
    with open(out_path, 'w') as out_file:
        for root, dirs, files in os.walk(in_path):
            for file_name in files:
                if '.root' not in file_name:
                    continue
                root_path = os.path.join(root, file_name)
                out_file.write(f'{root_path}\n')


def macro_caller():
    """ Still too slow """
    macro_name = 'ampt_identical_track_finder_single.cpp'
    macro_path = f'C:\\Users\\Dylan\\source\\sepos\\Dyn0402\\Root_Macros\\src\\{macro_name}'
    path = 'D:/Research/AMPT_Trees/'
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            os.system(f'root -l -q -b \'{macro_path}(\\"{root_path}\\")\'')


def uproot_finder():
    """ Looks like this runs just as fast as compiled ROOT code, just slow 2-p comparisons """
    start = datetime.now()
    print(f'Start {start}\n')

    # out_file_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_slim_most_central_new.txt'
    # out_file_path = 'F:/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    # out_file_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    out_file_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_central.txt'
    # path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/'
    # path = 'F:/Research/AMPT_Trees/min_bias/'
    # path = '/media/ucla/Research/AMPT_Trees/min_bias/'
    path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/string_melting/'
    threads = 8
    single_events = False
    tree_name = 'tree'
    write_mode = 'w'
    max_eta = 1
    ignore_pids = [313, 111]
    track_attributes = ['pid', 'px', 'py', 'pz']

    root_paths = []
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            root_paths.append(root_path)

    n_files = len(root_paths)
    print(f'{n_files} files found to be checked.')

    jobs = [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
            num, root_path in enumerate(root_paths)]

    jobs_start = datetime.now()
    print(f'Jobs Start {jobs_start}\n')
    bad_trees = []
    func = check_file
    if single_events:
        func = check_file_single_events
    with Pool(threads) as pool:
        for bad_tree in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
            bad_trees.append(bad_tree)

    # with Pool(threads) as pool:
    #     bad_trees = pool.starmap(check_file,
    #                              [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
    #                               num, root_path in enumerate(root_paths)])

    with open(out_file_path, write_mode) as file:
        for bad_tree in bad_trees:
            for bad_event in bad_tree:
                if bad_event is not None:
                    out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
                    file.write(f'{out_line}\n')

    end = datetime.now()
    print(f'\nEnd {end}')
    print(f'Run time: {end - start}')
    print(f'Job run time: {end - jobs_start}')


def uproot_finder_rcf():
    """ Looks like this runs just as fast as compiled ROOT code, just slow 2-p comparisons """
    start = datetime.now()
    print(f'Start {start}\n')

    # out_file_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_slim_most_central_new.txt'
    # out_file_path = 'F:/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    # out_file_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
    out_file_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_central.txt'
    # path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/'
    # path = 'F:/Research/AMPT_Trees/min_bias/'
    # path = '/media/ucla/Research/AMPT_Trees/min_bias/'
    path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/string_melting/'
    single_events = False
    tree_name = 'tree'
    write_mode = 'w'
    max_eta = 1
    ignore_pids = [313, 111]
    track_attributes = ['pid', 'px', 'py', 'pz']

    root_paths = []
    for root, dirs, files in os.walk(path):
        for file_name in files:
            if '.root' not in file_name:
                continue
            root_path = os.path.join(root, file_name)
            root_paths.append(root_path)

    # n_files = len(root_paths)
    # print(f'{n_files} files found to be checked.')
    #
    # jobs = [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
    #         num, root_path in enumerate(root_paths)]
    #
    # jobs_start = datetime.now()
    # print(f'Jobs Start {jobs_start}\n')
    # bad_trees = []
    # func = check_file
    # if single_events:
    #     func = check_file_single_events
    # with Pool(threads) as pool:
    #     for bad_tree in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
    #         bad_trees.append(bad_tree)

    jobs_start = datetime.now()
    print(f'Jobs Start {jobs_start}\n')
    bad_trees = []
    for root_path in root_paths:
        bad_trees.append(check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids))

    # with Pool(threads) as pool:
    #     bad_trees = pool.starmap(check_file,
    #                              [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
    #                               num, root_path in enumerate(root_paths)])

    with open(out_file_path, write_mode) as file:
        for bad_tree in bad_trees:
            for bad_event in bad_tree:
                if bad_event is not None:
                    out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
                    file.write(f'{out_line}\n')

    end = datetime.now()
    print(f'\nEnd {end}')
    print(f'Run time: {end - start}')
    print(f'Job run time: {end - jobs_start}')


def check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids, root_num=-1, root_num_total=-1):
    vector.register_awkward()
    bad_events = []
    # print(f'{root_num}/{root_num_total} {root_path} :\t{datetime.now()}')
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        tracks = tracks[abs(tracks.eta) < max_eta]
        for bad_pid in ignore_pids:
            tracks = tracks[abs(tracks['pid']) != bad_pid]
        track_a, track_b = ak.unzip(ak.combinations(tracks, 2))
        ident_p = track_a == track_b  # Check if momenta same
        ident_pid = ak.where(track_a['pid'] == track_b['pid'], True, False)
        ident = ident_p * ident_pid
        num_ident = [{'event_num': x + 1, 'num_identical': y} for x, y in enumerate(ak.sum(ident, axis=1)) if y > 0]
        if len(num_ident) > 0:
            print(f'{root_path} {datetime.now()}')
            id_tracks = track_a[ident]
            unique = np.unique(id_tracks[ak.num(id_tracks) > 0].pid, return_counts=True)
            for event in num_ident:
                event.update({'path': root_path, 'events': len(tracks), 'unique': unique})
                print(event)
                bad_events.append(event)

    return bad_events


def check_file_single_events(root_path, tree_name, track_attributes, max_eta, ignore_pids, root_num=-1,
                             root_num_total=-1):
    """
    Not fully implemented! Need to fix last part
    :param root_path:
    :param tree_name:
    :param track_attributes:
    :param max_eta:
    :param ignore_pids:
    :param root_num:
    :param root_num_total:
    :return:
    """
    vector.register_awkward()
    bad_events = []
    # print(f'{root_num}/{root_num_total} {root_path} :\t{datetime.now()}')
    with uproot.open(root_path) as file:
        tracks = file[tree_name].arrays(track_attributes)
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        tracks = tracks[abs(tracks.eta) < max_eta]
        for bad_pid in ignore_pids:
            tracks = tracks[abs(tracks['pid']) != bad_pid]
        num_ident = []
        for event_num, event_tracks in enumerate(tracks):
            # print(type(event_tracks))
            # print(type(tracks))
            # print(event_tracks)
            track_a, track_b = ak.unzip(ak.combinations(event_tracks, 2, axis=0))
            ident_p = track_a == track_b  # Check if momenta same
            ident_pid = ak.where(track_a['pid'] == track_b['pid'], True, False)
            ident = ident_p * ident_pid
            ident_sum = ak.sum(ident, axis=0)
            if ident_sum > 0:
                num_ident.append({'event_num': event_num, 'num_identical': ident_sum})
        if len(num_ident) > 0:
            print(f'{root_path} {datetime.now()}')
            id_tracks = track_a[ident]  # This doesn't work. Fix
            unique = np.unique(id_tracks[ak.num(id_tracks) > 0].pid, return_counts=True)
            for event in num_ident:
                event.update({'path': root_path, 'events': len(tracks), 'unique': unique})
                print(event)
                bad_events.append(event)

    return bad_events


if __name__ == '__main__':
    main()
