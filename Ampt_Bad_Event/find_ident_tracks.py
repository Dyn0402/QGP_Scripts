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
import os
import sys
from datetime import datetime
from multiprocessing import Pool

import tqdm
if __name__ == '__main__':
    try:
        import istarmap  # Needed for tqdm, not for rcf
    except ModuleNotFoundError:
        pass
# try:
#     import istarmap  # Needed for tqdm
# except ModuleNotFoundError:
#     import sys
#     print(f'Python path: {sys.path}')
#     sys.path.append('/star/u/dneff/git/QGP_Scripts/Analyzer')
#     import istarmap


def main():
    if len(sys.argv) == 2:
        uproot_finder_rcf(sys.argv[1])
    else:
        print('No input root files! Doing other things.')
        # file_writer()
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


# def pyroot_directly():
#     # Found: AMPT_Trees/min_bias/string_melting/7GeV/data_741821621.root EVENT:1617 - 1618
#     # path = '/home/dylan/Research/AMPT_Trees/min_bias/string_melting/7GeV/'
#     path = '/home/dylan/Research/Ampt_Bad_Event/'
#     max_rap = 1.0
#     exclude_pids = [111, 313]  # Neutral pions/K*
#     tree_num = 0
#     track_attributes = ['pid', 'px', 'py', 'pz']
#     for root_path in os.listdir(path):
#         file = ROOT.TFile(path + root_path, 'READ')
#         tree_num += 1
#         print(f'Tree #{tree_num}: {root_path}')
#         for event in file.tree:
#             for track_index_i in range(len(event.pid) - 1):
#                 vec_i = ROOT.TVector3(event.px[track_index_i], event.py[track_index_i], event.pz[track_index_i])
#                 if abs(vec_i.Eta()) > 1.0 or event.pid[track_index_i] in exclude_pids:
#                     continue
#                 for track_index_j in range(track_index_i + 1, len(event.pid)):
#                     track_eqs = [getattr(event, x)[track_index_i] == getattr(event, x)[track_index_j]
#                                  for x in track_attributes]
#                     if np.prod(track_eqs):
#                         print(f'Identical particles in event_num: {event.event}  -  '
#                               f'track indices: {track_index_i}, {track_index_j}')
#                         print(f'\ttrack {track_index_i}: '
#                               f'{[[x, getattr(event, x)[track_index_i]] for x in track_attributes]}')
#                         print(f'\ttrack {track_index_j}: '
#                               f'{[[x, getattr(event, x)[track_index_j]] for x in track_attributes]}')
#                         # event.Show()

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
    out_file_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_slim_min_bias_test.txt'
    # out_file_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_central.txt'
    # path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/'
    # path = 'F:/Research/AMPT_Trees/min_bias/'
    path = '/media/ucla/Research/AMPT_Trees/min_bias/'
    # path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/string_melting/'
    threads = 12
    event_chunks = True
    max_track_combos = 1e7  # Max number of track combinations for a chunk of events. Needed to keep memory from blowing
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

    jobs = [(root_path, tree_name, track_attributes, max_eta, ignore_pids) for root_path in root_paths]

    jobs_start = datetime.now()
    print(f'Jobs Start {jobs_start}\n')
    bad_trees = []
    func = check_file
    if event_chunks:
        func = check_file_chunks
        jobs = [(*job, max_track_combos) for job in jobs]
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


def uproot_finder_rcf(root_path_list):
    start = datetime.now()
    print(f'Start {start}\n')

    out_file_path = 'out.txt'
    tree_name = 'tree'
    write_mode = 'w'
    max_eta = 1
    ignore_pids = [313, 111]
    max_track_combos = 1e7  # Max number of track combinations for a chunk of events. Needed to keep memory from blowing
    track_attributes = ['pid', 'px', 'py', 'pz']

    with open(root_path_list, 'r') as file:
        root_paths = file.readlines()
    root_paths = [root_path.strip() for root_path in root_paths]

    jobs_start = datetime.now()
    print(f'Jobs Start {jobs_start}\n')
    bad_trees = []
    for root_path in root_paths:
        bad_trees.append(check_file_chunks(root_path, tree_name, track_attributes, max_eta, ignore_pids,
                                           max_track_combos))

    with open(out_file_path, write_mode) as file:
        for bad_tree in bad_trees:
            for bad_event in bad_tree:
                if bad_event is not None:
                    out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
                    file.write(f'{out_line}\n')
        file.write('\nFiles checked:\n')
        file.write('\n'.join(root_paths))

    end = datetime.now()
    print(f'\nEnd {end}')
    print(f'Run time: {end - start}')
    print(f'Job run time: {end - jobs_start}')
    print(f'Files run: {len(root_paths)}')
    print(f'Time per file: {(end - jobs_start) / len(root_paths)}')


# def uproot_finder_rcf():
#     """ Looks like this runs just as fast as compiled ROOT code, just slow 2-p comparisons """
#     start = datetime.now()
#     print(f'Start {start}\n')
#
#     # out_file_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_slim_most_central_new.txt'
#     # out_file_path = 'F:/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
#     # out_file_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
#     out_file_path = '/star/u/dneff/Ampt_Bad_Event/bad_ampt_events_central.txt'
#     # path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/'
#     # path = 'F:/Research/AMPT_Trees/min_bias/'
#     # path = '/media/ucla/Research/AMPT_Trees/min_bias/'
#     path = '/gpfs01/star/pwg/dneff/data/AMPT/most_central/string_melting/'
#     single_events = False
#     tree_name = 'tree'
#     write_mode = 'w'
#     max_eta = 1
#     ignore_pids = [313, 111]
#     track_attributes = ['pid', 'px', 'py', 'pz']
#
#     root_paths = []
#     for root, dirs, files in os.walk(path):
#         for file_name in files:
#             if '.root' not in file_name:
#                 continue
#             root_path = os.path.join(root, file_name)
#             root_paths.append(root_path)
#
#     # n_files = len(root_paths)
#     # print(f'{n_files} files found to be checked.')
#     #
#     # jobs = [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
#     #         num, root_path in enumerate(root_paths)]
#     #
#     # jobs_start = datetime.now()
#     # print(f'Jobs Start {jobs_start}\n')
#     # bad_trees = []
#     # func = check_file
#     # if single_events:
#     #     func = check_file_single_events
#     # with Pool(threads) as pool:
#     #     for bad_tree in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
#     #         bad_trees.append(bad_tree)
#
#     jobs_start = datetime.now()
#     print(f'Jobs Start {jobs_start}\n')
#     bad_trees = []
#     for root_path in root_paths:
#         bad_trees.append(check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids))
#
#     # with Pool(threads) as pool:
#     #     bad_trees = pool.starmap(check_file,
#     #                              [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
#     #                               num, root_path in enumerate(root_paths)])
#
#     with open(out_file_path, write_mode) as file:
#         for bad_tree in bad_trees:
#             for bad_event in bad_tree:
#                 if bad_event is not None:
#                     out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
#                     file.write(f'{out_line}\n')
#
#     end = datetime.now()
#     print(f'\nEnd {end}')
#     print(f'Run time: {end - start}')
#     print(f'Job run time: {end - jobs_start}')
#
#
# def uproot_finder_rcf2(root_paths):
#     """ Looks like this runs just as fast as compiled ROOT code, just slow 2-p comparisons """
#     start = datetime.now()
#     print(f'Start {start}\n')
#
#     # out_file_path = '/home/dylan/Research/Ampt_Bad_Event/bad_ampt_events_slim_most_central_new.txt'
#     # out_file_path = 'F:/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
#     # out_file_path = '/media/ucla/Research/Ampt_Bad_Event/bad_ampt_events_minbias.txt'
#     out_file_path = 'bad_ampt_events_central.txt'
#     # path = '/media/ucla/Research/AMPT_Trees/slim_most_central/string_melting/'
#     # path = 'F:/Research/AMPT_Trees/min_bias/'
#     # path = '/media/ucla/Research/AMPT_Trees/min_bias/'
#     single_events = False
#     tree_name = 'tree'
#     write_mode = 'w'
#     max_eta = 1
#     ignore_pids = [313, 111]
#     track_attributes = ['pid', 'px', 'py', 'pz']
#
#     # n_files = len(root_paths)
#     # print(f'{n_files} files found to be checked.')
#     #
#     # jobs = [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
#     #         num, root_path in enumerate(root_paths)]
#     #
#     # jobs_start = datetime.now()
#     # print(f'Jobs Start {jobs_start}\n')
#     # bad_trees = []
#     # func = check_file
#     # if single_events:
#     #     func = check_file_single_events
#     # with Pool(threads) as pool:
#     #     for bad_tree in tqdm.tqdm(pool.istarmap(func, jobs), total=len(jobs)):
#     #         bad_trees.append(bad_tree)
#
#     jobs_start = datetime.now()
#     print(f'Jobs Start {jobs_start}\n')
#     bad_trees = []
#     for root_path in root_paths:
#         bad_trees.append(check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids))
#
#     # with Pool(threads) as pool:
#     #     bad_trees = pool.starmap(check_file,
#     #                              [(root_path, tree_name, track_attributes, max_eta, ignore_pids, num, n_files) for
#     #                               num, root_path in enumerate(root_paths)])
#
#     with open(out_file_path, write_mode) as file:
#         for bad_tree in bad_trees:
#             for bad_event in bad_tree:
#                 if bad_event is not None:
#                     out_line = ''.join(f'{x}: {y}\t' for x, y in bad_event.items())
#                     file.write(f'{out_line}\n')
#
#     end = datetime.now()
#     print(f'\nEnd {end}')
#     print(f'Run time: {end - start}')
#     print(f'Job run time: {end - jobs_start}')


def check_file(root_path, tree_name, track_attributes, max_eta, ignore_pids):
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


def check_file_chunks(root_path, tree_name, track_attributes, max_eta, ignore_pids, max_track_combos):
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

        track_combos = [len(x)**2 for x in tracks]
        event_chunk_slices = get_event_chunk_indices(track_combos, max_track_combos)
        for event_slice in event_chunk_slices:
            tracks_chunk = tracks[event_slice]
            track_a, track_b = ak.unzip(ak.combinations(tracks_chunk, 2))
            ident_p = track_a == track_b  # Check if momenta same
            ident_pid = ak.where(track_a['pid'] == track_b['pid'], True, False)
            ident = ident_p * ident_pid
            num_ident = [{'event_num': x + 1 + event_slice.start, 'num_identical': y}
                         for x, y in enumerate(ak.sum(ident, axis=1)) if y > 0]
            if len(num_ident) > 0:
                print(f'{root_path} {datetime.now()}')
                id_tracks = track_a[ident]
                unique = np.unique(id_tracks[ak.num(id_tracks) > 0].pid, return_counts=True)
                for event in num_ident:
                    event.update({'path': root_path, 'events': len(tracks), 'unique': unique})
                    print(event)
                    bad_events.append(event)

    return bad_events


# def check_file_chunks(root_path, tree_name, track_attributes, max_eta, ignore_pids, max_tracks):
#     """
#     Pass
#     :param root_path:
#     :param tree_name:
#     :param track_attributes:
#     :param max_eta:
#     :param ignore_pids:
#     :param root_num:
#     :param root_num_total:
#     :return:
#     """
#     vector.register_awkward()
#     bad_events = []
#     # print(f'{root_num}/{root_num_total} {root_path} :\t{datetime.now()}')
#     with uproot.open(root_path) as file:
#         tracks = file[tree_name].arrays(track_attributes)
#         tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
#                         with_name='Momentum3D')
#         tracks = tracks[abs(tracks.eta) < max_eta]
#         for bad_pid in ignore_pids:
#             tracks = tracks[abs(tracks['pid']) != bad_pid]
#         num_ident = []
#         for event_num, event_tracks in enumerate(tracks):
#             # print(type(event_tracks))
#             # print(type(tracks))
#             # print(event_tracks)
#             track_a, track_b = ak.unzip(ak.combinations(event_tracks, 2, axis=0))
#             ident_p = track_a == track_b  # Check if momenta same
#             ident_pid = ak.where(track_a['pid'] == track_b['pid'], True, False)
#             ident = ident_p * ident_pid
#             ident_sum = ak.sum(ident, axis=0)
#             if ident_sum > 0:
#                 num_ident.append({'event_num': event_num, 'num_identical': ident_sum})
#         if len(num_ident) > 0:
#             print(f'{root_path} {datetime.now()}')
#             id_tracks = track_a[ident]  # This doesn't work. Fix
#             unique = np.unique(id_tracks[ak.num(id_tracks) > 0].pid, return_counts=True)
#             for event in num_ident:
#                 event.update({'path': root_path, 'events': len(tracks), 'unique': unique})
#                 print(event)
#                 bad_events.append(event)
#
#     return bad_events


def get_event_chunk_indices(track_counts, max_tracks):
    """
    Get indices to split input track_counts array into chunks no larger than max_tracks.
    Group these events together such that the sum of counts for each chunk is no larger than max_tracks.
    If a single event exceeds max_tracks, chunk that event on its own.
    :param track_counts: Array of track counts per event
    :param max_tracks: Max number of tracks for each chunk
    :return: Slices for each chunk of events
    """
    slices = []
    event_low, event_high = 0, 0
    while event_high < len(track_counts):
        tracks = 0
        while tracks < max_tracks:
            tracks += track_counts[event_high]
            event_high += 1
            if event_high >= len(track_counts):
                event_high += 1
                break
        if event_high - 1 > event_low:  # Only decrement if single event doesn't have more counts than max
            event_high -= 1
        slices.append(slice(event_low, event_high))
        event_low = event_high

    return slices


if __name__ == '__main__':
    main()
