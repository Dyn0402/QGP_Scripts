#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 30 3:30 PM 2022
Created in PyCharm
Created as QGP_Scripts/uproot_mem_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import vector
import uproot

from find_ampt_identical_tracks_uproot import get_event_chunk_indices, check_file_chunks


def main():
    vector.register_awkward()
    # test_path = 'C:/Users/Dylan/Desktop/test.root'
    test_path = 'C:/Users/Dyn04/Desktop/test.root'
    track_attributes = ['pid', 'px', 'py', 'pz']

    # mem_test(test_path, track_attributes)
    print(check_file_chunks(test_path, 'tree', track_attributes, 1, [313, 111], 1e7))


def mem_test(test_path, track_attributes):
    with uproot.open(test_path) as file:
        tracks = file['tree'].arrays(track_attributes)
        print(tracks)
        print(len(tracks))
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        print(len(ak.flatten(tracks[:200], axis=1)))

        track_counts = [len(x) for x in tracks]
        track_counts2 = [len(x)**2 for x in tracks]
        max_tracks = 6000**2

        event_slices = get_event_chunk_indices(track_counts2, max_tracks)
        print(event_slices)

        chunk_counts = [sum(track_counts[s]) for s in event_slices]
        print(chunk_counts)
        print(sum(chunk_counts))
        print(sum(track_counts))

        chunk_counts2 = [sum(np.array(track_counts[s], dtype=np.longlong)**2) for s in event_slices]
        print([np.array(track_counts[s]) ** 2 for s in event_slices])
        print(chunk_counts2)
        print(sum(chunk_counts2))
        print(sum(track_counts2))

        print(track_counts)
        print(track_counts2)
        # return
        # event_index_low, event_index_high = 0, 0
        #
        # while event_index_high < len(tracks):
        #     track_count = 0
        #     while track_count < max_tracks:
        #         track_count += len(tracks[event_index_high])
        #         if event_index_high + 1 < len(tracks):
        #             event_index_high += 1
        #         else:
        #             break

        # print(event_index, track_count, len(ak.flatten(tracks[:event_index])))
        # tracks_chunk = tracks[:event_index]
        # print('tracks: ', len(tracks), tracks)
        # print('tracks[0]: ', len(tracks[0]),  tracks[0])
        # print('tracks_chunk: ', tracks_chunk)
        # combos = ak.combinations(tracks_chunk, 2)
        # print('combos: ', len(combos), combos)
        # print('combos[0]', len(combos[0]), combos[0])
        # track_a, track_b = ak.unzip(combos)
        # print('track_a: ', len(track_a), track_a)
        # print('track_b: ', len(track_b), track_b)
        # print(tracks[0])
        # print(type(tracks))
        # input()
        # tracks = tracks[abs(tracks.eta) < max_eta]

    print('donzo')


if __name__ == '__main__':
    main()
