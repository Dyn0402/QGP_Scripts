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


def main():
    vector.register_awkward()
    test_path = 'C:/Users/Dylan/Desktop/test.root'
    track_attributes = ['pid', 'px', 'py', 'pz']
    with uproot.open(test_path) as file:
        tracks = file['tree'].arrays(track_attributes)
        print(tracks)
        print(len(tracks))
        tracks = ak.zip({'pid': tracks['pid'], 'px': tracks['px'], 'py': tracks['py'], 'pz': tracks['pz']},
                        with_name='Momentum3D')
        print(len(ak.flatten(tracks[:200], axis=1)))

        track_chunk = tracks[:]

        track_counts = [len(x) for x in tracks]
        max_tracks = 10000

        event_slices = get_event_chunk_indices(track_counts, max_tracks)
        print(event_slices)
        print([sum(track_counts[s]) for s in event_slices])

        print([len(x) for x in tracks])
        return
        event_index_low, event_index_high = 0, 0

        while event_index_high < len(tracks):
            track_count = 0
            while track_count < max_tracks:
                track_count += len(tracks[event_index_high])
                if event_index_high + 1 < len(tracks):
                    event_index_high += 1
                else:
                    break

        # print(event_index, track_count, len(ak.flatten(tracks[:event_index])))
        tracks_chunk = tracks[:event_index]
        print('tracks: ', len(tracks), tracks)
        print('tracks[0]: ', len(tracks[0]),  tracks[0])
        print('tracks_chunk: ', tracks_chunk)
        combos = ak.combinations(tracks_chunk, 2)
        print('combos: ', len(combos), combos)
        print('combos[0]', len(combos[0]), combos[0])
        track_a, track_b = ak.unzip(combos)
        print('track_a: ', len(track_a), track_a)
        print('track_b: ', len(track_b), track_b)
        print(tracks[0])
        print(type(tracks))
        # input()
        # tracks = tracks[abs(tracks.eta) < max_eta]
    print('donzo')


def get_event_chunk_indices(track_counts, max_tracks):
    slices = []
    event_low, event_high = 0, 0
    while event_high < len(track_counts):
        tracks = 0
        while tracks < max_tracks:
            tracks += track_counts[event_high]
            event_high += 1
            if event_high >= len(track_counts):
                break
        # Forgot event_low
        slices.append(slice(event_low, event_high))

    return slices



if __name__ == '__main__':
    main()
