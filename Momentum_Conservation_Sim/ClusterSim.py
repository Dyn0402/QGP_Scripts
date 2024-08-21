#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 04 11:12 AM 2024
Created in PyCharm
Created as QGP_Scripts/ClusterSim.py

@author: Dylan Neff, Dylan
"""

import numpy as np

from PConsSim import rotate_vector


class ClusterSim:
    def __init__(self, tracks, p_clust=0, frac_rotate=0.9, rng=None):
        self.tracks = tracks.copy()
        self.n_tracks = len(tracks)
        self.dim = len(tracks[0])
        self.p_clust = p_clust
        self.frac_rotate = frac_rotate

        if rng is None:
            self.rng = np.random.default_rng()
        else:
            self.rng = rng
        self.initial_rng_state = self.rng.bit_generator.state

    def cluster_tracks(self):
        cluster_vec = self.rng.uniform(-1, 1, size=self.dim)
        for track_i, track in enumerate(self.tracks):
            if self.rng.random() < self.p_clust:
                self.tracks[track_i] = rotate_vector(track, cluster_vec, self.frac_rotate)
