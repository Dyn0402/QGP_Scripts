#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 19 6:15 PM 2022
Created in PyCharm
Created as QGP_Scripts/PConsSim

@author: Dylan Neff, Dylan
"""


import numpy as np
from scipy.stats import norm
# from scipy.spatial.transform import Rotation as R

from momentum_conservation_model import rotate_vector


def main():
    PConsSim(10, 5)
    print('donzo')


class PConsSim:
    def __init__(self, energy, n_particles, dimension=3):
        self.energy = energy
        self.n_particles = n_particles
        self.dim = dimension

        self.angle_frac = 0.8 / np.sqrt(self.n_particles)
        self.iterations = 4

        self.momenta = np.random.uniform(-self.energy, self.energy, size=(self.n_particles, self.dim))
        self.net_momentum_vec = np.sum(self.momenta, axis=0)
        self.init_net_momentum = np.linalg.norm(self.net_momentum_vec)
        self.net_momentum = self.init_net_momentum

    def rotate_tracks(self):
        for i in range(self.iterations):
            for j in range(len(self.momenta)):
                angle_frac = self.angle_frac * self.net_momentum / self.init_net_momentum
                self.momenta[j] = rotate_vector(self.momenta[j], -self.net_momentum_vec, angle_frac)

            self.net_momentum_vec = np.sum(self.momenta, axis=0)
            self.net_momentum = np.linalg.norm(self.net_momentum_vec)

    def get_m_tracks(self, m):
        random_track_indices = np.random.choice(self.momenta.shape[0], size=m, replace=False)
        return self.momenta[random_track_indices, :]


if __name__ == '__main__':
    main()
