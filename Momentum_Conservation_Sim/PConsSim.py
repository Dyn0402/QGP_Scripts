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


def main():
    PConsSim(10, 5)
    print('donzo')


class PConsSim:
    def __init__(self, energy, n_particles, dimension=3, rng=None):
        self.energy = energy
        self.n_particles = n_particles
        self.dim = dimension

        self.alpha = 0.6
        self.alpha_upper_bound = 0.9
        self.alpha_lower_bound = 0.01
        self.alphas = []
        self.ratios = []
        self.iterations = 4
        self.max_iterations = 200
        self.convergence_momentum = None
        self.adaptive_alpha = False

        if rng is None:
            self.rng = np.random.default_rng()
        else:
            self.rng = rng
        self.initial_rng_state = self.rng.bit_generator.state

        self.momenta = None
        self.net_momentum_vec = None
        self.init_net_momentum = None
        self.net_momentum = None
        self.net_momentum_iterations = None

        self.gen_tracks()

    def gen_tracks(self):
        self.momenta = self.rng.uniform(-self.energy, self.energy, size=(self.n_particles, self.dim))
        self.net_momentum_vec = np.sum(self.momenta, axis=0)
        self.init_net_momentum = np.linalg.norm(self.net_momentum_vec)
        self.net_momentum = self.init_net_momentum
        self.net_momentum_iterations = [self.net_momentum]

    def rotate_tracks(self):
        iteration = 0
        stop = False
        while not stop:
            angle_frac = self.alpha * self.net_momentum / self.init_net_momentum / np.sqrt(self.n_particles)
            for j in range(len(self.momenta)):
                self.momenta[j] = rotate_vector(self.momenta[j], -self.net_momentum_vec, angle_frac)

            self.net_momentum_vec = np.sum(self.momenta, axis=0)
            self.net_momentum = np.linalg.norm(self.net_momentum_vec)
            self.net_momentum_iterations.append(self.net_momentum)
            if self.adaptive_alpha:
                success = self.update_alpha()
                if not success:
                    break
            iteration += 1
            if self.convergence_momentum is not None:
                if self.net_momentum <= self.convergence_momentum:
                    stop = True
                elif iteration > self.max_iterations:
                    stop = True
                    print(f'Event Convergence Failed, net_p={self.net_momentum}')
                    print(self.alphas)
            else:
                if iteration >= self.iterations:
                    stop = True

    def update_alpha(self):
        ratio = self.net_momentum_iterations[-1] / self.net_momentum_iterations[-2]
        if len(self.ratios) == 0:
            if ratio >= 1:
                new_alpha = self.alpha_lower_bound
            else:
                new_alpha = self.alpha - 0.1
        else:
            if self.alpha <= self.alpha_lower_bound:
                new_alpha = self.alpha + 0.01
            elif self.alpha >= self.alpha_upper_bound:
                new_alpha = self.alpha - 0.05
            else:
                dalpha = self.alpha - self.alphas[-1]
                dratio = ratio - self.ratios[-1]
                new_alpha = self.alpha - 0.05 * np.sign(dratio) * np.sign(dalpha)
            if new_alpha < self.alpha_lower_bound:
                new_alpha = self.alpha_lower_bound
            elif new_alpha > self.alpha_upper_bound:
                new_alpha = self.alpha_upper_bound
        self.ratios.append(ratio)
        self.alphas.append(self.alpha)
        self.alpha = new_alpha
        return True

    def rotate_tracks_debug(self):
        for i in range(self.iterations):
            angle_frac = self.angle_frac * self.net_momentum / self.init_net_momentum
            # angle_frac = 0.1
            print(f'\n\nangle_frac: {angle_frac}')
            for j in range(len(self.momenta)):
                print(f'Pre: momentum {j}: {self.momenta[j]}, neg_net_p: {-self.net_momentum_vec}')
                dot = np.dot(self.momenta[j] / np.linalg.norm(self.momenta[j]),
                             -self.net_momentum_vec / np.linalg.norm(-self.net_momentum_vec))
                print(f'Pre-dot: {dot}, angle: {np.rad2deg(np.arccos(dot))}')
                self.momenta[j] = rotate_vector_debug(self.momenta[j], -self.net_momentum_vec, angle_frac)
                print(f'Post: momentum {j}: {self.momenta[j] / np.linalg.norm(self.momenta[j])}, '
                      f'{-self.net_momentum_vec / np.linalg.norm(self.net_momentum_vec)}')
                dot = np.dot(self.momenta[j] / np.linalg.norm(self.momenta[j]),
                             -self.net_momentum_vec / np.linalg.norm(-self.net_momentum_vec))
                print(f'Post-dot: {dot}, angle: {np.rad2deg(np.arccos(dot))}\n')

            self.net_momentum_vec = np.sum(self.momenta, axis=0)
            self.net_momentum = np.linalg.norm(self.net_momentum_vec)
            self.net_momentum_iterations.append(self.net_momentum)
            # print(self.net_momentum)

    def get_m_tracks(self, m):
        random_track_indices = self.rng.choice(self.momenta.shape[0], size=m, replace=False)
        return self.momenta[random_track_indices, :]


def rotate_vector(vec, vec_target, angle_fraction=None):
    vec_unit = vec / np.linalg.norm(vec)
    vec_target_unit = vec_target / np.linalg.norm(vec_target)

    # Find the angle and axis of rotation
    angle = np.arccos(np.dot(vec_unit, vec_target_unit))
    angle = angle * angle_fraction if angle_fraction is not None else angle
    axis = np.cross(vec_unit, vec_target_unit)
    axis = axis / np.linalg.norm(axis)

    # Create the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    rot_matrix = np.array(
        [[axis[0] ** 2 * C + c, axis[0] * axis[1] * C - axis[2] * s, axis[0] * axis[2] * C + axis[1] * s],
         [axis[1] * axis[0] * C + axis[2] * s, axis[1] ** 2 * C + c, axis[1] * axis[2] * C - axis[0] * s],
         [axis[2] * axis[0] * C - axis[1] * s, axis[2] * axis[1] * C + axis[0] * s, axis[2] ** 2 * C + c]])

    # Apply the rotation matrix to vec1
    rotated_vec1 = np.dot(rot_matrix, vec)

    return rotated_vec1


def rotate_vector_debug(vec, vec_target, angle_fraction=None):
    vec_unit = vec / np.linalg.norm(vec)
    vec_target_unit = vec_target / np.linalg.norm(vec_target)

    # Find the angle and axis of rotation
    angle = np.arccos(np.dot(vec_unit, vec_target_unit))
    angle = angle * angle_fraction if angle_fraction is not None else angle
    axis = np.cross(vec_unit, vec_target_unit)
    axis = axis / np.linalg.norm(axis)
    print(f'angle={angle}, axis={axis}')

    # Create the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    rot_matrix = np.array(
        [[axis[0] ** 2 * C + c, axis[0] * axis[1] * C - axis[2] * s, axis[0] * axis[2] * C + axis[1] * s],
         [axis[1] * axis[0] * C + axis[2] * s, axis[1] ** 2 * C + c, axis[1] * axis[2] * C - axis[0] * s],
         [axis[2] * axis[0] * C - axis[1] * s, axis[2] * axis[1] * C + axis[0] * s, axis[2] ** 2 * C + c]])

    # Apply the rotation matrix to vec1
    rotated_vec1 = np.dot(rot_matrix, vec)

    return rotated_vec1


if __name__ == '__main__':
    main()
