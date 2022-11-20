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
from scipy.spatial.transform import Rotation as R


def main():
    PConsSim(10, 5, 3)
    print('donzo')


class PConsSim:
    def __init__(self, energy, particles, dimension=3):
        self.energy = energy
        self.particles = particles
        self.dim = dimension

        energies = np.random.random(particles)
        energies = energies / sum(energies) * energy

        p_unit_vectors = norm.rvs(size=(particles, self.dim))  # Need to check for 0,0,0 s
        p_unit_norms = np.sqrt(np.sum(p_unit_vectors**2, axis=1))
        p_unit_vectors /= np.tile(np.array([p_unit_norms]).transpose(), (1, dimension))

        print(energies)
        print(sum(energies))
        print(p_unit_vectors)
        print(np.sqrt(np.sum(p_unit_vectors**2, axis=1)))

        ps = p_unit_vectors * np.tile(np.array([energies]).transpose(), (1, dimension))
        print(ps)

        p_net = np.sum(ps, axis=0)
        print(p_net)

        p_net_uv = p_net / np.sqrt(np.sum(p_net**2))

        for i in range(particles):
            particle = i  # np.random.randint(0, particles)
            print(f'\nparticle {i}:')
            p, p_uv = ps[particle], p_unit_vectors[particle]
            e = energies[particle]
            theta = np.arccos(np.dot(p_net, p) / (np.linalg.norm(p_net) * np.linalg.norm(p)))
            p_perp_uv = np.cross(p_net_uv, p_uv)
            print([np.cos(theta / 2), *(np.sin(theta / 2) * p_perp_uv)])
            rnp = R.from_quat([*(np.sin(theta / 2) * p_perp_uv), np.cos(theta / 2)])
            r = rotation_matrix_3d(*p_perp_uv, theta)
            print(p_net, p, e, np.rad2deg(theta))
            p_rot = np.dot(p.transpose(), r)
            p_rotnp = rnp.apply(p)
            print('rotated p: ',  p_rot)
            print('rotated pnp: ',  p_rotnp)
            theta_new = np.arccos(np.dot(p_net, p_rot) / (np.linalg.norm(p_net) * np.linalg.norm(p_rot)))
            theta_newnp = np.arccos(np.dot(p_net, p_rotnp) / (np.linalg.norm(p_net) * np.linalg.norm(p_rotnp)))
            print(np.rad2deg(theta_new))
            print(np.rad2deg(theta_newnp))


def rotation_matrix_3d(ux, uy, uz, theta):
    c, s, omc = np.cos(theta), np.sin(theta), 1 - np.cos(theta)
    r1 = (c + ux**2 * omc, ux * uy * omc - uz * s, ux * uz * omc + uy * s)
    r2 = (uy * ux * omc + uz * s, c + uy**2 * omc, uy * ux * omc - ux * s)
    r3 = (uz * ux * omc - uy * s, uz * uy * omc + ux * s, c + uz**2 * omc)

    return np.array([r1, r2, r3])


if __name__ == '__main__':
    main()
