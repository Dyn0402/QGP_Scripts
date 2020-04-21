#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 13 8:49 PM 2020
Created in PyCharm
Created as QGP_Scripts/ampt_particle_counter.py

@author: Dylan Neff, dylan
"""


from particletools.tables import PYTHIAParticleData as PythiaData
import os


def main():
    particles = read_particles()
    plot_particles(particles)
    print('donzo')


def read_particles():
    path = '/home/dylan/Research/Ampt_Particle_Counts/62GeV/'
    particles = {}
    pdg = PythiaData()
    files = [x for x in os.listdir(path) if '.out' in x]
    for file in files:
        with open(path+file, 'r') as f:
            lines = f.readlines()
            flag = False
            for line in lines:
                if 'making .root file from .dat file...' in line:
                    flag = True
                elif not flag:
                    continue
                if 'Tree written succesfully' in line:
                    break
                line = line.strip().split(' ')
                if len(line) == 2:
                    name = pdg.name(int(line[0]))
                    if name in particles.keys():
                        particles[name] += int(line[1])
                    else:
                        particles.update({name: int(line[1])})

    return particles


def plot_particles(particles):
    print(particles)
    particles_swap = {value: key for key, value in particles.items()}
    total = sum(particles.values())
    for key in sorted(particles_swap.keys(), reverse=True):
        print(f'{particles_swap[key]}'.ljust(13) + f'{key}'.ljust(9) + f'{float(key)/total*100:.4f}%')


if __name__ == '__main__':
    main()
