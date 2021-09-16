#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 06 7:27 PM 2020
Created in PyCharm
Created as QGP_Scripts/AzimuthBinData.py

@author: Dylan Neff, dylan
"""


class AzimuthBinData:
    def __init__(self, div=0, path=''):
        """

        :param div: Azimuthal division size in degrees
        :param path: Path to bin data text file
        """
        self.path = path
        self.div = div

        self.data = {}
        self.ratio_dist = {}
        self.pull_dist = {}

        if path != '':
            self.read_data()

    def read_data(self):
        with open(self.path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip().split('\t')
                total_particles = int(line[0])
                bin_particles_string = line[1].split(' ')
                if len(bin_particles_string) > 0:
                    self.data[total_particles] = []
                i = 0
                for entry in bin_particles_string:
                    entry = entry.split(':')
                    while int(entry[0]) > i:
                        self.data[total_particles].append(0)
                        i += 1
                    self.data[total_particles].append(int(entry[1]))
                    i += 1

    def ratio_trans(self):
        for total_particles in self.data:
            for bin_particles in range(len(self.data[total_particles])):
                ratio = bin_particles / total_particles
                if ratio in self.ratio_dist:
                    self.ratio_dist[ratio] += self.data[total_particles][bin_particles]
                else:
                    self.ratio_dist.update({ratio: self.data[total_particles][bin_particles]})

    def pull_trans(self):
        for total_particles in self.data:
            for bin_particles in range(len(self.data[total_particles])):
                pull = bin_particles - total_particles * (self.div / 360)
                if pull in self.ratio_dist:
                    self.pull_dist[pull] += self.data[total_particles][bin_particles]
                else:
                    self.pull_dist.update({pull: self.data[total_particles][bin_particles]})

    def get_ratio_dist(self, binning=None, norm=False):
        self.ratio_trans()
        if binning is None:
            return self.ratio_dist
        else:
            hist = [0 for i in range(len(binning) - 1)]
            for bin_index in range(len(hist)):
                for val, count in self.ratio_dist.items():
                    if binning[bin_index] <= val < binning[bin_index + 1]:
                        hist[bin_index] += count
            if norm:
                total = sum(hist)
                for bin_index in range(len(hist)):
                    hist[bin_index] /= total
            return hist

    def get_pull_dist(self, binning=None, norm=False):
        self.pull_trans()
        if binning is None:
            return self.pull_dist
        else:
            hist = [0 for i in range(len(binning) - 1)]
            for bin_index in range(len(hist)):
                for val, count in self.pull_dist.items():
                    if binning[bin_index] <= val < binning[bin_index + 1]:
                        hist[bin_index] += count
            if norm:
                total = sum(hist)
                for bin_index in range(len(hist)):
                    hist[bin_index] /= total
            return hist

    def get_dist(self):
        return self.data

    def print_dist(self):
        for total_particles in self.data:
            out_line = f'{total_particles}\t'
            for bin_particles in range(len(self.data[total_particles])):
                out_line += f'{bin_particles}:{self.data[total_particles][bin_particles]} '
            print(out_line)
