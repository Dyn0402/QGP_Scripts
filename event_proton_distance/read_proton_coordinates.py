#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 03 12:38 PM 2023
Created in PyCharm
Created as QGP_Scripts/read_proton_coordinates.py

@author: Dylan Neff, Dylan
"""
import os

import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    Example of reading all sample data files (each a different energy) from the given directory path and binning into
    un-normalized 2D histograms. Example of an 62GeV event plotted to make sure things look right.
    :return:
    """
    # path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Jared_Tejes_Sample_Data/BES1/'
    path = 'N:/UCLA_Microsoft/OneDrive - personalmicrosoftsoftware.ucla.edu/Research/UCLA/Jared_Tejes_Sample_Data' \
               '/BES1/'
    rapidity_bin_edges = np.linspace(-0.5, 0.5, 50, endpoint=True)
    phi_bin_edges = np.linspace(0, 2 * np.pi, 50, endpoint=True)
    binning = [rapidity_bin_edges, phi_bin_edges]

    sample_data = read_files(path, binning)  # Read all files in path and return dictionary {energy: binned_events}
    for energy, binned_protons in sample_data.items():
        # Do math things
        print(f'{energy}GeV: ')

    plot_protons_per_event(sample_data)
    plot_event(sample_data[11][26], binning)
    plt.show()
    print('donzo')


def read_files(path, binning):
    """
    Read all files in path directory into a dictionary.
    :param path: Directory containing sample data
    :param binning: [rapidity_bin_edges, phi_bin_edges]
    :return: Dictionary {energy: binned_events}
    """
    sample_data = {}
    for file_name in os.listdir(path):
        energy = int(file_name.strip('GeV.txt'))
        unbinned_event_list = read_file(f'{path}{file_name}')
        binned_events = [bin_event(event, binning) for event in unbinned_event_list]
        sample_data.update({energy: binned_events})

    return sample_data


def read_file(file_path):
    """
    Read proton distributions from file into a list of events, each of which is a list of [rapidity, phi] pairs
    :param file_path: Path of file
    :return: List of events, each of which is a list of [rapidity, phi] pairs
    """
    events = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        line_index = 0
        event = []
        while line_index < len(lines):
            if 'Event ' in lines[line_index]:
                if len(event) > 1:
                    events.append(event)
                event = []
            else:
                line = lines[line_index].strip().split(', ')
                if len(line) == 2:
                    try:
                        event.append([float(line[0]), float(line[1])])
                    except ValueError:
                        pass
            line_index += 1

    return events


def bin_event(event, binning):
    """
    Given event (list) of [rapidity, phi] pairs, create 2D numpy array based on input binning
    :param event: List of [rapidity, phi] pairs
    :param binning: [rapidity_bin_edges, phi_bin_edges]
    :return: Numpy 2D array of number of protons per bin (not normalized)
    """
    return np.histogram2d(*list(zip(*event)), bins=binning)[0]


def plot_protons_per_event(sample_data):
    """
    Plot protons per event distribution for all energies
    :param sample_data: {energy: binned_events}
    :return:
    """
    plt.figure()
    for energy, binned_events in sample_data.items():
        protons_per_event = [np.sum(event) for event in binned_events]
        plt.hist(protons_per_event, bins=np.arange(0, max(protons_per_event) + 1),
                 histtype='step', label=f'{energy}GeV')
    plt.xlabel('Number of Protons per Event')
    plt.legend()
    plt.tight_layout()


def plot_event(event, binning):
    """
    Plot single event
    :param event: Binned event --> 2D array of number of protons per bin
    :param binning: [rapidity_bin_edges, phi_bin_edges]
    :return:
    """
    plt.figure()
    rapid_min, rapid_max, phi_min, phi_max = binning[0].min(), binning[0].max(), binning[1].min(), binning[1].max()
    plt.imshow(event, origin='lower', cmap='jet', extent=[phi_min, phi_max, rapid_min, rapid_max], aspect='auto')
    plt.colorbar()
    plt.xlabel('Phi')
    plt.ylabel('Rapidity')


if __name__ == '__main__':
    main()
