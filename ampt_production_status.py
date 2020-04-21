#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 20 5:33 PM 2020
Created in PyCharm
Created as QGP_Scripts/ampt_production_status.py

@author: Dylan Neff, dylan
"""

import os
import subprocess as sp
import matplotlib.pyplot as plt


def main():
    trees_path = '/gpfs01/star/pwg/dneff/scratch/ampt/output/'
    energies = [7, 11, 19, 27, 39, 62]
    event_time_data = {'total': [[], []]}
    for energy in energies:
        events = []
        times = []
        path = f'{trees_path}{energy}GeV/'
        files = os.listdir(path)
        for file in files:
            if '.root' in file:
                events.append(get_events(path+file))
                times.append(os.path.getmtime(path+file))
        events, times = zip(*sorted(zip(times, events)))
        event_time_data.update({energy: [list(times), list(events)]})
        event_time_data['total'][0] += times
        event_time_data['total'][1] += events

    plot_event_time_data(event_time_data, energies)

    # print(get_events('/media/dylan/SSD_Storage/Research/Trees_Ampt/7.root'))
    print('donzo')


def get_events(path):
    # root_path = '/home/dylan/git/Research/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    root_path = '/star/u/dneff/git/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    command = f'{root_path}("{path}")'
    res = sp.check_output(['root', '-q', '-b', command]).decode('utf-8')
    events = -1
    for line in res.strip().split('\n'):
        if 'Number of events in tree:' in line:
            events = int(line.strip().split(': ')[-1])

    return events


def plot_event_time_data(data, energies):
    data_total = data['total']
    events_total = get_cumulative(data_total[1])
    plt.plot(data_total[0], events_total)
    plt.show()


def get_cumulative(x):
    x_cum = []
    x_last = 0
    for xi in x:
        x_last += xi
        x_cum.append(x_last)

    return x_cum


if __name__ == '__main__':
    main()
