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
from matplotlib.dates import date2num
from datetime import datetime


def main():
    trees_path = '/gpfs01/star/pwg/dneff/scratch/ampt/output/'
    energies = [11]  # [7, 11, 19, 27, 39, 62]
    event_time_data = {'total': [[], []]}
    for energy in energies:
        events = []
        times = []
        path = f'{trees_path}{energy}GeV/'
        files = os.listdir(path)
        print(f'Reading ~{len(files)} {energy}GeV trees...')
        print(get_all_events(path))
        # for file in files:
        #     if '.root' in file:
        #         print(f'  Opening {path+file}')
        #         events.append(get_events(path+file))
        #         times.append(datetime.fromtimestamp(os.path.getmtime(path+file)))
        # times, events = zip(*sorted(zip(times, events)))
        # event_time_data.update({energy: [list(times), list(events)]})
        # event_time_data['total'][0] += times
        # event_time_data['total'][1] += events

    # plot_event_time_data(event_time_data, energies)

    # print(get_events('/media/dylan/SSD_Storage/Research/Trees_Ampt/7.root'))
    print('donzo')


def get_events(path):
    # root_path = '/home/dylan/git/Research/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    root_path = '/star/u/dneff/git/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    command = f'{root_path}("{path}")'
    res = sp.check_output(['root', '-q', '-b', command], timeout=10).decode('utf-8')
    events = 0
    for line in res.strip().split('\n'):
        if 'Number of events in tree:' in line:
            events = int(line.strip().split(': ')[-1])

    return events


def get_all_events(path):
    # root_path = '/home/dylan/git/Research/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    root_path = '/star/u/dneff/git/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    command = f'{root_path}("{path}")'
    res = sp.check_output(['root', '-q', '-b', command], timeout=10).decode('utf-8')
    event_file_data = [[], []]
    for line in res.strip().split('\n'):
        if 'Number of events in tree:' in line:
            events = int(line.strip().split(': ')[-1])
            file = line.strip().split(': ')[-2]
            event_file_data[0].append(file)
            event_file_data[1].append(events)

    return event_file_data


def plot_event_time_data(data, energies):
    print('Plotting total')
    data_total = data['total']
    events_total = get_cumulative(data_total[1])
    dates_total = date2num(data_total[0])
    plt.plot(dates_total, events_total)
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
