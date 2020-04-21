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
from matplotlib.dates import date2num, DateFormatter, HourLocator
from datetime import datetime


def main():
    trees_path = '/gpfs01/star/pwg/dneff/scratch/ampt/output/'
    energies = [11, 39]  # [7, 11, 19, 27, 39, 62]
    event_time_data = {'total': [[], []]}
    print()
    for energy in energies:
        events = []
        times = []
        path = f'{trees_path}{energy}GeV/'
        files = os.listdir(path)
        print(f'Reading ~{len(files)} {energy}GeV trees...')
        path_events_data = get_events(path)
        print(f'  Getting times for each {energy}GeV tree...')
        for file in files:
            if '.root' in file:
                try:
                    events.append(path_events_data[file])
                    times.append(datetime.fromtimestamp(os.path.getmtime(path + file)))
                except KeyError:
                    pass
        times, events = zip(*sorted(zip(times, events)))
        event_time_data.update({energy: [list(times), list(events)]})
        event_time_data['total'][0] += times
        event_time_data['total'][1] += events

    times, events = zip(*sorted(zip(event_time_data['total'][0], event_time_data['total'][1])))
    event_time_data[0] = times
    event_time_data[1] = events

    plot_event_time_data(event_time_data, energies)

    # print(get_events('/media/dylan/SSD_Storage/Research/Trees_Ampt/7.root'))
    print('donzo')


def get_events(path):
    # root_path = '/home/dylan/git/Research/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    root_path = '/star/u/dneff/git/Ampt_Runner/Root_Macros/get_tree_events.cpp'
    command = f'{root_path}("{path}")'
    res = sp.check_output(['root', '-q', '-b', command]).decode('utf-8')
    event_file_data = {}
    for line in res.strip().split('\n'):
        if 'Number of events in tree:' in line:
            line = line.strip().split(': ')[-1]
            events = int(line.split(' ')[-1])
            file = line.split(' ')[-2].split('/')[-1]
            event_file_data.update({file: events})

    return event_file_data


def plot_event_time_data(data, energies):
    formatter = DateFormatter('%d %h')
    loc = HourLocator(interval=12)

    plt.figure(1)
    data_total = data['total']
    events_total = get_cumulative(data_total[1])
    dates_total = date2num(data_total[0])
    plt.plot_date(dates_total, events_total, fmt='-')
    plt.ylabel('Events Produced')
    plt.title('Total AMPT Events Produced')
    plt.grid()
    plt.gca().xaxis.set_major_locator(loc)
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().xaxis.set_tick_params(rotation=30, labelsize=10)

    plt.figure(2)
    colors = {7: '#1f77b4', 11: '#ff7f0e', 19: '#2ca02c', 27: '#d62728', 39: '#9467bd', 62: '#8c564b'}
    plt.title('AMPT Events Produced by Energy')
    plt.ylabel('Events Produced')
    plt.grid()
    plt.gca().xaxis.set_major_locator(loc)
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().xaxis.set_tick_params(rotation=30, labelsize=10)
    for energy in energies:
        energy_events = get_cumulative(data[energy][1])
        energy_dates = date2num(data[energy][0])
        plt.plot_date(energy_dates, energy_events, label=f'{energy}GeV', c=colors[energy], fmt='-')
    plt.legend()

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
