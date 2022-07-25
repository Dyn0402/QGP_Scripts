#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 25 10:57 AM 2022
Created in PyCharm
Created as QGP_Scripts/ampt_production_status_uproot.py

@author: Dylan Neff, Dylan
"""

import os
from datetime import datetime

import uproot

from ampt_production_status import plot_event_time_data


def main():
    trees_path = '/gpfs01/star/pwg/dneff/data/AMPT/dylan_run/output/'
    energies = [7, 11, 19, 27, 39, 62, '2-7TeV_PbPb']
    energies_found = []
    event_time_data = {}
    total_events = []
    total_times = []
    print()
    for energy in energies:
        if type(energy) == int:
            energy = f'{energy}GeV'
        path = f'{trees_path}{energy}/'
        files = os.listdir(path)
        if len(files) <= 0:
            print(f'{len(files)} {energy} trees, skipping energy')
            continue
        energies_found.append(energy)
        print(f'Reading ~{len(files)} {energy} trees...')
        events, times = get_events(path)
        times, events = zip(*sorted(zip(times, events)))
        event_time_data.update({energy: [list(times), list(events)]})

        total_events.extend(events)
        total_times.extend(times)

    total_times, total_events = zip(*sorted(zip(total_times, total_events)))
    event_time_data.update({'total': [list(total_times), list(total_events)]})

    plot_event_time_data(event_time_data, energies_found)

    print('donzo')


def get_events(path):
    events, times = [], []
    for file_name in os.listdir(path):
        with uproot.open(path) as file:
            events.append(len(file['tree'].arrays('event')['event']))
            times.append(datetime.fromtimestamp(os.path.getmtime(path + file_name)))

    return events, times


if __name__ == '__main__':
    main()
