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

from multiprocessing import Pool
import tqdm
try:
    import istarmap
except ModuleNotFoundError:
    try:
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'Analyzer'))
        import istarmap
    except ModuleNotFoundError:
        print('Can\'t find istarmap!')


def main():
    trees_path = '/gpfs01/star/pwg/dneff/data/AMPT/dylan_run/output/'
    out_path = '/star/u/dneff/'
    energies = [7, 11, 19, 27, 39, 62, '2-7TeV_PbPb']
    threads = 8

    energies_found = []
    event_time_data = {}
    total_events = []
    total_times = []

    print()
    for energy in energies:
        if type(energy) == int:
            energy = f'{energy}GeV'
        path = f'{trees_path}{energy}/'
        if not os.path.exists(path):
            print(f'{path} doesn\'t exist, skipping energy')
            continue
        files = os.listdir(path)
        if len(files) <= 0:
            print(f'{len(files)} {energy} trees, skipping energy')
            continue
        energies_found.append(energy)
        print(f'Reading ~{len(files)} {energy} trees...')
        events, times = get_events(path, threads)
        times, events = zip(*sorted(zip(times, events)))
        event_time_data.update({energy: [list(times), list(events)]})

        total_events.extend(events)
        total_times.extend(times)

    total_times, total_events = zip(*sorted(zip(total_times, total_events)))
    event_time_data.update({'total': [list(total_times), list(total_events)]})

    plot_event_time_data(event_time_data, energies_found, out_path)

    print('donzo')


def get_events(path, threads=1):
    events, times = [], []

    jobs = [(path + file_name,) for file_name in os.listdir(path) if '.root' in file_name]
    with Pool(threads) as pool:
        for num_events, time in tqdm.tqdm(pool.istarmap(get_events_file, jobs), total=len(jobs)):
            events.append(num_events)
            times.append(time)

    return events, times


def get_events_file(path):
    with uproot.open(path) as file:
        num_events = len(file['tree'].arrays('event')['event'])
        time = datetime.fromtimestamp(os.path.getmtime(path))
    return num_events, time


if __name__ == '__main__':
    main()
