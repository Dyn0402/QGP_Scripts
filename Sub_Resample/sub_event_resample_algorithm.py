#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 03 1:14 PM 2021
Created in PyCharm
Created as QGP_Scripts/sub_event_resample_algorithm.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import math

from sub_sample_test import get_hist


def main():
    """
    Copy and transcribe C++ algorithm to test and visualize
    :return:
    """
    # animation_function_test()
    # return
    # angles = [0.5, 1.5]
    # angles = list(np.deg2rad([20, 50, 51, 53, 55, 60, 70, 137, 187, 220, 224, 228, 273, 310, 354]))
    angles = list(np.deg2rad([20, 50, 55, 145, 195, 340]))
    bin_width = np.deg2rad(120)  # 2.09
    samples = 360
    samples_list = np.arange(1, 361)
    fps = 10
    # gif_path = '/home/dylan/Research/Results/Presentations/12-21-21/6particles_360samples.gif'
    gif_path = 'D:/Research/Results/Resample_POC/ani_resamp_test.gif'
    # angles = list(np.deg2rad(angles))
    # print(get_resamples(np.asarray(angles), bin_width, samples))
    # print(get_resamples3(angles, bin_width, samples))
    # angles = list(np.rad2deg(angles))
    # print(get_hist(angles, np.rad2deg(bin_width), samples))
    # hist = plot_resamples3(angles, bin_width, samples, plot='event')
    # animate_resamples3(angles, bin_width, samples, gif_path, fps)
    animate_nsamples_resamples3(angles, bin_width, samples_list, gif_path, fps=10)
    # print(hist)
    # plt.hist(hist, bins=np.arange(-0.5, len(angles) + 0.5, 1))
    # plt.show()

    print('donzo')


def get_resamples3(angles_in, bin_width, samples):
    angles = list(angles_in.copy())
    if bin_width > 2 * np.pi or bin_width <= 0:
        print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
        bin_width = 2 * np.pi
    if samples < 0:
        print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
        samples = abs(samples)

    # sort here

    hist = np.zeros(samples, dtype=int)  # Changed from "empty" to "zeros", think this is better?
    if samples == 0:
        return hist
    bin_low = 0
    bin_high = bin_width
    dphi = 2 * np.pi / samples

    num_angles = len(angles)
    for i in range(num_angles):
        if angles[i] >= 2 * np.pi or angles[i] < 0:
            print('bad angle range')
        if angles[i + 1] < angles[i]:
            print('angles unsorted')
        angles.append(angles[i] + 2 * np.pi)
    angles.append(4 * np.pi)  # Last goalpost

    low_index = 0
    high_index = 0
    for sample_i in range(samples):
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        hist[sample_i] = high_index - low_index
        bin_low += dphi
        bin_high += dphi

    return hist


def get_resamples(angles_in, bin_width, samples):
    # Expecting angles is numpy array
    # angles = angles_in.copy()
    # if bin_width > 2 * np.pi or bin_width <= 0:
    #     print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
    #     bin_width = 2 * np.pi
    # if samples < 0:
    #     print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
    #     samples = abs(samples)

    # sort here

    bin_low = 0
    bin_high = bin_width
    dphi = 2 * np.pi / samples

    num_angles = angles_in.size
    hist = np.zeros(num_angles + 1, dtype=int)
    if samples == 0:
        return hist

    # for i in range(num_angles):
        # if angles[i] >= 2 * np.pi or angles[i] < 0:
        #     print('bad angle range')
        # if angles[i + 1] < angles[i]:
        #     print('angles unsorted')
    angles = np.append(angles_in, np.append(angles_in + 2 * np.pi, 4 * np.pi))

    low_index = 0
    high_index = 0
    sample_i = 0
    while sample_i < samples:
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        step = math.ceil((angles[low_index] - bin_low) / dphi)
        step_high = math.ceil((angles[high_index] - bin_high) / dphi)
        if step_high < step:
            step = step_high
        hist[high_index - low_index] += step
        sample_i += step
        bin_low += step * dphi
        bin_high += step * dphi

    hist[high_index - low_index] -= sample_i - samples

    return hist


def plot_resamples3(angles_in, bin_width, samples, plot='hist'):
    angles = angles_in.copy()
    if bin_width > 2 * np.pi or bin_width <= 0:
        print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
        bin_width = 2 * np.pi
    if samples < 0:
        print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
        samples = abs(samples)

    hist = np.empty(samples, dtype=int)
    if samples == 0:
        return hist
    bin_low = 0
    bin_high = bin_width
    dphi = 2 * np.pi / samples

    num_angles = len(angles)
    for i in range(num_angles):
        if angles[i] >= 2 * np.pi or angles[i] < 0:
            print('bad angle range')
        if angles[i + 1] < angles[i]:
            print('angles unsorted')
        angles.append(angles[i] + 2 * np.pi)
    angles.append(4 * np.pi)  # Last goalpost

    low_index = 0
    high_index = 0
    for sample_i in range(samples):
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        hist[sample_i] = high_index - low_index
        if plot == 'hist':
            plot_binning(angles[:num_angles], bin_low, bin_high, dphi, bin_width, hist[sample_i], hist[:sample_i + 1])
        elif plot == 'event':
            plot_event(angles[:num_angles], bin_low, bin_high, bin_width, hist[sample_i])
        bin_low += dphi
        bin_high += dphi

    return hist


def animate_resamples3(angles_in, bin_width, samples, gif_path, fps=10):
    angles = angles_in.copy()
    if bin_width > 2 * np.pi or bin_width <= 0:
        print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
        bin_width = 2 * np.pi
    if samples < 0:
        print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
        samples = abs(samples)

    hist = np.empty(samples, dtype=int)
    if samples == 0:
        return hist
    bin_low = 0
    bin_high = bin_width
    dphi = 2 * np.pi / samples

    num_angles = len(angles)
    for i in range(num_angles):
        if angles[i] >= 2 * np.pi or angles[i] < 0:
            print('bad angle range')
        if angles[i + 1] < angles[i]:
            print('angles unsorted')
        angles.append(angles[i] + 2 * np.pi)
    angles.append(4 * np.pi)  # Last goalpost

    low_index = 0
    high_index = 0
    for sample_i in range(samples):
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        hist[sample_i] = high_index - low_index
        bin_low += dphi
        bin_high += dphi

    fig = plt.figure(figsize=(10, 5))
    ax = plt.subplot(121, projection='polar')
    ax_hist = plt.subplot(122)

    ani = FuncAnimation(fig, ani_func, frames=samples, interval=1.0 / fps * 1000, repeat_delay=5000, repeat=False,
                        fargs=(angles[:num_angles], dphi, bin_width, hist, ax, ax_hist))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))

    return hist


def animate_nsamples_resamples3(angles_in, bin_width, samples, gif_path, fps=10):
    angles = angles_in.copy()
    if bin_width > 2 * np.pi or bin_width <= 0:
        print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
        bin_width = 2 * np.pi

    num_angles = len(angles)
    hists = []
    for sample in samples:
        if sample < 0:
            print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
            samples = abs(samples)
        hist = np.empty(sample, dtype=int)
        if sample == 0:
            return hist
        bin_low = 0
        bin_high = bin_width
        dphi = 2 * np.pi / sample

        for i in range(num_angles):
            if angles[i] >= 2 * np.pi or angles[i] < 0:
                print('bad angle range')
            if angles[i + 1] < angles[i]:
                print('angles unsorted')
            angles.append(angles[i] + 2 * np.pi)
        angles.append(4 * np.pi)  # Last goalpost

        low_index = 0
        high_index = 0
        for sample_i in range(sample):
            while angles[low_index] < bin_low:
                low_index += 1
            while angles[high_index] < bin_high:
                high_index += 1
            hist[sample_i] = high_index - low_index
            bin_low += dphi
            bin_high += dphi
        hists.append(hist)

    fig, ax = plt.subplots(figsize=(10, 5))

    print(samples)
    print(hists)
    print(len(samples), len(hists))
    ani = FuncAnimation(fig, ani_nsamples_func, frames=list(zip(samples, hists)), interval=1.0 / fps * 1000,
                        repeat_delay=5000, repeat=False, fargs=(num_angles, ax))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))


def plot_binning(angles, bin_low, bin_high, dphi, bin_width, counts, hist):
    fig = plt.figure(figsize=(10, 5))
    ax = plt.subplot(121, projection='polar')
    ax_hist = plt.subplot(122)
    ax.vlines(angles, 0, 1, color='red', label='tracks')
    bw_deg = int(bin_width / np.pi * 180)
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray', label=f'{bw_deg}° bin')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    ax.text(-0.12, -0.05, f'Tracks in \nbin:   {counts}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.text(-0.12, 1,
            f'Samples:  {int(np.pi * 2 / dphi + 0.5)}\nPhi Step:  {dphi / np.pi * 180:.1f}°\nTracks:     {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    ax_hist.hist(hist, bins=np.arange(-0.5, len(angles) + 1.5), color='red', label='new')
    ax_hist.hist(hist[:-1], bins=np.arange(-0.5, len(angles) + 1.5), color='blue')
    ax_hist.legend()
    ax_hist.set_xlabel('Tracks in Bin')
    fig.tight_layout()
    plt.show()


def plot_event(angles, bin_low, bin_high, bin_width, counts):
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(111, projection='polar')
    ax.vlines(angles, 0, 1, color='red', label='tracks')
    bw_deg = int(bin_width / np.pi * 180)
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray', label=f'{bw_deg}° bin')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    ax.text(-0.05, -0.05, f'Tracks in \nbin:  {counts}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.text(-0.05, 0.97, f'Tracks in \nevent: {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    fig.tight_layout()
    plt.show()


def ani_func(sample_i, angles, dphi, bin_width, hist, ax, ax_hist):
    print(sample_i)
    bin_low = dphi * sample_i
    bin_high = bin_low + bin_width
    plot_binning_ani(angles, bin_low, bin_high, bin_width, dphi, hist[sample_i], hist[:sample_i + 1], ax, ax_hist)


def plot_binning_ani(angles, bin_low, bin_high, bin_width, dphi, counts, hist, ax, ax_hist):
    ax.clear()
    ax_hist.clear()
    ax.vlines(angles, 0, 1, color='red', label='tracks')
    bw_deg = int(bin_width / np.pi * 180)
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray', label=f'{bw_deg}° bin')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    ax.text(-0.12, -0.05, f'Tracks in \nbin: {counts}', horizontalalignment='left', transform=ax.transAxes, size='large')
    ax.text(-0.12, 1,
            f'Samples:  {int(np.pi * 2 / dphi + 0.5)}\nPhi Step:  {dphi / np.pi * 180:.1f}°\nTracks:     {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    ax_hist.hist(hist, bins=np.arange(-0.5, len(angles) + 1.5), color='red', label='new')
    ax_hist.hist(hist[:-1], bins=np.arange(-0.5, len(angles) + 1.5), color='blue')
    ax_hist.legend()
    ax_hist.set_xlabel('Tracks in Bin')
    plt.tight_layout()


def ani_nsamples_func(frame_data, num_angles, ax):
    n_sample, hist = frame_data
    print(n_sample)
    plot_binning_nsamples_ani(num_angles, n_sample, hist, ax)


def plot_binning_nsamples_ani(num_angles, n_sample, hist, ax):
    ax.clear()
    ax.hist(hist, bins=np.arange(-0.5, num_angles + 1.5), density=True, color='blue')
    ax.text(0.75, 0.85, f'Number of samples: {n_sample}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.set_xlabel('Tracks in Bin')
    plt.tight_layout()


def polar_plot_test():
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # ax.set_theta_zero_location()
    # fig = plt.figure(figsize=(6, 6))
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar')
    ax.vlines([0.5 * np.pi, 0.6 * np.pi], 0, 1, color='green', label='tracks')
    ax.fill_between(np.linspace(0, 0.4 * np.pi, 1000), 0, 1, label='bin')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    plt.show()


def animation_function_test():
    fig, ax = plt.subplots()
    ims = []
    phases = np.linspace(0, 2, 120)
    x = np.linspace(0, 10, 100)
    # y = np.cos(x + phases[0])
    # plt.plot(x, y)
    # plt.show()
    # for phase in phases:
    #     ims.append(plt.plot(x, np.cos(x + phase)))

    fps = 60
    interval = 1.0 / fps * 1000  # ms / frame
    print(interval * len(phases))
    print(len(phases) / fps)
    ani = FuncAnimation(fig, ani_func_test, frames=phases, interval=interval * 8, fargs=(x, ax))
    ani.save('/home/dylan/Desktop/gif_test.gif', dpi=100, writer=PillowWriter(fps=fps / 2))


def ani_func_test(phases, x, ax):
    ax.clear()
    plt.plot(x, np.cos(x + phases))


if __name__ == '__main__':
    main()
