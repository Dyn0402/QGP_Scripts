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
import math

from sub_sample_test import get_hist


def main():
    """
    Copy and transcribe C++ algorithm to test and visualize
    :return:
    """
    angles = [0.5, 1.5]
    # angles = list(np.deg2rad([20, 50, 55, 60, 187, 310]))
    bin_width = np.deg2rad(120)  # 2.09
    samples = 3
    # angles = list(np.deg2rad(angles))
    print(get_resamples(np.asarray(angles), bin_width, samples))
    print(get_resamples3(angles, bin_width, samples))
    # angles = list(np.rad2deg(angles))
    # print(get_hist(angles, np.rad2deg(bin_width), samples))
    hist = plot_resamples3(angles, bin_width, samples)
    # print(hist)
    # plt.hist(hist, bins=np.arange(-0.5, len(angles) + 0.5, 1))
    # plt.show()

    print('donzo')


def get_resamples3(angles_in, bin_width, samples):
    angles = angles_in.copy()
    if bin_width > 2 * np.pi or bin_width <= 0:
        print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
        bin_width = 2 * np.pi
    if samples < 0:
        print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
        samples = abs(samples)

    # sort here

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


def plot_resamples3(angles_in, bin_width, samples):
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
    # print(angles)
    for sample_i in range(samples):
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        hist[sample_i] = high_index - low_index
        plot_binning(angles[:num_angles], bin_low, bin_high, hist[sample_i], hist[:sample_i + 1])
        bin_low += dphi
        bin_high += dphi

    return hist


def plot_binning(angles, bin_low, bin_high, counts, hist):
    fig = plt.figure(figsize=(10, 5))
    ax = plt.subplot(121, projection='polar')
    ax_hist = plt.subplot(122)
    # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.vlines(angles, 0, 1, color='red', label='tracks')
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray', label='bin')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    ax.text(-0.1, -0.02, f'Tracks in \nbin: {counts}', horizontalalignment='left', transform=ax.transAxes)
    # print(hist)
    print(hist, counts)
    ax_hist.hist(hist, bins=np.arange(-0.5, len(angles) + 1.5), color='red', label='new')
    ax_hist.hist(hist[:-1], bins=np.arange(-0.5, len(angles) + 1.5), color='blue')
    ax_hist.legend()
    ax_hist.set_xlabel('Tracks in Bin')
    fig.tight_layout()
    plt.show()


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


if __name__ == '__main__':
    main()
