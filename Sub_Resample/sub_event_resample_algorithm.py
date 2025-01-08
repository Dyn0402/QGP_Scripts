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
from matplotlib.patches import FancyArrowPatch, Circle
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.lines import Line2D
import math

from DistStats import DistStats


def main():
    """
    Copy and transcribe C++ algorithm to test and visualize
    :return:
    """
    # test_multi_alg4()
    # return
    # animation_function_test()
    # return
    # angles = [0.5, 1.5]
    # angles = np.deg2rad([20, 50, 51, 53, 55, 60, 70, 137, 187, 220, 224, 228, 273, 310, 354])
    angles = np.deg2rad([20, 50, 55, 145, 195, 340])
    bin_width = np.deg2rad(120)  # 2.09
    samples = 180
    samples_list = np.arange(1, 180)
    fps = 5
    # gif_path = '/home/dylan/Research/Results/Presentations/12-21-21/6particles_360samples.gif'
    gif_path = 'E:/Transfer/Research/Resample_POC/Visualizations/animations/alg4_stats.gif'
    # hist = plot_resamples3(angles, bin_width, samples, plot='event')
    # animate_resamples2(angles, bin_width, samples, gif_path, fps)
    # animate_resamples4(angles, bin_width, samples, gif_path, fps)
    # plot_event(angles, 0, bin_width, bin_width, 3)
    # plot_method_paper_event()
    # plot_star_paper_event()
    plot_realistic_event()
    # plot_event_nobin(angles)
    # animate_nsamples_resamples2(angles, bin_width, samples_list, gif_path, fps=fps)
    # animate_nsamples_resamples4(angles, bin_width, samples_list, gif_path, fps=fps)
    # print(hist)
    # plt.hist(hist, bins=np.arange(-0.5, len(angles) + 0.5, 1))
    plt.show()

    print('donzo')


def plot_method_paper_event():
    angles = np.deg2rad([20, 50, 55, 145, 195, 340])
    bin_width = np.deg2rad(120)  # 2.09
    plot_event(angles, 0, bin_width, bin_width, 3)


def plot_star_paper_event():
    angles = np.deg2rad([20, 50, 55, 145, 195, 340])
    bin_width = np.deg2rad(120)  # 2.09
    plot_event(angles, 0, bin_width, bin_width, 3)


def plot_realistic_event():
    # Angles for protons and other particles
    angles = np.deg2rad([20, 50, 55, 145, 195, 340])  # Proton directions
    non_proton_angles = np.deg2rad([80, 110, 160, 250, 290, 310])  # Other particles
    bin_width = np.deg2rad(120)  # Partition width

    # Proton and particle circle radii
    proton_circle_radius = 0.2
    non_proton_circle_radius = 0.18

    # Create the figure
    fig = plt.figure(figsize=(6, 6), dpi=144)
    ax = plt.subplot(111, projection='polar')

    # Background for QGP
    ax.set_facecolor('black')  # Black background
    qgp_background_circle = Circle((0.5, 0.5), 0.45, transform=ax.transAxes, color='white', alpha=1.0, zorder=0)
    qgp_circle = Circle((0.5, 0.5), 0.45, transform=ax.transAxes, color='blue', alpha=0.05, zorder=1)
    ax.add_artist(qgp_background_circle)
    ax.add_artist(qgp_circle)

    # Add protons as circles with arrows
    for angle in angles:
        # Proton circle
        ax.plot(angle, proton_circle_radius, 'o', color='red', alpha=0.6, markersize=6,
                label='Proton' if angle == angles[0] else "", zorder=5)
        # Proton arrow
        arrow = FancyArrowPatch(
            posA=(angle, proton_circle_radius), posB=(angle, 1),
            arrowstyle='-|>', color='red', mutation_scale=20, lw=1.5, zorder=5
        )
        ax.add_artist(arrow)

    # Add other particles as circles with arrows
    for angle in non_proton_angles:
        # Non-proton circle
        ax.plot(angle, non_proton_circle_radius, 'o', color='blue', alpha=0.4, markersize=3,
                label='Other Particle' if angle == non_proton_angles[0] else "", zorder=2)
        # Non-proton arrow
        arrow = FancyArrowPatch(
            posA=(angle, non_proton_circle_radius), posB=(angle, 0.7),
            arrowstyle='-|>', color='blue', alpha=0.3, mutation_scale=15, lw=1, zorder=2
        )
        ax.add_artist(arrow)

    # Collision point at the center
    # ax.plot(0, 0, marker='o', color='yellow', markersize=23, label='Hot Medium', zorder=4)

    # Insert the image at the collision point
    image = mpimg.imread('qgp_medium.png')
    image_box = OffsetImage(image, zoom=0.2)  # Adjust zoom for size
    image_ab = AnnotationBbox(image_box, (0, 0), frameon=False, box_alignment=(0.5, 0.5), zorder=1)
    ax.add_artist(image_ab)

    # Create a custom legend entry for the "Hot Medium"
    image_handle = OffsetImage(image, zoom=0.08)  # Adjust zoom for legend size
    image_handle_box = AnnotationBbox(image_handle, (0, 0), frameon=False, box_alignment=(0.5, 0.5))

    # Partition shading
    bin_low = 0
    bin_high = bin_width
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray',
                    label='120° Partition', zorder=3)

    # Styling and axis settings
    ax.grid(False)
    ax.set_yticklabels([])  # Remove radial labels
    ax.set_ylim((0, 1))  # Set radial limits
    ax.set_theta_zero_location("E")  # 0° at the right

    # Get all handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # Separate handles for tracks and medium/partition
    track_handles = [h for h, l in zip(handles, labels) if 'Proton' in l or 'Particle' in l]
    collision_handles = [h for h, l in zip(handles, labels) if 'Partition' in l or 'Medium' in l]

    # Create the first legend (track)
    legend1 = plt.figlegend(handles=track_handles, loc='upper right',
                            bbox_to_anchor=(0.97, 0.12), fontsize=12, frameon=False)

    # Create the second legend (medium/partition)
    legend2 = plt.figlegend(handles=collision_handles, loc='upper right',
                            bbox_to_anchor=(0.97, 0.99), fontsize=12, frameon=False)

    ax.text(-0.05, 1.03, f'Protons in\nevent: {len(angles)}', transform=ax.transAxes, color='black', fontsize=14,
            ha='left', va='top')
    ax.text(-0.05, 0.05, f'Protons in\npartition: 3', transform=ax.transAxes, color='black', fontsize=14, ha='left',
            va='top')

    fig.tight_layout()


def test_single_alg4():
    """
    For the case of equally spaced bins like in previous algorithm 3, compare algorithm 4 to make sure results are same
    :return:
    """
    angles = np.deg2rad([20, 50, 55, 145, 195, 340])
    bin_width = np.deg2rad(120)  # 2.09
    n_samples = 15510
    dphi = 2 * np.pi / n_samples
    bin_lows = np.arange(0, 2 * np.pi, dphi)
    hist3 = get_resamples3(angles, bin_width, n_samples)
    print(f'hist3: {hist3}')
    print(f'hist3 sum: {np.sum(hist3)}')

    hist4 = get_resamples4_testing(angles, bin_width, n_samples, bin_lows)
    print(f'hist4: {hist4}')
    print(f'hist4 sum: {np.sum(hist4)}')

    print(f'hist3 equal to hist4? {np.all(hist3 == hist4)}')


def test_multi_alg4():
    """
    For the case of equally spaced bins like in previous algorithm 3, compare algorithm 4 to make sure results are same.
    Run many tests with random number of tracks, angles, bin_width, and number of samples. Count number of tests in
    which algorithms produce different results.
    Ran 100,000 tests and no differences found. Confident Algorithm 4 produces identical results for equally spaced bins
    :return:
    """
    n_tests = 100000
    n_bad_tests = 0
    for n_test in range(n_tests):
        n_tracks = np.random.randint(1, 100)
        angles = np.deg2rad(360 * np.random.random(n_tracks))
        bin_width = np.deg2rad(360 * np.random.random())
        n_samples = np.random.randint(1, 10000)

        dphi = 2 * np.pi / n_samples
        bin_lows = np.arange(0, 2 * np.pi, dphi)
        hist3 = get_resamples3(angles, bin_width, n_samples)
        hist4 = get_resamples4_testing(angles, bin_width, n_samples, bin_lows)

        if not np.all(hist3 == hist4):
            print(f'Bad test #{n_test}')
            n_bad_tests += 1

    print('Algorithm 4 testing finished')
    print(f'{n_bad_tests} bad tests')


def get_resamples2(angles_in, bin_width, samples):
    # Same as Algorithm 3 except that single steps of dphi taken each time instead of multiple.
    # Good for illustrating event spaced algorithms without complication of skipping steps in Algorithm 3
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


def get_resamples4(angles_in, bin_width, samples, rng=None):
    # angles = list(angles_in.copy())
    # if bin_width > 2 * np.pi or bin_width <= 0:
    #     print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
    #     bin_width = 2 * np.pi
    # if samples < 0:
    #     print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
    #     samples = abs(samples)

    # sort here

    num_angles = angles_in.size
    hist = np.zeros(num_angles + 1, dtype=int)
    if samples == 0:
        return hist

    # Generate samples random numbers on azimuth
    if rng is None:
        rng = np.random.default_rng()
    bin_lows = np.sort(rng.random(samples)) * 2 * np.pi
    bin_lows = np.append(bin_lows, 6 * np.pi)  # Dummy index for end of algorithm

    # Append duplicate set +2pi to end for bins that wrap around the azimuth
    angles = np.append(angles_in, np.append(angles_in + 2 * np.pi, 4 * np.pi))

    low_index = 0
    high_index = 0
    sample_i = 0
    bin_low = bin_lows[sample_i]
    bin_high = bin_low + bin_width
    while sample_i < samples:
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        sample_i_start = sample_i
        while sample_i < samples and angles[low_index] >= bin_low and angles[high_index] >= bin_high:
            sample_i += 1
            bin_low = bin_lows[sample_i]
            bin_high = bin_low + bin_width
        hist[high_index - low_index] += sample_i - sample_i_start

    return hist


def get_resamples4_testing(angles_in, bin_width, samples, bin_lows):
    # angles = list(angles_in.copy())
    # if bin_width > 2 * np.pi or bin_width <= 0:
    #     print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
    #     bin_width = 2 * np.pi
    # if samples < 0:
    #     print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
    #     samples = abs(samples)

    # sort here

    num_angles = angles_in.size
    hist = np.zeros(num_angles + 1, dtype=int)
    if samples == 0:
        return hist

    # Generate samples random numbers on azimuth
    # if rng is None:
    #     rng = np.random.default_rng()
    # bin_lows = np.sort(rng.rand(samples)) * 2 * np.pi

    # bin_low = 0
    # bin_high = bin_width
    # dphi = 2 * np.pi / samples

    # Append duplicate set +2pi to end for bins that wrap around the azimuth
    angles = np.append(angles_in, np.append(angles_in + 2 * np.pi, 4 * np.pi))
    bin_lows = np.append(bin_lows, 6 * np.pi)  # Dummy index for end of algorithm

    low_index = 0
    high_index = 0
    sample_i = 0
    bin_low = bin_lows[sample_i]
    bin_high = bin_low + bin_width
    while sample_i < samples:
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        sample_i_start = sample_i
        while sample_i < samples and angles[low_index] >= bin_low and angles[high_index] >= bin_high:
            sample_i += 1
            bin_low = bin_lows[sample_i]
            bin_high = bin_low + bin_width
        hist[high_index - low_index] += sample_i - sample_i_start

    return hist


def get_resamples3(angles_in, bin_width, samples):
    # Main evenly spaced algorithm. Calculates number of dphi steps to take before new track added or removed from bin
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


def plot_resamples2(angles_in, bin_width, samples, plot='hist'):
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


def animate_resamples2(angles_in, bin_width, samples, gif_path, fps=10):
    angles = list(angles_in.copy())
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


def animate_resamples4(angles_in, bin_width, samples, gif_path, fps=10, rng=None):
    # angles = angles_in.copy()
    # if bin_width > 2 * np.pi or bin_width <= 0:
    #     print(f'get_resamples bin_width {bin_width} out of range, setting to 2_PI')
    #     bin_width = 2 * np.pi
    # if samples < 0:
    #     print(f'get_resamples samples {samples} less than 0, taking absolute value: {samples} --> {abs(samples)}')
    #     samples = abs(samples)

    num_angles = angles_in.size
    hist = np.zeros(samples, dtype=int)
    if samples == 0:
        return hist

    # Generate samples random numbers on azimuth
    if rng is None:
        rng = np.random.default_rng()
    bin_lows = np.sort(rng.random(samples)) * 2 * np.pi
    bin_lows = np.append(bin_lows, 6 * np.pi)  # Dummy index for end of algorithm

    # Append duplicate set +2pi to end for bins that wrap around the azimuth
    angles = np.append(angles_in, np.append(angles_in + 2 * np.pi, 4 * np.pi))

    low_index = 0
    high_index = 0
    sample_i = 0
    bin_low = bin_lows[sample_i]
    bin_high = bin_low + bin_width
    for sample_i in range(samples):
        while angles[low_index] < bin_low:
            low_index += 1
        while angles[high_index] < bin_high:
            high_index += 1
        hist[sample_i] = high_index - low_index
        bin_low = bin_lows[sample_i]
        bin_high = bin_low + bin_width

    fig = plt.figure(figsize=(10, 5))
    ax = plt.subplot(121, projection='polar')
    ax_hist = plt.subplot(122)

    angles = list(angles)
    # print((samples, angles[:num_angles], bin_width, hist, ax, ax_hist))
    ani = FuncAnimation(fig, ani_func_rand, frames=samples, interval=1.0 / fps * 1000, repeat_delay=5000, repeat=False,
                        fargs=(samples, angles[:num_angles], bin_lows, bin_width, hist, ax, ax_hist))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))

    return hist


def animate_nsamples_resamples2(angles_in, bin_width, samples, gif_path, fps=10):
    angles = list(angles_in.copy())
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

    # hist_stats = []
    # for hist in hists:
    #     stats = DistStats(hist, unbinned=True)
    #     hist_stats.append({'sd': stats.get_sd().val, 'skew': stats.get_skewness().val,
    #                        'kurt': stats.get_kurtosis().val})

    fig = plt.figure(figsize=(10, 5))
    axs = [fig.add_subplot(1, 2, 1),
           fig.add_subplot(3, 2, 2),
           fig.add_subplot(3, 2, 4),
           fig.add_subplot(3, 2, 6)]
    plt.setp(axs[-3].get_xticklabels(), visible=False)
    plt.setp(axs[-2].get_xticklabels(), visible=False)
    axs[-1].set_xlabel('Number of Samples')
    axs[1].get_shared_x_axes().join(axs[1], axs[2])
    axs[2].get_shared_x_axes().join(axs[2], axs[3])
    axs[1].set_xlim((0, max(samples)))
    axs[1].set_ylabel('Standard Deviation')
    axs[2].set_ylabel('Skewness')
    axs[3].set_ylabel('Kurtosis')
    axs[1].grid()
    axs[2].grid()
    axs[3].grid()
    fig.tight_layout()

    # print(samples)
    # print(hists)
    # print(len(samples), len(hists))
    ani = FuncAnimation(fig, ani_nsamples_func, frames=list(zip(samples, hists)), interval=1.0 / fps * 1000,
                        repeat_delay=5000, repeat=False, fargs=(num_angles, axs))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))
    # plt.show()


def animate_nsamples_resamples4(angles_in, bin_width, samples, gif_path, fps=10, rng=None):
    num_angles = angles_in.size

    # Generate samples random numbers on azimuth
    if rng is None:
        rng = np.random.default_rng()

    hists = []
    for sample in samples:
        if sample < 0:
            print(f'get_resamples samples {sample} less than 0, taking absolute value: {sample} --> {abs(sample)}')
            sample = abs(sample)
        hist = np.zeros(sample, dtype=int)
        if sample == 0:
            return hist
        bin_lows = np.sort(rng.random(sample)) * 2 * np.pi
        bin_lows = np.append(bin_lows, 6 * np.pi)  # Dummy index for end of algorithm

        # Append duplicate set +2pi to end for bins that wrap around the azimuth
        angles = np.append(angles_in, np.append(angles_in + 2 * np.pi, 4 * np.pi))

        low_index = 0
        high_index = 0
        sample_i = 0
        bin_low = bin_lows[sample_i]
        bin_high = bin_low + bin_width
        for sample_i in range(sample):
            while angles[low_index] < bin_low:
                low_index += 1
            while angles[high_index] < bin_high:
                high_index += 1
            hist[sample_i] = high_index - low_index
            bin_low = bin_lows[sample_i]
            bin_high = bin_low + bin_width
        hists.append(hist)

    # hist_stats = []
    # for hist in hists:
    #     stats = DistStats(hist, unbinned=True)
    #     hist_stats.append({'sd': stats.get_sd().val, 'skew': stats.get_skewness().val,
    #                        'kurt': stats.get_kurtosis().val})

    fig = plt.figure(figsize=(10, 5))
    axs = [fig.add_subplot(1, 2, 1),
           fig.add_subplot(3, 2, 2),
           fig.add_subplot(3, 2, 4),
           fig.add_subplot(3, 2, 6)]
    plt.setp(axs[-3].get_xticklabels(), visible=False)
    plt.setp(axs[-2].get_xticklabels(), visible=False)
    axs[-1].set_xlabel('Number of Samples')
    axs[1].get_shared_x_axes().join(axs[1], axs[2])
    axs[2].get_shared_x_axes().join(axs[2], axs[3])
    axs[1].set_xlim((0, max(samples)))
    axs[1].set_ylabel('Standard Deviation')
    axs[2].set_ylabel('Skewness')
    axs[3].set_ylabel('Kurtosis')
    axs[1].grid()
    axs[2].grid()
    axs[3].grid()
    fig.tight_layout()

    # print(samples)
    # print(hists)
    # print(len(samples), len(hists))
    ani = FuncAnimation(fig, ani_nsamples_func, frames=list(zip(samples, hists)), interval=1.0 / fps * 1000,
                        repeat_delay=5000, repeat=False, fargs=(num_angles, axs))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))
    # plt.show()


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
    fig = plt.figure(figsize=(4.5, 4.5), dpi=144)
    ax = plt.subplot(111, projection='polar')
    # ax.vlines(angles, 0, 1, ls='--', color='red', label='Tracks')
    ax.plot([], [], color='red', ls='--', label='Protons')  # Just for legend
    for angle in angles:
        arrow = FancyArrowPatch(posA=(angle, 0), posB=(angle, 1), arrowstyle='-|>', color='red', ls='--',
                                mutation_scale=20, shrinkA=0, shrinkB=0)
        ax.add_artist(arrow)
    # ax.arrow(angles[0], 0, angles[0], 1, color='red', width=0.1)
    bw_deg = int(bin_width / np.pi * 180)
    ax.fill_between(np.linspace(bin_low, bin_high, 1000), 0, 1, alpha=0.5, color='gray', label=f'{bw_deg}° Partition')
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(290)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
    ax.text(-0.05, -0.05, f'Protons in \npartition: {counts}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.text(-0.05, 0.97, f'Protons in \nevent: {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    fig.tight_layout()
    # plt.show()


def plot_event_nobin(angles):
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(111, projection='polar')
    ax.plot([], [], color='red', ls='--', label='Tracks')  # Just for legend
    for angle in angles:
        arrow = FancyArrowPatch(posA=(angle, 0), posB=(angle, 1), arrowstyle='-|>', color='red', ls='--',
                                mutation_scale=20, shrinkA=0, shrinkB=0)
        ax.add_artist(arrow)
    ax.grid(False)
    ax.set_yticklabels([])
    ax.set_ylim((0, 1))
    leg_angle = np.deg2rad(300)
    ax.legend(loc="upper left", bbox_to_anchor=(.5 + np.cos(leg_angle) / 2, .5 + np.sin(leg_angle) / 2))
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
    ax.text(-0.12, -0.05, f'Tracks in \nbin: {counts}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.text(-0.12, 1,
            f'Samples:  {int(np.pi * 2 / dphi + 0.5)}\nPhi Step:  {dphi / np.pi * 180:.1f}°\nTracks:     {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    ax_hist.hist(hist, bins=np.arange(-0.5, len(angles) + 1.5), color='red', label='new')
    ax_hist.hist(hist[:-1], bins=np.arange(-0.5, len(angles) + 1.5), color='blue')
    ax_hist.legend()
    ax_hist.set_xlabel('Tracks in Bin')
    plt.tight_layout()


def ani_func_rand(sample_i, n_samples, angles, bin_lows, bin_width, hist, ax, ax_hist):
    print(sample_i)
    bin_low = bin_lows[sample_i]
    bin_high = bin_low + bin_width
    plot_binning_ani_rand(angles, n_samples, bin_low, bin_high, bin_width, hist[sample_i], hist[:sample_i + 1], ax,
                          ax_hist)


def plot_binning_ani_rand(angles, n_samples, bin_low, bin_high, bin_width, counts, hist, ax, ax_hist):
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
    ax.text(-0.12, -0.05, f'Tracks in \nbin: {counts}', horizontalalignment='left', transform=ax.transAxes,
            size='large')
    ax.text(-0.12, 1,
            f'Samples:  {n_samples}\nTracks:     {len(angles)}',
            horizontalalignment='left', transform=ax.transAxes, size='large')
    ax_hist.hist(hist, bins=np.arange(-0.5, len(angles) + 1.5), color='red', label='new')
    ax_hist.hist(hist[:-1], bins=np.arange(-0.5, len(angles) + 1.5), color='blue')
    ax_hist.legend()
    ax_hist.set_xlabel('Tracks in Bin')
    plt.tight_layout()


def ani_nsamples_func(frame_data, num_angles, axs):
    n_sample, hist = frame_data
    print(n_sample)
    plot_binning_nsamples_ani(num_angles, n_sample, hist, axs)


def plot_binning_nsamples_ani(num_angles, n_sample, hist, axs):
    stats = DistStats(hist, unbinned=True)
    sd, skew, kurt = stats.get_sd(), stats.get_skewness(), stats.get_kurtosis()
    # for ax in axs:
    #     ax.clear()
    axs[0].clear()

    axs[0].hist(hist, bins=np.arange(-0.5, num_angles + 1.5), density=True, color='blue')
    axs[0].text(0.65, 0.9, f'Number of \nsamples: {n_sample}', horizontalalignment='left', transform=axs[0].transAxes,
                size='large')

    axs[0].set_xlabel('Tracks in Bin')

    axs[1].scatter(n_sample, sd.val, color='blue', zorder=10)
    axs[2].scatter(n_sample, skew.val, color='blue', zorder=10)
    axs[3].scatter(n_sample, kurt.val, color='blue', zorder=10)

    # plt.setp(axs[-3].get_xticklabels(), visible=False)
    # plt.setp(axs[-2].get_xticklabels(), visible=False)
    # axs[-1].set_xlabel('Number of Samples')
    # axs[1].get_shared_x_axes().join(axs[1], axs[2])
    # axs[2].get_shared_x_axes().join(axs[2], axs[3])
    # axs[1].set_xlim((0, max(samples)))
    # axs[1].set_ylabel('Standard Deviation')
    # axs[2].set_ylabel('Skewness')
    # axs[3].set_ylabel('Kurtosis')
    # axs[1].grid()
    # axs[2].grid()
    # axs[3].grid()
    # fig.tight_layout()

    plt.subplots_adjust(wspace=0.15)


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
