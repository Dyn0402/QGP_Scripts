#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 27 11:18 AM 2023
Created in PyCharm
Created as QGP_Scripts/momentum_conservation_model

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.optimize import curve_fit as cf

from multiprocessing import Pool
import tqdm
import istarmap  # Needed for tqdm

import uproot
import awkward as ak
import vector

from Analysis_POCs.poc_functions import bin_experiment, bin_experiment_no_bs
from PConsSim import PConsSim, rotate_vector
from DistStats import DistStats
from Measure import Measure


def main():
    # chat_gpt_test()
    # single_event_test()
    # momentum_test()
    # full_test()
    # rotation_test()
    plot_full_test_from_file()
    print('donzo')


def full_test():
    ss = np.random.SeedSequence(42)

    n_tracks = [50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 240, 280, 320, 390, 500, 640, 800, 1000]
    n_protons_frac = 0.4
    n_events = 300000
    energy = 2  # Currently just the range of momenta
    y_max, pt_max, p_max = 0.5, 2, 2

    bin_width = np.deg2rad(120)
    resamples = 72

    threads = 16

    n_rndms = (n_events + int(max(n_tracks) * n_protons_frac)) * len(n_tracks)
    rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_rndms)])

    dsig_avgs, dsig_avg_errs = [], []
    for n_i, n_track in enumerate(n_tracks):
        print(f'\nStarting {n_track} track events {n_i + 1}/{len(n_tracks)}:')
        n_protons = int(n_track * n_protons_frac + 0.5)
        jobs = [(n_track, next(rngs), energy, n_protons, y_max, pt_max, p_max) for i in range(n_events)]
        experiment_tracks = {}
        with Pool(threads) as pool:
            for tracks in tqdm.tqdm(pool.istarmap(sim_event, jobs), total=len(jobs)):
                if len(tracks) not in experiment_tracks:
                    experiment_tracks.update({len(tracks): [tracks]})
                else:
                    experiment_tracks[len(tracks)].append(tracks)
        # for event_i, rng in zip(range(n_events), rngs):
            # if event_i % 100 == 0:
            #     print(f'Event #{event_i} n_track={n_track} {time.time() - start_time:.2f}s')
            #
            # event = PConsSim(energy, n_track, rng=rng)
            # event.rotate_tracks()
            #
            # while event.init_net_momentum <= event.net_momentum:
            #     angle_frac = event.angle_frac
            #     print(f' Event #{event_i}: net_p {event.init_net_momentum:.1f}->{event.net_momentum:.1f}')
            #     # print(f'Old init_momentum={event.init_net_momentum}, old final_momentum={event.net_momentum}')
            #     event = PConsSim(energy, n_track, rng=rng)
            #     event.angle_frac = angle_frac / 2
            #     # print(f'Old angle_frac={angle_frac}, new angle_frac={event.angle_frac}')
            #     event.rotate_tracks()
            #
            # tracks = event.get_m_tracks(n_protons)
            #
            # tracks = ak.Array({'px': tracks[:, 0], 'py': tracks[:, 1], 'pz': tracks[:, 2], 'M': tracks[:, 0] * 0},
            #                   with_name='Momentum4D')
            # # print(tracks.p)
            # tracks = tracks[(abs(tracks.rapidity) < y_max) & (tracks.pt < pt_max) & (tracks.p < p_max)]
            # # print(len(tracks))
            # tracks = np.array(tracks.phi)
            # tracks[tracks < 0] += 2 * np.pi
            # # if len(tracks) == 27:
            # #     plot_momenta(event.get_m_tracks(n_track), event.net_momentum_vec)
            # #     plot_momenta(event.get_m_tracks(n_track), event.net_momentum_vec)
            # # tracks = np.array(np.random.uniform(0, 2 * np.pi, size=len(tracks)))
            # if len(tracks) not in experiment_tracks:
            #     experiment_tracks.update({len(tracks): [tracks]})
            # else:
            #     experiment_tracks[len(tracks)].append(tracks)

        n_protons_list = sorted([x for x in experiment_tracks.keys() if x > 1])
        dsig_2_list, dsig_2_err_list = [], []
        for n_proton in n_protons_list:
            p = bin_width / (2 * np.pi)
            binom_var = n_proton * p * (1 - p)
            data = bin_experiment_no_bs(experiment_tracks[n_proton], n_proton, bin_width, resamples, next(rngs),
                                        alg=3, sort=True)
            stats = DistStats(data, unbinned=False)
            delta_sig_2 = (stats.get_k_stat(2) - binom_var) / (n_proton * (n_proton - 1))
            dsig_2_list.append(delta_sig_2.val)
            dsig_2_err_list.append(delta_sig_2.err * np.sqrt(resamples))
            # print(
            #     f'{n_proton} protons: mean: {stats.get_mean()}, variance: {stats.get_variance()}, delta_sig_2: {delta_sig_2}')

        print(f'Total Particles {n_track}')
        print(n_protons_list)
        print(dsig_2_list)
        print(dsig_2_err_list)
        wavg, wavg_err = plot_vs_total_protons(n_protons_list, dsig_2_list, dsig_2_err_list, n_track=n_track)
        dsig_avgs.append(wavg)
        dsig_avg_errs.append(wavg_err)

    print(r'$\langle \Delta \sigma^2 \rangle$')
    print(n_tracks)
    print(dsig_avgs)
    print(dsig_avg_errs)
    plot_vs_total_particles(n_tracks, dsig_avgs, dsig_avg_errs)
    plt.show()


def plot_full_test_from_file():
    # path = 'N:/UCLA_Microsoft/OneDrive - personalmicrosoftsoftware.ucla.edu/Research/UCLA/Results/mom_cons_new_pc.txt'
    path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Results/mom_cons_new_pc.txt'
    total_tracks, total_protons, dsigs, dsig_errs = [], [], [], []
    with open(path, 'r') as file:
        lines = file.readlines()
    rolling = False
    for line in lines:
        if rolling == 1:
            total_protons.append([int(x) for x in line.strip('\n[]').split(', ')])
            rolling += 1
        elif rolling == 2:
            dsigs.append([float(x) for x in line.strip('\n[]').split(', ')])
            rolling += 1
        elif rolling == 3:
            dsig_errs.append([float(x) for x in line.strip('\n[]').split(', ')])
            rolling = 0
        if 'Total Particles ' in line:
            total_tracks.append(int(line.strip().replace('Total Particles ', '')))
            rolling = 1

    # print(total_tracks)
    # print(total_protons)
    # print(dsigs)
    # print(dsig_errs)

    dsig_avgs, dsig_avg_errs = [], []
    for tracks, protons, dsigs_i, dsig_errs_i in zip(total_tracks, total_protons, dsigs, dsig_errs):
        wavg, wavg_err = plot_vs_total_protons(protons, dsigs_i, dsig_errs_i, tracks, plot=False)
        dsig_avgs.append(wavg)
        dsig_avg_errs.append(wavg_err)
    plot_vs_total_particles(total_tracks, dsig_avgs, dsig_avg_errs)
    plot_fits(total_tracks, dsig_avgs, dsig_avg_errs)

    plt.show()


def plot_fits(total_tracks, dsig_avgs, dsig_avg_errs):
    fits = {
        # r'$\frac{a}{\sqrt{N}} + b$': sqrt_n_fit,
        # r'$\frac{a}{N} + b$': n_fit,
        r'$\frac{a}{N^c} + b$': n_exp_fit,
    }
    xs = np.linspace(min(total_tracks), max(total_tracks), 1000)
    par_letters = ['a', 'b', 'c']
    for fit_lab, fit in fits.items():
        popt, pcov = cf(fit, total_tracks, dsig_avgs, sigma=dsig_avg_errs, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        fit_values = {par: Measure(val, err) for par, val, err in zip(par_letters, popt, perr)}

        plt.plot(xs, fit(xs, *popt), label=fit_lab)

        # Create the textbox
        x_pos, y_pos = 400, -0.001  # Set the position of the textbox
        textbox_content = '\n'.join([f"{key} = {value}" for key, value in fit_values.items()])
        textbox = plt.text(x_pos, y_pos, textbox_content, fontsize=12, va='center', ha='left',
                           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
        print('\n', fit_lab)
        print(', '.join([f'{par} = {meas}' for par, meas in fit_values.items()]))
    plt.legend(loc='lower right', fontsize=20)
    plt.tight_layout()


def single_event_test():
    n_part = 50
    seed = 47

    bad_event = False
    while not bad_event:
        seed += 1
        rng = np.random.default_rng(seed=seed)
        event = PConsSim(5, n_part, rng=rng)
        event.rotate_tracks()
        p_net_0, p_net_f = event.init_net_momentum, event.net_momentum_iterations[-1]
        print(f'Seed = {seed}: p_net_0 = {p_net_0}, p_net_f = {p_net_f}')
        print(f'{event.net_momentum_iterations}')
        if event.net_momentum_iterations[-1] > event.init_net_momentum:
            bad_event = True
    angle_frac_init = event.angle_frac
    fig, ax = plt.subplots()
    ax.grid()
    rng = np.random.default_rng(seed=seed)
    for x in range(1, 8):
        event = PConsSim(5, n_part, rng=rng)
        event.angle_frac = angle_frac_init / x
        # print(event.get_m_tracks(20))
        # plot_momenta(event.get_m_tracks(n_part), event.net_momentum_vec)
        # event.rotate_tracks_debug()
        event.rotate_tracks()
        # print(event.get_m_tracks(20))
        # plot_momenta(event.get_m_tracks(n_part), event.net_momentum_vec)
        ax.plot(range(len(event.net_momentum_iterations)), event.net_momentum_iterations,
                label=f'{event.angle_frac:.4f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Event Net Momentum')
    plt.legend()
    plt.tight_layout()
    plt.show()


def momentum_test():
    n_particles = 50
    angle_frac = 0.8 / np.sqrt(n_particles)
    iterations = 4
    momenta = np.random.uniform(-5, 5, size=(n_particles, 3))

    net_momentum_vec = np.sum(momenta, axis=0)
    init_net_momentum = np.linalg.norm(net_momentum_vec)
    print(f'Initial: Net-Momentum={init_net_momentum}')
    plot_momenta(momenta, net_momentum_vec)
    net_momentum = init_net_momentum

    net_momentum_list = [net_momentum]
    for i in range(iterations):
        for j in range(len(momenta)):
            momenta[j] = rotate_vector(momenta[j], -net_momentum_vec, angle_frac * net_momentum / init_net_momentum)

        net_momentum_vec = np.sum(momenta, axis=0)
        net_momentum = np.linalg.norm(net_momentum_vec)
        print(f'Iteration {i}: Net-Momentum={net_momentum}')
        plot_momenta(momenta, net_momentum_vec)
        net_momentum_list.append(net_momentum)

    fig, ax = plt.subplots(dpi=144)
    ax.plot(range(len(net_momentum_list)), net_momentum_list)
    ax.axhline(0, color='black')
    ax.set_xlabel('Number of Iterations')
    ax.set_ylabel('Net Momentum Magnitude')
    fig.tight_layout()
    plt.show()


def rotation_test():
    v1 = np.array([1, 0, 0])  # Initial unit vector
    v2 = np.array([0, 1, 0])  # Target unit vector

    rotated_v1 = rotate_vector(v1, v2, angle_fraction=0.3)
    print("Rotated v1:", rotated_v1)


def chat_gpt_test():
    # Define the vectors
    # vec1 = np.array([1, 1, 0.5])  # Initial unit vector to be rotated
    vec1 = np.array([-1.12300128, 0.9577349, 0.13800545])  # Initial unit vector to be rotated
    # vec1 /= np.linalg.norm(vec1)
    # vec2 = np.array([0.1, 1, 0.2])  # Target unit vector
    vec2 = np.array([-1.89403927, 8.96250786, 14.5447655])  # Target unit vector
    # vec2 /= np.linalg.norm(vec2)

    rotated_vec1 = rotate_vector(vec1, vec2, angle_fraction=0.2)

    # Create a 3D figure and add axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set limits for the axes
    ax.set_xlim([0, np.linalg.norm(vec2)])
    ax.set_ylim([0, np.linalg.norm(vec2)])
    ax.set_zlim([0, np.linalg.norm(vec2)])

    # Set the origin point
    ax.quiver(0, 0, 0, 0, 0, 0, color='black')

    # Plot the vectors
    ax.quiver(0, 0, 0, *(vec1 / np.linalg.norm(vec1)), color='red', length=np.linalg.norm(vec1), label='vec1')
    ax.quiver(0, 0, 0, *(vec2 / np.linalg.norm(vec2)), color='blue', length=np.linalg.norm(vec2), label='vec2')
    ax.quiver(0, 0, 0, *(rotated_vec1 / np.linalg.norm(rotated_vec1)), color='green',
              length=np.linalg.norm(rotated_vec1), label='rotated vec1')
    # ax.quiver(0, 0, 0, vec1[0], vec1[1], vec1[2], color='red', length=np.linalg.norm(vec1), label='vec1')
    # ax.quiver(0, 0, 0, vec2[0], vec2[1], vec2[2], color='blue', length=np.linalg.norm(vec2), label='vec2')
    # ax.quiver(0, 0, 0, rotated_vec1[0], rotated_vec1[1], rotated_vec1[2], color='green',
    #           length=np.linalg.norm(rotated_vec1), label='rotated vec1')
    print(vec1, np.linalg.norm(vec1))
    print(vec2)
    print(rotated_vec1, np.linalg.norm(rotated_vec1))

    # Add labels for the vectors
    ax.text(vec1[0], vec1[1], vec1[2], "vec1", color='red')
    ax.text(vec2[0], vec2[1], vec2[2], "vec2", color='blue')

    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.legend()

    # Show the plot
    plt.show()


def sim_event(n_track, rng, energy, n_protons, y_max, pt_max, p_max):
    vector.register_awkward()
    event = PConsSim(energy, n_track, rng=rng)
    event.rotate_tracks()

    while event.init_net_momentum <= event.net_momentum:
        angle_frac = event.angle_frac
        event = PConsSim(energy, n_track, rng=rng)
        event.angle_frac = angle_frac / 2
        event.rotate_tracks()

    tracks = event.get_m_tracks(n_protons)

    tracks = ak.Array({'px': tracks[:, 0], 'py': tracks[:, 1], 'pz': tracks[:, 2], 'M': tracks[:, 0] * 0},
                      with_name='Momentum4D')
    tracks = tracks[(abs(tracks.rapidity) < y_max) & (tracks.pt < pt_max) & (tracks.p < p_max)]
    tracks = np.array(tracks.phi)
    tracks[tracks < 0] += 2 * np.pi

    return tracks


def plot_vs_total_protons(n_protons, dsig2, dsig2_err, n_track='-', plot=True):
    dsig_2_list, dsig_2_err_list = np.array(dsig2), np.array(dsig2_err)
    weight_avg = np.average(dsig_2_list, weights=1 / dsig_2_err_list ** 2)
    weight_avg_err = np.sqrt(1 / np.sum(1 / dsig_2_err_list ** 2))
    if plot:
        plt.figure()
        plt.axhline(0, color='black')
        plt.grid()
        plt.title(f'Total Particles {n_track}')

        plt.errorbar(n_protons, dsig_2_list, yerr=dsig_2_err_list, marker='o', ls='none')
        plt.xlabel('Total Protons Within Acceptance')
        plt.ylabel(r'$\Delta\sigma^2$')
        plt.axhline(weight_avg)
        plt.tight_layout()

    return weight_avg, weight_avg_err


def plot_vs_total_particles(n_tracks, dsig_avgs, dsig_avg_errs):
    plt.figure()
    plt.grid()
    plt.axhline(0, color='black')
    plt.errorbar(n_tracks, dsig_avgs, yerr=dsig_avg_errs, ls='none', marker='o', alpha=0.8)
    plt.xlabel('Total Particles in Event')
    plt.ylabel(r'$\langle\Delta\sigma^2\rangle$')
    plt.tight_layout()


def plot_vectors(vectors):
    # Create a 3D figure and add axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set the origin point
    ax.quiver(0, 0, 0, 0, 0, 0, color='black')

    # Plot the vectors
    max_axis = 0
    for vec in vectors:
        vec_len = np.linalg.norm(vec)
        ax.quiver(0, 0, 0, *(vec / vec_len), length=vec_len)
        for vec_i in vec:
            max_axis = vec_i if vec_i > max_axis else max_axis

    ax.set_xlim([-max_axis, max_axis])
    ax.set_ylim([-max_axis, max_axis])
    ax.set_zlim([-max_axis, max_axis])

    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def plot_momenta(momenta, net_momentum_vec):
    plot_vectors(momenta)
    net_momentum = np.linalg.norm(net_momentum_vec)
    plt.quiver(0, 0, 0, *(net_momentum_vec / net_momentum), length=net_momentum, color='red', label='Net Momentum')
    plt.quiver(0, 0, 0, *(-net_momentum_vec / net_momentum), length=net_momentum, color='orange', label='-Net Momentum')
    plt.legend()
    plt.tight_layout()


def rapidity(px, py, pz, m):
    e = np.sqrt(m ** 2 + px ** 2 + py ** 2 + pz ** 2)
    return np.arctanh(pz / e)


def sqrt_n_fit(n, a, b):
    return a / np.sqrt(n) + b


def n_fit(n, a, b):
    return a / n + b


def n_exp_fit(n, a, b, c):
    return a / n**c + b


if __name__ == '__main__':
    main()
