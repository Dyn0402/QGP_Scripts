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
import math
import os
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
from analyze_binom_slices import vn_divs


def main():
    # chat_gpt_test()
    # single_event_test()
    # check_convergence_vs_alpha()
    # test_alpha_adaptive()
    # momentum_test()
    # full_test()
    # rotation_test()
    plot_full_test_from_file()
    # plot_from_files()
    plt.show()
    print('donzo')


def full_test():
    out_path = 'C:/Users/Dylan/Desktop/test/'
    # out_path = '/home/dylan/mom_cons_model/'
    ss = np.random.SeedSequence(42)

    n_tracks = [30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 240, 280, 320, 390, 500, 640, 800]
    # n_tracks = [30, 40, 50]
    n_tracks.reverse()
    # n_protons_fracs = [0.4, 0.6]
    n_protons_fracs = [0.4]
    n_events = 100000
    # convergence_momenta = [0.001, 1]
    convergence_momenta = [0.001]
    # energies = [2, 4]  # Currently just the range of momenta
    energies = [2]  # Currently just the range of momenta
    y_max, pt_max, p_max = 0.5, 2, 2

    bin_widths = np.deg2rad([60, 72, 90, 120, 180, 288])
    resamples = 72

    threads = 15

    variations = len(n_protons_fracs) * len(energies) * len(convergence_momenta)

    variation_num = 1
    for n_protons_frac in n_protons_fracs:
        for convergence_momentum in convergence_momenta:
            for energy in energies:
                var_dir = (f'{out_path}pfrac{int(n_protons_frac * 100 + 0.1)}_'
                           f'convp{int(convergence_momentum * 1000 + 0.1)}_energy{energy}_nevent{n_events/1000}k/')
                if os.path.exists(var_dir):
                    os.rmdir(var_dir)
                os.mkdir(var_dir)
                print(f'\n\nStarting variation {variation_num}/{variations} ({var_dir.split("/")[-2]}):\n')
                for n_i, n_track in enumerate(n_tracks):
                    print(f'Starting {n_track} track events {n_i + 1}/{len(n_tracks)}:')
                    n_protons = int(n_track * n_protons_frac + 0.5)
                    # n_rndms = n_events + int(n_track * n_protons_frac * len(bin_widths))
                    # rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_rndms)])
                    rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_events)])
                    jobs = [(n_track, next(rngs), energy, n_protons, y_max, pt_max, p_max, convergence_momentum)
                            for i in range(n_events)]
                    experiment_tracks = {}
                    with Pool(threads) as pool:
                        for tracks in tqdm.tqdm(pool.istarmap(sim_event, jobs), total=len(jobs)):
                            if len(tracks) not in experiment_tracks:
                                experiment_tracks.update({len(tracks): [tracks]})
                            else:
                                experiment_tracks[len(tracks)].append(tracks)

                    n_protons_list = sorted([x for x in experiment_tracks.keys() if x > 1])
                    n_rndms = len(n_protons_list) * len(bin_widths)
                    rngs = iter([np.random.default_rng(s) for s in ss.spawn(n_rndms)])
                    for bin_width in bin_widths:
                        dsig_2_list, dsig_2_err_list, v1_list, v2_list, v3_list = [], [], [], [], []
                        for n_proton in n_protons_list:
                            p = bin_width / (2 * np.pi)
                            binom_var = n_proton * p * (1 - p)
                            data = bin_experiment_no_bs(experiment_tracks[n_proton], n_proton, bin_width, resamples,
                                                        next(rngs), alg=3, sort=True)
                            stats = DistStats(data, unbinned=False)
                            delta_sig_2 = (stats.get_k_stat(2) - binom_var) / (n_proton * (n_proton - 1))
                            dsig_2_list.append(delta_sig_2.val)
                            dsig_2_err_list.append(delta_sig_2.err * np.sqrt(resamples))
                            track_combos = ak.combinations(ak.Array(experiment_tracks[n_proton]), 2)
                            combo_a, combo_b = ak.unzip(track_combos)
                            combos_dphi = combo_a - combo_b
                            v1_list.append(ak.mean(np.cos(combos_dphi)))
                            v2_list.append(ak.mean(np.cos(2 * combos_dphi)))
                            v3_list.append(ak.mean(np.cos(3 * combos_dphi)))

                        with open(f'{var_dir}bin_width{int(np.rad2deg(bin_width) + 0.1)}.txt', 'a') as file:
                            file.write(f'Total Particles {n_track}\n')
                            file.write(', '.join([str(x) for x in n_protons_list]) + '\n')
                            file.write(', '.join([str(x) for x in dsig_2_list]) + '\n')
                            file.write(', '.join([str(x) for x in dsig_2_err_list]) + '\n')
                            file.write(', '.join([str(x) for x in v1_list]) + '\n')
                            file.write(', '.join([str(x) for x in v2_list]) + '\n')
                            file.write(', '.join([str(x) for x in v3_list]) + '\n')
                            file.write('\n')


def plot_full_test_from_file():
    # path = 'N:/UCLA_Microsoft/OneDrive - personalmicrosoftsoftware.ucla.edu/Research/UCLA/Results/momentum_conservation_model/mom_cons_new_pc.txt'
    # path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Results/momentum_conservation_model/mom_cons_new_pc.txt'
    bin_width = 120
    # path = f'C:/Users/Dylan/Desktop/test/pfrac40_convp1_energy2_nevent250.0k/bin_width{bin_width}.txt'
    path = f'C:/Users/Dylan/Desktop/test/pfrac40_convp1_energy2_nevent100.0k/bin_width{bin_width}.txt'
    bin_width_rad = np.deg2rad(bin_width)
    total_tracks, total_protons, dsigs, dsig_errs, v1s, v2s, v3s = [], [], [], [], [], [], []
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
            rolling += 1
        elif rolling == 4:
            v1s.append([float(x) for x in line.strip('\n[]').split(', ')])
            rolling += 1
        elif rolling == 5:
            v2s.append([float(x) for x in line.strip('\n[]').split(', ')])
            rolling += 1
        elif rolling == 6:
            v3s.append([float(x) for x in line.strip('\n[]').split(', ')])
            rolling = 0
        if 'Total Particles ' in line:
            total_tracks.append(int(line.strip().replace('Total Particles ', '')))
            rolling = 1

    dsig_avgs, dsig_avg_errs = [], []
    for tracks, protons, dsigs_i, dsig_errs_i in zip(total_tracks, total_protons, dsigs, dsig_errs):
        if len(dsigs_i) != len(dsig_errs_i):
            print('Error')
            print(tracks, protons, dsigs_i, dsig_errs_i)
        wavg, wavg_err = plot_vs_total_protons(protons, dsigs_i, dsig_errs_i, tracks, plot=False)
        dsig_avgs.append(wavg)
        dsig_avg_errs.append(wavg_err)
    plot_vs_total_particles(total_tracks, dsig_avgs, dsig_avg_errs)
    # plot_fits(total_tracks, dsig_avgs, dsig_avg_errs)
    fig = plt.gcf()
    fig.canvas.manager.set_window_title('dsig2_vs_total_particles_fit')
    ax = fig.gca()
    ax.grid(False)
    fig.subplots_adjust(left=0.165, right=0.995, top=0.995, bottom=0.12)

    dsig_avgs, dsig_avg_errs, dsig_v1, dsig_v1_v2, dsig_v1_v2_v3 = [], [], [], [], []
    for tracks, protons, dsigs_i, dsig_errs_i, v1s_i, v2s_i, v3s_i \
            in zip(total_tracks, total_protons, dsigs, dsig_errs, v1s, v2s, v3s):
        wavg, wavg_err = plot_vs_total_protons(protons, dsigs_i, dsig_errs_i, tracks, plot=False)
        dsig_avgs.append(wavg)
        dsig_avg_errs.append(wavg_err)
        fig2, ax2 = plt.subplots()
        ax2.grid()
        ax2.axhline(0, color='black')
        ax2.set_title(f'{tracks} Total Particles')
        ax2.scatter(protons, v1s_i, label=r'$v_1^2$', zorder=10, color='blue')
        ax2.axhline(np.mean(v1s_i), color='blue', ls='--')
        ax2.scatter(protons, v2s_i, label=r'$v_2^2$', zorder=10, color='orange')
        ax2.axhline(np.mean(v2s_i), color='orange', ls='--')
        ax2.scatter(protons, v3s_i, label=r'$v_3^2$', zorder=10, color='green')
        ax2.axhline(np.mean(v3s_i), color='green', ls='--')
        ax2.set_xlabel('Total Protons')
        ax2.set_ylabel(r'$v_n^2$')
        ax2.legend()
        fig2.tight_layout()
        mean_v1 = np.mean(v1s_i)
        sign_mean_v1 = np.sign(mean_v1)
        mean_v1 = np.sqrt(abs(mean_v1))
        wavg -= vn_divs(bin_width_rad, sign_mean_v1 * mean_v1, 1)
        dsig_v1.append(wavg)
        mean_v2 = np.mean(v2s_i)
        sign_mean_v2 = np.sign(mean_v2)
        mean_v2 = np.sqrt(abs(mean_v2))
        wavg -= vn_divs(bin_width_rad, sign_mean_v2 * mean_v2, 2)
        dsig_v1_v2.append(wavg)
        mean_v3 = np.mean(v3s_i)
        sign_mean_v3 = np.sign(mean_v3)
        mean_v3 = np.sqrt(abs(mean_v3))
        wavg -= vn_divs(bin_width_rad, sign_mean_v3 * mean_v3, 3)
        dsig_v1_v2_v3.append(wavg)

    plot_vs_total_particles(total_tracks, dsig_avgs, dsig_avg_errs, label='No Flow Subtraction', color='black')
    # plot_fits(total_tracks, dsig_avgs, dsig_avg_errs)
    fig = plt.gcf()
    plot_vs_total_particles(total_tracks, dsig_v1, dsig_avg_errs, label=f'$v_1$ Subtraction', color='blue', fig=fig)
    plot_vs_total_particles(total_tracks, dsig_v1_v2, dsig_avg_errs, label=f'$v_1$, $v_2$ Subtraction', color='orange',
                            fig=fig)
    plot_vs_total_particles(total_tracks, dsig_v1_v2_v3, dsig_avg_errs, label=f'$v_1$, $v_2$, $v_3$ Subtraction',
                            color='green', fig=fig)
    fig.canvas.manager.set_window_title('dsig2_vs_total_particles_fit')
    ax = fig.gca()
    ax.set_title(f'{bin_width} Bin Width')
    ax.legend()
    ax.grid(False)
    fig.subplots_adjust(left=0.165, right=0.995, top=0.94, bottom=0.12)


def plot_from_files():
    # path = ('N:/UCLA_Microsoft/OneDrive - personalmicrosoftsoftware.ucla.edu/Research/UCLA/Results/'
    #         'momentum_conservation_model/')
    path = 'C:/Users/Dylan/OneDrive - UCLA IT Services/Research/UCLA/Results/momentum_conservation_model/'
    for dir_name in os.listdir(path):
        dir_path = os.path.join(path, dir_name)
        if os.path.isdir(dir_path):
            for file_name in os.listdir(dir_path):
                file_path = os.path.join(dir_path, file_name)
                if 'bin_width90' not in file_path or 'energy4' in file_path:
                    continue
                print(dir_name, file_name)
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                total_tracks, total_protons, dsigs, dsig_errs = [], [], [], []
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

                dsig_avgs, dsig_avg_errs = [], []
                for tracks, protons, dsigs_i, dsig_errs_i in zip(total_tracks, total_protons, dsigs, dsig_errs):
                    protons, dsigs_i, dsig_errs_i = remove_inf_nan_indices(protons, dsigs_i, dsig_errs_i)

                    wavg, wavg_err = plot_vs_total_protons(protons, dsigs_i, dsig_errs_i, tracks, plot=False)
                    dsig_avgs.append(wavg)
                    dsig_avg_errs.append(wavg_err)
                total_tracks, dsig_avgs, dsig_avg_errs = remove_inf_nan_indices(total_tracks, dsig_avgs, dsig_avg_errs)
                plot_vs_total_particles(total_tracks, dsig_avgs, dsig_avg_errs)
                plot_fits(total_tracks, dsig_avgs, dsig_avg_errs)
                plt.title(f'{dir_name} {file_name}')
                plt.tight_layout()

    plt.show()


def plot_fits(total_tracks, dsig_avgs, dsig_avg_errs):
    fits = {
        # r'$\frac{a}{\sqrt{N}} + b$': sqrt_n_fit,
        # r'$\frac{a}{N} + b$': n_fit,
        r'$\frac{a}{M^c} + b$': n_exp_fit,
    }
    xs = np.linspace(min(total_tracks), max(total_tracks), 1000)
    par_letters = ['a', 'b', 'c']
    for fit_lab, fit in fits.items():
        try:
            popt, pcov = cf(fit, total_tracks, dsig_avgs, sigma=dsig_avg_errs, absolute_sigma=True)
        except ValueError:
            print(total_tracks, dsig_avgs, dsig_avg_errs)
        perr = np.sqrt(np.diag(pcov))
        fit_values = {par: Measure(val, err) for par, val, err in zip(par_letters, popt, perr)}

        plt.plot(xs, fit(xs, *popt), label=fit_lab)

        # Create the textbox
        x_pos, y_pos = 0.4, 0.4  # Set the position of the textbox
        textbox_content = '\n'.join([f"{key} = {value}" for key, value in fit_values.items()])
        textbox = plt.text(x_pos, y_pos, textbox_content, fontsize=12, va='center', ha='left',
                           bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),
                           transform=plt.gca().transAxes)
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
    alpha_init = event.alpha
    fig, ax = plt.subplots()
    ax.grid()
    for x in range(1, 8):
        rng = np.random.default_rng(seed=seed)
        event = PConsSim(5, n_part, rng=rng)
        event.alpha = alpha_init / x
        # print(event.get_m_tracks(20))
        # plot_momenta(event.get_m_tracks(n_part), event.net_momentum_vec)
        # event.rotate_tracks_debug()
        event.rotate_tracks()
        # print(event.get_m_tracks(20))
        # plot_momenta(event.get_m_tracks(n_part), event.net_momentum_vec)
        ax.plot(range(len(event.net_momentum_iterations)), event.net_momentum_iterations,
                label=f'{event.alpha:.2f}')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Event Net Momentum')
    plt.legend()
    plt.tight_layout()
    plt.show()


def check_convergence_vs_alpha():
    n_part = 20
    seed = 588

    rng = np.random.default_rng(seed=seed)
    event = PConsSim(2, n_part, rng=rng)
    plot_momenta(event.momenta, event.net_momentum_vec)
    # plt.show()

    while True:
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        ax.grid()
        ax2.grid()
        fractions = {}
        fraction_alphas = {}
        alphas = np.arange(0.01, 1.0, 0.1)
        for alpha in alphas:
            rng = np.random.default_rng(seed=seed)
            event = PConsSim(2, n_part, rng=rng)
            event.adaptive_alpha = True
            event.max_iterations = 50
            event.convergence_momentum = 0.001
            event.alpha = alpha
            event.rotate_tracks()
            net_ps = event.net_momentum_iterations
            fracs = [net_ps[i] / net_ps[i - 1] if net_ps[i - 1] != 0 else float('inf') for i in range(1, len(net_ps))]
            for iteration, frac in enumerate(fracs):
                if iteration in fractions:
                    fractions[iteration].append(frac)
                    fraction_alphas[iteration].append(alpha)
                else:
                    fractions[iteration] = [frac]
                    fraction_alphas[iteration] = [alpha]
            # if len(fractions) == 0:
            #     fractions = [[frac] for frac in fracs]
            # else:
            #     for i in range(len(fractions)):
            #         fractions[i].append(fracs[i])
            ax2.plot(range(1, len(event.net_momentum_iterations)), fracs, label=f'{alpha:.2f}')
            ax.plot(range(len(event.net_momentum_iterations)), event.net_momentum_iterations,
                    label=f'{alpha:.2f}')
        ax.set_xlabel('Iteration')
        ax2.set_xlabel('Iteration')
        ax.set_ylabel('Event Net Momentum')
        ax2.set_ylabel('Fractional Change in Net Momentum')
        ax.legend()
        ax2.legend()
        fig.tight_layout()
        fig2.tight_layout()

        fig3, ax3 = plt.subplots()
        ax3.grid()
        for it, fracs in fractions.items():
            ax3.plot(fraction_alphas[it], fracs, label=f'Iteration {it}')
        ax3.set_xlabel('Alpha')
        ax3.set_ylabel('Event Net Momentum Fractional Reduction')
        ax3.legend()
        fig3.tight_layout()

        plt.show()
        seed += 1


def test_alpha_adaptive():
    n_part = 20
    seed = 47

    while True:
        if seed % 100 == 0:
            print(seed)
        rng = np.random.default_rng(seed=seed)
        event = PConsSim(2, n_part, rng=rng)
        event.adaptive_alpha = True
        event.convergence_momentum = 0.001
        event.rotate_tracks()
        if len(event.alphas) >= event.max_iterations:
            print('Bad one')
            print(seed)
            fig, ax = plt.subplots()
            ax.grid()
            ax.plot(range(len(event.net_momentum_iterations)), event.net_momentum_iterations)
            plt.show()
        seed += 1


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


def sim_event(n_track, rng, energy, n_protons, y_max, pt_max, p_max, convergence_momentum):
    vector.register_awkward()
    event = PConsSim(energy, n_track, rng=rng)
    event.convergence_momentum = convergence_momentum
    event.adaptive_alpha = True
    event.rotate_tracks()

    while event.init_net_momentum <= event.net_momentum:
        alpha = event.alpha
        event = PConsSim(energy, n_track, rng=rng)
        event.alpha = alpha / 2
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
    if len(dsig_2_list) != len(dsig_2_err_list):
        print('Error')
        print(n_protons, dsig_2_list, dsig_2_err_list)
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


def plot_vs_total_particles(n_tracks, dsig_avgs, dsig_avg_errs, fig=None, label=None, color=None):
    if fig is None:
        plt.figure(figsize=(6, 4), dpi=144)
        plt.grid()
        plt.axhline(0, color='black')
        plt.xlabel(r'Total Particles in Event ($M$)', fontsize=14)
        plt.ylabel(r'$\langle\Delta\sigma^2\rangle$', fontsize=14)
    if label is not None:
        if color is not None:
            plt.errorbar(n_tracks, dsig_avgs, yerr=dsig_avg_errs, ls='none', marker='o', alpha=0.8, label=label,
                         color=color)
        else:
            plt.errorbar(n_tracks, dsig_avgs, yerr=dsig_avg_errs, ls='none', marker='o', alpha=0.8, label=label)
    else:
        if color is not None:
            plt.errorbar(n_tracks, dsig_avgs, yerr=dsig_avg_errs, ls='none', marker='o', alpha=0.8, color=color)
        else:
            plt.errorbar(n_tracks, dsig_avgs, yerr=dsig_avg_errs, ls='none', marker='o', alpha=0.8)
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


def remove_inf_nan_indices(*lists):
    # Check if all input lists have the same length
    list_lengths = [len(lst) for lst in lists]
    if len(set(list_lengths)) != 1:
        raise ValueError("Input lists must have equal lengths")

    # Create a list to store the indices to remove
    indices_to_remove = set()

    # Iterate through the lists and find indices with inf or nan values
    for i in range(list_lengths[0]):
        values_at_index = [lst[i] for lst in lists]
        if any(math.isnan(value) or math.isinf(value) for value in values_at_index):
            indices_to_remove.add(i)

    # Remove the corresponding indices from each list
    cleaned_lists = []
    for lst in lists:
        cleaned_list = [value for index, value in enumerate(lst) if index not in indices_to_remove]
        cleaned_lists.append(cleaned_list)

    return cleaned_lists


def sqrt_n_fit(n, a, b):
    return a / np.sqrt(n) + b


def n_fit(n, a, b):
    return a / n + b


def n_exp_fit(n, a, b, c):
    return a / n**c + b


if __name__ == '__main__':
    main()
