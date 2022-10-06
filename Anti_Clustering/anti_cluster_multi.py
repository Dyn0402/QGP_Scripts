#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 01 1:15 PM 2021
Created in PyCharm
Created as QGP_Scripts/anti_cluster_multi

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import norm
from scipy.stats import uniform
from scipy.interpolate import CubicSpline

from ClustDist import ClustDist, ClustDist_slow, Testrv


def main():
    # specific_vals()
    # random_tracks()
    # clustered_tracks()
    clustered_tracks_ani()
    # rv_test()
    print('donzo')


def specific_vals():
    means = [0.7, 3.5, 2.4]
    sds = [0.2, 0.2, 0.2]
    amps = 1 / (np.sqrt(2 * np.pi) * np.asarray(sds))
    aclust_amp = 0.4
    dists = [norm(mean, sd) for mean, sd in zip(means, sds)]
    x = np.linspace(0, 2 * np.pi, 1000)
    base = 1
    base_int = base * (x[-1] - x[0])
    p_clust = np.ones(1000)
    p_aclust = np.ones(1000)

    fig_ind_clust, (ax_ind_clust, ax_ind_aclust) = plt.subplots(2, 1, sharex=True)
    ax_ind_aclust.set_xlabel('Phi Angle')
    for dist, amp, sd in zip(dists, amps, sds):
        # clust_pdf = base + clust_amp(aclust_amp, sd, base_int) * dist.pdf(x) / amp
        clust_pdf = base + aclust_amp * dist.pdf(x) / amp
        # clust_pdf /= np.sum(clust_pdf) * (x[-1] - x[0]) / len(x)
        aclust_pdf = base - aclust_amp * dist.pdf(x) / amp
        # aclust_pdf /= np.sum(aclust_pdf) * (x[-1] - x[0]) / len(x)
        ax_ind_clust.plot(x, clust_pdf)
        ax_ind_aclust.plot(x, aclust_pdf)
        ax_ind_clust.set_ylim(bottom=0)
        ax_ind_aclust.set_ylim(bottom=0)
        p_clust *= clust_pdf
        p_aclust *= aclust_pdf

    fig_clust, (ax_clust, ax_aclust) = plt.subplots(2, 1, sharex=True, figsize=(6.5, 6.8), dpi=144)
    ax_aclust.set_xlabel('Phi Angle')
    ax_clust.set_ylabel('Next Track Probability')
    ax_aclust.set_ylabel('Next Track Probability')
    ax_clust.axhline(1, color='gray', linewidth=0.5, alpha=0.5)
    ax_aclust.axhline(1, color='gray', linewidth=0.5, alpha=0.5)
    ax_clust.vlines(means, ymin=0.0, ymax=max(p_clust), color='black', ls='--', lw=1.0, alpha=0.4,
                    label='Previous Tracks')
    ax_aclust.vlines(means, ymin=0.0, ymax=max(p_aclust), color='black', ls='--', lw=1.0, alpha=0.4)
    ax_clust.plot(x, p_clust)
    ax_aclust.plot(x, p_aclust)
    # ax_clust.text(4.5, 0.6, 'Attractive', size='x-large')
    # ax_aclust.text(4.5, 0.4, 'Repulsive', size='x-large')
    ax_clust.set_ylim(bottom=0)
    ax_aclust.set_ylim(bottom=0)
    ax_clust.legend(loc='upper right')

    fig_ind_clust.tight_layout()
    fig_clust.tight_layout()

    plt.show()


def random_tracks():
    n_tracks = 50
    track_means = uniform.rvs(size=n_tracks) * 2 * np.pi
    sd = 0.2
    amp = 1 / (np.sqrt(2 * np.pi) * sd)
    aclust_amp = 0.4
    dists = [[norm(mean, sd), norm(mean - 2 * np.pi, sd), norm(mean + 2 * np.pi, sd)] for mean in track_means]
    x = np.linspace(0, 2 * np.pi, 1000)
    base = 1
    p_clust = np.ones(1000)
    p_aclust = np.ones(1000)

    fig_ind_clust, ax_ind_clust = plt.subplots()
    fig_ind_aclust, ax_ind_aclust = plt.subplots()
    ax_ind_clust.set_xlabel('Phi Angle')
    ax_ind_aclust.set_xlabel('Phi Angle')
    for dist in dists:
        clust_pdf = base
        aclust_pdf = base
        for dist_sub in dist:
            clust_pdf += aclust_amp * dist_sub.pdf(x) / amp
            aclust_pdf -= aclust_amp * dist_sub.pdf(x) / amp
        ax_ind_clust.plot(x, clust_pdf)
        ax_ind_aclust.plot(x, aclust_pdf)
        ax_ind_clust.set_ylim(bottom=0)
        ax_ind_aclust.set_ylim(bottom=0)
        p_clust *= clust_pdf
        p_aclust *= aclust_pdf

    fig_clust, ax_clust = plt.subplots()
    fig_aclust, ax_aclust = plt.subplots()
    ax_clust.set_xlabel('Phi Angle')
    ax_aclust.set_xlabel('Phi Angle')
    ax_clust.plot(x, p_clust)
    ax_aclust.plot(x, p_aclust)
    ax_clust.set_ylim(bottom=0)
    ax_aclust.set_ylim(bottom=0)

    fig_ind_clust.tight_layout()
    fig_ind_aclust.tight_layout()
    fig_clust.tight_layout()
    fig_aclust.tight_layout()

    plt.show()


def clustered_tracks():
    n_tracks = 80
    sd = 0.1
    wrap_num = 5
    wrap_ns = [5]
    cl_amp = -0.4
    x = np.linspace(0, 2 * np.pi, 1000)

    # fig_ind_clust, ax_ind_clust = plt.subplots()
    fig_ind_aclust, ax_ind_aclust = plt.subplots()
    # ax_ind_clust.set_xlabel('Phi Angle')
    ax_ind_aclust.set_xlabel('Phi Angle')
    phis = []
    prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2*np.pi, wrap_num=wrap_num)
    while len(phis) < n_tracks:
        phis.append(prob_dist.rvs())
        prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2*np.pi, wrap_num=wrap_num)
        # ax_ind_aclust.plot(x, aclust_pdf)
        ax_ind_aclust.axvline(phis[-1], color='black', ls='--', alpha=0.7)
        # ax_ind_clust.set_ylim(bottom=0)
        ax_ind_aclust.set_ylim(bottom=0)
    # for i in range(100):
    #     ax_ind_aclust.axvline(prob_dist.rvs(), color='black', ls='--', alpha=0.7)
    ax_ind_aclust.set_xlim(left=0, right=2*np.pi)

    # fig_clust, ax_clust = plt.subplots()
    # ax_clust.set_xlabel('Phi Angle')
    # ax_clust.plot(x, p_clust)
    # ax_clust.set_ylim(bottom=0)

    fig_aclust, ax_aclust = plt.subplots()
    ax_aclust.set_xlabel('Phi Angle')
    for wrap_n in wrap_ns:
        prob_dist = ClustDist(phis, sd, cl_amp, a=0, b=2 * np.pi, wrap_num=wrap_n)
        ax_aclust.plot(x, prob_dist.pdf(x), label=f'Wrap_num = {wrap_n}')
        print(prob_dist.pdf(x))
    ax_aclust.set_ylim(bottom=0, top=1.1*max(prob_dist.pdf(x)))
    ax_aclust.legend()
    fig_aclust.tight_layout()

    # fig_ind_clust.tight_layout()
    fig_ind_aclust.tight_layout()
    # fig_clust.tight_layout()
    fig_aclust.tight_layout()

    plt.show()


def clustered_tracks_ani():
    # gif_dir = 'C:/Users/Dylan/Desktop/pdf_test3.gif'
    gif_dir = 'F:/Research/Results/Presentations/10-3-22/cluster_example.gif'
    # gif_dir = 'C:/Users/Dyn04/Desktop/pdf_test.gif'
    fps = 0.8
    n_tracks = 25
    sd = 1.0  # np.pi
    wrap_num = 5
    cl_amp = 0.5
    x = np.linspace(0, 2 * np.pi, 1000)
    gif_path = f'{gif_dir}'

    phis = []
    prob_dists = [ClustDist(phis, sd, cl_amp, a=0, b=2*np.pi, wrap_num=wrap_num)]
    y_maxes = []
    while len(phis) < n_tracks:
        phis.append(prob_dists[-1].rvs())
        prob_dists.append(ClustDist(phis, sd, cl_amp, a=0, b=2*np.pi, wrap_num=wrap_num))
        y_maxes.append(max(prob_dists[-1].pdf(x)))
    fig = plt.figure(figsize=(6.67, 4), dpi=144)
    ax = plt.subplot()

    y_lim = (0, max(y_maxes) * 1.2)

    ani = FuncAnimation(fig, ani_pdf, frames=range(len(phis)), interval=1.0 / fps * 1000, repeat_delay=5000,
                        repeat=False, fargs=(phis, prob_dists, x, y_lim, {'sd': sd, 'amp': -cl_amp}, ax))
    ani.save(gif_path, dpi=100, writer=PillowWriter(fps=fps))


def ani_pdf(phi_index, phis, prob_dists, x_vals, y_lim, sim_pars, ax):
    print(phi_index)
    plot_pdf_ani(phis[:phi_index], prob_dists[phi_index], x_vals, y_lim, sim_pars, ax)


def plot_pdf_ani(phis, pdf, x_vals, y_lim, sim_pars, ax):
    ax.clear()
    pdf_vals = pdf.pdf(x_vals)
    ax.plot(x_vals, pdf_vals, label='PDF for Next Track')
    ax.vlines(phis[:-1], 0, y_lim[-1] * 0.5, color='black', ls='--', label='Tracks')
    if len(phis) > 0:
        ax.axvline(phis[-1], 0, y_lim[-1] * 0.5, color='red', ls='--', label='New Track')
    ax.set_title(f'Probability Distribution for Track #{len(phis)}  A={sim_pars["amp"]}, σ={sim_pars["sd"]:.1f}')
    ax.set_xlabel('ϕ')
    ax.set_ylabel('P(ϕ)')
    ax.set_ylim(y_lim)
    ax.legend()
    plt.tight_layout()


def rv_test():
    dist = Testrv(a=0, b=2)
    x = np.linspace(-3, 3, 1000)
    plt.plot(x, dist.pdf(x))
    plt.hist(dist.rvs(size=10000), density=True)
    plt.show()


def clust_amp(aclust_amp, sigma, base_int):
    return 1 / (1 / aclust_amp - 2 * sigma * np.sqrt(2 * np.pi) / base_int)


if __name__ == '__main__':
    main()
