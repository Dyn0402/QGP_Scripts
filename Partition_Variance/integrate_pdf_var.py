#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 17 7:23 PM 2023
Created in PyCharm
Created as QGP_Scripts/integrate_pdf_var

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, nquad
from scipy.stats import norm
from scipy.optimize import curve_fit as cf
from scipy.fft import fft, ifft, fftfreq

import uproot
from scipy.interpolate import CubicSpline

from Binom_Slices.analyze_binom_slices import *


def main():
    # vn_test()
    # gaus_test()
    # gaus_fit()
    # gaus_test_fits()
    # gaus_fourier_decomp()
    # flow_fourier_decomp()
    # gaus_fit_dependence()
    # convo_test()
    # eff_v2_combo()
    # eff_v2_combo2()
    # eff_gaus_combo()
    # eff_gaus_combo2()
    gaus_v2_combo()
    # eff_plotting()
    # v2_plotting()
    # vn_analytic_plotting()
    # recursive_integration_test()
    print('donzo')


def vn_test():
    func = vn_pdf
    n = 2
    v = 0.07
    psi = np.pi / 3
    func_args = (psi, v, n)

    plot_variance(func, func_args)


def gaus_test():
    mu = np.pi
    sigma = 0.1
    amp = 0.5

    func = base_gaus_pdf
    func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))

    plot_variance(func, func_args)


def recursive_integration_test():
    mu = np.pi
    sigma = 0.8
    amp = 0.4
    func1 = base_gaus_pdf_wrap
    func1_args = [mu, sigma, amp, 1. / get_norm_scipy(func1, (mu, sigma, amp, 1))]

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = [psi, v2, n]

    func = lambda x, mu_, sigma_, amp_, c_, v2_, psi_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf_wrap(x, mu_, sigma_, amp_, c_) * vn_pdf(x, psi_, v2_, n_)
    func_args = [*func1_args, *func2_args, 1. / get_norm_scipy(func, (*func1_args, *func2_args, 1))]

    func2 = lambda x, psi_1_, psi_2_, sigma_, amp_, c_, v2_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf_wrap(x, psi_1_, sigma_, amp_, c_) * vn_pdf(x, psi_2_, v2_, n_)
    func2_args = [mu, psi, sigma, amp, func1_args[-1], v2, n, func_args[-1]]

    width = np.deg2rad(201)

    def avg_orient(psi_in, func_in, func_args_in, width_in):
        func_args_in[5] = psi_in
        func_args_in[-1] = 1
        func_args_in[-1] = 1. / get_norm_scipy(func_in, func_args_in)
        return get_partition_variance_scipy(func_in, tuple(func_args_in), width_in)[0]

    print(quad(avg_orient, 0, 2 * np.pi, args=(func, func_args, width))[0] / (2 * np.pi))

    print(nquad(integrate_partition, [[0, 2 * np.pi], [0, 2 * np.pi]], args=((func2, width, func2_args),))[0])


def gaus_fit():
    mu = np.pi
    sigma = 0.3
    amp = 0.5

    # func = norm(mu, sigma).pdf
    # func_args = ()
    func = base_gaus_pdf
    func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))
    bounds = (0, 2 * np.pi)

    plot_pdf(func, func_args)
    widths = np.linspace(*bounds, 100)
    pp_minus_p2_terms, pp_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds)
        pp_minus_p2_terms.append(pp_minus_p2)
        pp_terms.append(pp)
    fig, ax = plt.subplots()
    ax.axhline(0, color='black')
    ax.plot(widths, pp_minus_p2_terms)
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')

    width_fit_low, width_fit_high = np.deg2rad([60, 300])
    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    width_filter = (width_fit_low < widths) & (widths < width_fit_high)
    # print(widths[width_filter], np.array(lin_terms)[width_filter])
    popt_quad, pcov_quad = cf(quad_180_rad, widths[width_filter], np.array(pp_minus_p2_terms)[width_filter])
    xs_fit_plot = np.linspace(min(widths), max(widths), 1000)
    ax.plot(xs_fit_plot, quad_180_rad(xs_fit_plot, *popt_quad), color='red', label='Quadratic')
    ax.axvline(width_fit_low, ls='--', color='orange', label='Fit Range')
    ax.axvline(width_fit_high, ls='--', color='orange')
    ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)

    p0 = (0.075, 0.2)
    popt_cos, pcov = cf(cos_pi, widths[width_filter], np.array(pp_minus_p2_terms)[width_filter], p0=p0)
    # ax.plot(xs_fit_plot, cos_pi(xs_fit_plot, *p0), color='purple', alpha=0.3, label='Sine p0')
    ax.plot(xs_fit_plot, cos_pi(xs_fit_plot, *popt_cos), color='purple', label='Sine')
    # ax.plot(xs_fit_plot, 5 * cos_pi(xs_fit_plot, *popt)**2, color='green', label='Sine**2')

    # p0 = (*popt_quad, *popt_cos, 0.5)
    p0 = (-3e-4, 0.002, 0.0004, 0.5)
    popt, pcov = cf(quad_cos, widths[width_filter], np.array(pp_minus_p2_terms)[width_filter], p0=p0)
    ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *p0), color='olive', alpha=0.4, label='Combo Initial Guess')
    ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *popt), color='olive', label='Combo')
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    fig.tight_layout()

    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    width_fit_low, width_fit_high = np.deg2rad([0, 360])
    width_filter = (width_fit_low < widths) & (widths < width_fit_high)
    ax.axvline(width_fit_low, ls='--', color='orange', label='Fit Range')
    ax.axvline(width_fit_high, ls='--', color='orange')
    ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)
    p0 = (-3e-4, 0.002, 0.0004, 0.5)
    popt, pcov = cf(quad_cos, widths[width_filter], np.array(pp_minus_p2_terms)[width_filter], p0=p0)
    ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *p0), color='olive', alpha=0.7, ls='--', label='Combo Initial Guess')
    ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *popt), color='olive', label='Combo')
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    fig.tight_layout()

    plt.show()


def gaus_fourier_decomp():
    mu = np.pi
    sigmas = [0.3, 0.8]
    amp = 0.5

    bounds = (0, 2 * np.pi)
    width_points = 100
    integration_points = 100000  # 100000

    fig, ax = plt.subplots(dpi=144)
    fig2, ax2 = plt.subplots(dpi=144)
    fig3, ax3 = plt.subplots(dpi=144)
    colors = iter(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    for sigma in sigmas:
        func = base_gaus_pdf
        func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))

        widths = np.linspace(*bounds, width_points)
        pp_minus_p2_terms, pp_terms = [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds, points=integration_points)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_terms.append(pp)
        pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)

        # pp_minus_p2_norm = pp_minus_p2_terms - pp_minus_p2_terms.max() / 2
        # pp_minus_p2_norm = np.int16((pp_minus_p2_norm / pp_minus_p2_norm.max()) * 32767)
        pp_minus_p2_norm = pp_minus_p2_terms

        n_points = (width_points - 1)

        yf = fft(pp_minus_p2_norm[:-1])
        xf = fftfreq(n_points, widths[1] - widths[0])
        xt = 1 / xf

        def fit_f(x, a, b):
            return a * (np.exp(b * abs(x)) - 1)

        p0 = [100, 1.2]
        x_filter = (abs(xt) > 1.5) & np.isfinite(abs(xt))
        popt, pcov = cf(fit_f, xt[x_filter], abs(yf[x_filter]), p0=p0)

        c = next(colors)
        ax.scatter(xt, np.abs(yf), label=f'sigma={sigma}', color=c)
        xs_fit = np.linspace(-7, 7, 1000)
        ax.plot(xs_fit, fit_f(xs_fit, *popt), ls='--', color=c)
        print(sigma, popt)
        # ax.plot(xs_fit, fit_f(xs_fit, *p0), ls='-', alpha=0.3, color=c)

        def fit_f(x, a, b):
            return a * np.exp(-(x / b)**2)

        p0 = [0.003, 0.1]
        # x_filter = (abs(xf) > 1.5) & np.isfinite(abs(xt))
        popt, pcov = cf(fit_f, xf, abs(yf), p0=p0)

        c = next(colors)
        ax2.scatter(xf, np.abs(yf), label=f'sigma={sigma}', color=c)
        xs_fit = np.linspace(-1, 1, 1000)
        ax2.plot(xs_fit, fit_f(xs_fit, *popt), ls='--', color=c)
        ax2.plot(xs_fit, fit_f(xs_fit, *p0), ls='-', alpha=0.2, color=c)
        print(sigma, popt)

        ax3.plot(xf, np.abs(yf) - fit_f(xf, *popt), marker='o', alpha=0.6, color=c)

    ax.set_xlim(-7, 7)
    ax.legend()
    fig.tight_layout()

    ax2.set_xlim(-1, 1)
    ax2.legend()
    fig2.tight_layout()

    plt.show()


def flow_fourier_decomp():
    n = 3
    v = 0.07
    psi = np.pi / 3

    bounds = (0, 2 * np.pi)
    width_points = 100

    func = vn_pdf
    func_args = (psi, v, n)

    widths = np.linspace(*bounds, width_points)
    pp_minus_p2_terms, pp_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds, points=100000)
        pp_minus_p2_terms.append(pp_minus_p2)
        pp_terms.append(pp)
    pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)

    pp_minus_p2_norm = pp_minus_p2_terms - pp_minus_p2_terms.max() / 2
    pp_minus_p2_norm = np.int16((pp_minus_p2_norm / pp_minus_p2_norm.max()) * 32767)

    fig, ax = plt.subplots(dpi=144)

    n_points = width_points - 1  # == len(widths[:-1]) == len(pp_minus_p2_norm[:-1])

    # fig, ax = plt.subplots(dpi=144)
    # ax.plot(rep_widths, rep_pp_m_p2)

    yf = fft(pp_minus_p2_norm[:-1])
    xf = fftfreq(n_points, widths[1] - widths[0])

    # fig, ax = plt.subplots(dpi=144)
    # ax.plot(xf, np.abs(yf))

    # fig, ax = plt.subplots(dpi=144)
    ax.scatter(1 / xf, np.abs(yf))

    ax.set_xlim(-7, 7)
    ax.legend()
    fig.tight_layout()

    plt.show()


def gaus_test_fits():
    mu = np.pi
    sigma = 0.3
    amp = 0.5

    func = base_gaus_pdf
    func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))
    bounds = (0, 2 * np.pi)

    widths = np.linspace(*bounds, 100)
    pp_minus_p2_terms, pp_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds, points=100000)
        pp_minus_p2_terms.append(pp_minus_p2)
        pp_terms.append(pp)
    pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)
    print(', '.join([str(x) for x in widths]))
    print(', '.join([str(x) for x in pp_minus_p2_terms]))

    xs_fit_plot = np.linspace(min(widths), max(widths), 1000)

    # p0 = [-3e-4, 0.002, 0.0004, 0.5]  # quad_cos
    # p0 = [0.0003, 0.1, 5]  # cos_sin_gaus
    # p0 = [0.0003]  # cos_pi_fixed
    p0 = [0.0003, 0.2, 0.02, 0.002]  # cos_pi_fixed
    func_fit = cos_sin4

    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)
    # popt, pcov = cf(quad_cos, widths, np.array(pp_minus_p2_terms), p0=p0)

    popt, pcov = cf(func_fit, widths, pp_minus_p2_terms, p0=p0)
    ax.plot(xs_fit_plot, func_fit(xs_fit_plot, *p0), color='olive', alpha=0.7, ls='--', label='Combo Initial Guess')
    ax.plot(xs_fit_plot, func_fit(xs_fit_plot, *popt), color='olive', label='Combo')
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    fig.tight_layout()

    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms - func_fit(widths, *popt))

    pp_minus_p2_terms, pp_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds, points=10000)
        pp_minus_p2_terms.append(pp_minus_p2)
        pp_terms.append(pp)
    pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)
    ax.scatter(widths, pp_minus_p2_terms - func_fit(widths, *popt))

    pp_minus_p2_terms, pp_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds, points=1000)
        pp_minus_p2_terms.append(pp_minus_p2)
        pp_terms.append(pp)
    pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)
    ax.scatter(widths, pp_minus_p2_terms - func_fit(widths, *popt))

    print(popt)

    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel('Fit Deviation')
    ax.grid()
    fig.tight_layout()

    plt.show()


def gaus_fit_dependence():
    mu = np.pi
    # sigmas = np.linspace(0.3, 1.0, 10)
    sigmas = [0.3]
    amp = 0.5
    bounds = (0, 2 * np.pi)
    func = base_gaus_pdf
    widths = np.linspace(*bounds, 100)
    xs_fit_plot = np.linspace(min(widths), max(widths), 1000)
    plot_pars = False

    quad_mags, quad_consts, cos_mags, gaus_sigs = [], [], [], []
    quad_mags_err, quad_consts_err, cos_mags_err, gaus_sigs_err = [], [], [], []
    for sigma in sigmas:
        func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))

        pp_minus_p2_terms, pp_terms = [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_terms.append(pp)

        p0 = (-1e-4, 0.0006, 0.004, 3)
        pp_minus_p2_terms = np.nan_to_num(pp_minus_p2_terms)
        popt, pcov = cf(quad_cos, widths, np.array(pp_minus_p2_terms), p0=p0)
        perr = np.sqrt(np.diag(pcov))
        quad_mags.append(popt[0])
        quad_consts.append(popt[1])
        cos_mags.append(popt[2])
        gaus_sigs.append(popt[3])
        quad_mags_err.append(perr[0])
        quad_consts_err.append(perr[1])
        cos_mags_err.append(perr[2])
        gaus_sigs_err.append(perr[3])

        fig, ax = plt.subplots(dpi=144)
        ax.axhline(0, color='black')
        ax.scatter(widths, pp_minus_p2_terms)
        ax.set_xlabel('Partition Width (w)')
        ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
        ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)
        ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *p0), color='olive', alpha=0.7, ls='--', label='Combo Initial Guess')
        ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *popt), color='olive', label='Combo')
        ax.set_title(f'Fit for sigma={sigma:.2f}')
        ax.legend()
        fig.tight_layout()

    if plot_pars:
        fig_quad_mag, ax_quad_mag = plt.subplots(dpi=144)
        fig_quad_consts, ax_quad_consts = plt.subplots(dpi=144)
        fig_cos_mags, ax_cos_mags = plt.subplots(dpi=144)
        fig_gaus_sigs, ax_gaus_sigs = plt.subplots(dpi=144)
        ax_quad_mag.errorbar(sigmas, quad_mags, yerr=quad_mags_err, ls='none', marker='o')
        ax_quad_consts.errorbar(sigmas, quad_consts, yerr=quad_consts_err, ls='none', marker='o')
        ax_cos_mags.errorbar(sigmas, cos_mags, yerr=cos_mags_err, ls='none', marker='o')
        ax_gaus_sigs.errorbar(sigmas, gaus_sigs, yerr=gaus_sigs_err, ls='none', marker='o')
        ax_quad_mag.set_xlabel('Gaussian Cluster Sigma')
        ax_quad_consts.set_xlabel('Gaussian Cluster Sigma')
        ax_cos_mags.set_xlabel('Gaussian Cluster Sigma')
        ax_gaus_sigs.set_xlabel('Gaussian Cluster Sigma')
        ax_quad_mag.set_ylabel('Quadratic Magnitude')
        ax_quad_consts.set_ylabel('Quadratic Constant')
        ax_cos_mags.set_ylabel('Cosine Magnitude')
        ax_gaus_sigs.set_ylabel('Mixing Gaussian Sigma')
        fig_quad_mag.tight_layout()
        fig_quad_consts.tight_layout()
        fig_cos_mags.tight_layout()
        fig_gaus_sigs.tight_layout()

    plt.show()


def convo_test():
    mu = np.pi
    sigma = 0.1
    amp = 0.5
    func1 = base_gaus_pdf
    func1_args = (mu, sigma, amp, 1. / get_norm(func1, (mu, sigma, amp, 1)))
    func1_name = 'Gaussian Cluster'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = (psi, v2, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, mu_, sigma_, amp_, c_, v2_, psi_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf(x, mu_, sigma_, amp_, c_) * vn_pdf(x, psi_, v2_, n_)
    func3_args = (*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1)))
    func3_name = 'Convolution'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pp, ax_pp = plt.subplots(dpi=144, figsize=(6, 3))
    fig_norm, ax_norm = plt.subplots(dpi=144, figsize=(7, 3))
    fig_norm_convo, ax_norm_convo = plt.subplots(dpi=144, figsize=(7, 3))
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax_norm.axhline(0, color='black')
    ax_norm_convo.axhline(0, color='black')
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_minus_p2_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        pp_terms, pp_minus_p2_terms, pp_minus_p2_norm_terms = [], [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width)
            p = width / (2 * np.pi)
            pp_terms.append(pp)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_minus_p2_norm_terms.append(pp_minus_p2 / (p * (1 - p)))
        ax_norm.plot(widths, pp_minus_p2_norm_terms, label=func_name)
        if func_name == 'Convolution':
            ax_norm_convo.plot(widths, pp_minus_p2_norm_terms, color='green', label=func_name)
        ax.plot(widths, pp_minus_p2_terms, label=func_name)
        ax_pp.plot(widths, pp_terms, label=func_name)
        pp_minus_p2_dict.update({func_name: pp_minus_p2_terms})
        # pp_minus_p2_dict.update({func_name: pp_terms})
    ax_norm.set_xlabel('Partition Width (w)')
    ax_norm.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    ax_pp.legend()
    fig_norm.tight_layout()
    fig_norm_convo.tight_layout()
    fig.tight_layout()
    fig_pp.tight_layout()

    fig_ft, ax_ft = plt.subplots(dpi=144)
    yf_dict, yf_plt_dict = {}, {}
    for func_name, y in pp_minus_p2_dict.items():
        # print(f'{func_name}:\n{y}')
        y = np.nan_to_num(y)
        # print(f'{func_name}:\n{y}')
        n = len(y)
        yf = fft(y)
        xf = fftfreq(n, 2 * np.pi / n)[:n//2]
        ax_ft.plot(xf, 2.0 / n * np.abs(yf[0:n//2]), label=func_name)
        yf_dict.update({func_name: yf})
        yf_plt_dict.update({func_name: 2.0 / n * np.abs(yf[0:n//2])})
    # ax_ft.plot(xf, 2.0 / n * np.abs((yf_dict['Convolution'] / yf_dict['Elliptic Flow'])[0:n//2]),
    #            label='Convolution / Flow')
    # ax_ft.plot(xf, yf_plt_dict['Convolution'] / yf_plt_dict['Elliptic Flow'], label='Convolution / Flow')
    ax_ft.legend()
    ax_ft.grid()
    ax_ft.set_xlabel('Frequency')
    ax_ft.set_ylabel('Power')
    fig_ft.tight_layout()

    plt.show()


def eff_v2_combo():
    func1_pre = get_efficiency_pdf(62)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x: func1_pre(x) / func1_norm
    func1_args = ()
    func1_name = 'Efficiency'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = (psi, v2, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, v2_, psi_, n_, c_combo_: c_combo_ * func1(x) * vn_pdf(x, psi_, v2_, n_)
    func3_args = (*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1)))
    func3_name = 'Combination'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pp, ax_pp = plt.subplots(dpi=144, figsize=(6, 3))
    fig_norm, ax_norm = plt.subplots(dpi=144, figsize=(7, 3))
    fig_norm_convo, ax_norm_convo = plt.subplots(dpi=144, figsize=(7, 3))
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax_norm.axhline(0, color='black')
    ax_norm_convo.axhline(0, color='black')
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_minus_p2_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        pp_terms, pp_minus_p2_terms, pp_minus_p2_norm_terms = [], [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width)
            p = width / (2 * np.pi)
            pp_terms.append(pp)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_minus_p2_norm_terms.append(pp_minus_p2 / (p * (1 - p)))
        ax_norm.plot(widths, pp_minus_p2_norm_terms, label=func_name)
        if func_name == 'Combination':
            ax_norm_convo.plot(widths, pp_minus_p2_norm_terms, color='green', label=func_name)
        ax.plot(widths, pp_minus_p2_terms, label=func_name)
        ax_pp.plot(widths, pp_terms, label=func_name)
        pp_minus_p2_dict.update({func_name: pp_minus_p2_terms})
        # pp_minus_p2_dict.update({func_name: pp_terms})
    ax_norm.set_xlabel('Partition Width (w)')
    ax_norm.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    ax_pp.legend()
    fig_norm.tight_layout()
    fig_norm_convo.tight_layout()
    fig.tight_layout()
    fig_pp.tight_layout()

    fig_v2_div_eff, ax_v2_div_eff = plt.subplots(dpi=144, figsize=(7, 3))
    ax_v2_div_eff.axhline(0, color='black')
    for func_name, y in pp_minus_p2_dict.items():
        ax_v2_div_eff.plot(widths, y, label=func_name)
    combo_div_eff = np.array(pp_minus_p2_dict['Combination']) - np.array(pp_minus_p2_dict['Efficiency'])
    ax_v2_div_eff.plot(widths, combo_div_eff, label='Combo - Eff')
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_v2_div_eff.legend()
    fig_v2_div_eff.tight_layout()

    fig_ft, ax_ft = plt.subplots(dpi=144)
    yf_dict, yf_plt_dict = {}, {}
    for func_name, y in pp_minus_p2_dict.items():
        # print(f'{func_name}:\n{y}')
        y = np.nan_to_num(y)
        # print(f'{func_name}:\n{y}')
        n = len(y)
        yf = fft(y)
        xf = fftfreq(n, 2 * np.pi / n)[:n//2]
        ax_ft.plot(xf, 2.0 / n * np.abs(yf[0:n//2]), label=func_name)
        yf_dict.update({func_name: yf})
        yf_plt_dict.update({func_name: 2.0 / n * np.abs(yf[0:n//2])})
    # ax_ft.plot(xf, 2.0 / n * np.abs((yf_dict['Convolution'] / yf_dict['Elliptic Flow'])[0:n//2]),
    #            label='Convolution / Flow')
    # ax_ft.plot(xf, yf_plt_dict['Convolution'] / yf_plt_dict['Elliptic Flow'], label='Convolution / Flow')
    ax_ft.legend()
    ax_ft.grid()
    ax_ft.set_xlabel('Frequency')
    ax_ft.set_ylabel('Power')
    fig_ft.tight_layout()

    plt.show()


def eff_v2_combo2():
    func1_pre = get_efficiency_pdf(62)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x: func1_pre(x) / func1_norm
    func1_args = ()
    func1_name = 'Efficiency'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = (psi, v2, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, v2_, psi_, n_, c_combo_: c_combo_ * func1(x) * vn_pdf(x, psi_, v2_, n_)
    func3_args = [*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1))]
    func3_name = 'Combination'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pp, ax_pp = plt.subplots(dpi=144, figsize=(6, 3))
    fig_norm, ax_norm = plt.subplots(dpi=144, figsize=(7, 3))
    fig_norm_convo, ax_norm_convo = plt.subplots(dpi=144, figsize=(7, 3))
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax_norm.axhline(0, color='black')
    ax_norm_convo.axhline(0, color='black')
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_minus_p2_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        pp_terms, pp_minus_p2_terms, pp_minus_p2_norm_terms = [], [], []
        for width in widths:
            if func_name == 'Combination':
                pp_minus_p2_psi, pp_psi = [], []
                for psi in np.linspace(0, 2 * np.pi, 100):
                    func_args[1] = psi
                    pp_minus_p2, pp = get_partition_variance(func, func_args, width)
                    pp_minus_p2_psi.append(pp_minus_p2)
                    pp_psi.append(pp)
                pp_minus_p2 = np.mean(pp_minus_p2_psi)
                pp = np.mean(pp)
            else:
                pp_minus_p2, pp = get_partition_variance(func, func_args, width)
            p = width / (2 * np.pi)
            pp_terms.append(pp)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_minus_p2_norm_terms.append(pp_minus_p2 / (p * (1 - p)))
        ax_norm.plot(widths, pp_minus_p2_norm_terms, label=func_name)
        if func_name == 'Combination':
            ax_norm_convo.plot(widths, pp_minus_p2_norm_terms, color='green', label=func_name)
        ax.plot(widths, pp_minus_p2_terms, label=func_name)
        ax_pp.plot(widths, pp_terms, label=func_name)
        pp_minus_p2_dict.update({func_name: pp_minus_p2_terms})
        # pp_minus_p2_dict.update({func_name: pp_terms})
    ax_norm.set_xlabel('Partition Width (w)')
    ax_norm.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    ax_pp.legend()
    fig_norm.tight_layout()
    fig_norm_convo.tight_layout()
    fig.tight_layout()
    fig_pp.tight_layout()

    fig_v2_div_eff, ax_v2_div_eff = plt.subplots(dpi=144, figsize=(7, 3))
    ax_v2_div_eff.axhline(0, color='black')
    for func_name, y in pp_minus_p2_dict.items():
        ax_v2_div_eff.plot(widths, y, label=func_name)
    combo_div_eff = np.array(pp_minus_p2_dict['Combination']) - np.array(pp_minus_p2_dict['Efficiency'])
    ax_v2_div_eff.plot(widths, combo_div_eff, label='Combo - Eff')
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_v2_div_eff.legend()
    fig_v2_div_eff.tight_layout()

    fig_ft, ax_ft = plt.subplots(dpi=144)
    yf_dict, yf_plt_dict = {}, {}
    for func_name, y in pp_minus_p2_dict.items():
        # print(f'{func_name}:\n{y}')
        y = np.nan_to_num(y)
        # print(f'{func_name}:\n{y}')
        n = len(y)
        yf = fft(y)
        xf = fftfreq(n, 2 * np.pi / n)[:n//2]
        ax_ft.plot(xf, 2.0 / n * np.abs(yf[0:n//2]), label=func_name)
        yf_dict.update({func_name: yf})
        yf_plt_dict.update({func_name: 2.0 / n * np.abs(yf[0:n//2])})
    # ax_ft.plot(xf, 2.0 / n * np.abs((yf_dict['Convolution'] / yf_dict['Elliptic Flow'])[0:n//2]),
    #            label='Convolution / Flow')
    # ax_ft.plot(xf, yf_plt_dict['Convolution'] / yf_plt_dict['Elliptic Flow'], label='Convolution / Flow')
    ax_ft.legend()
    ax_ft.grid()
    ax_ft.set_xlabel('Frequency')
    ax_ft.set_ylabel('Power')
    fig_ft.tight_layout()

    plt.show()


def eff_gaus_combo():
    func1_pre = get_efficiency_pdf(62)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x: func1_pre(x) / func1_norm
    func1_args = ()
    func1_name = 'Efficiency'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    mu = np.pi
    sigma = 0.1
    amp = 0.5
    func2 = base_gaus_pdf
    func2_args = (mu, sigma, amp, 1. / get_norm(func2, (mu, sigma, amp, 1)))
    func2_name = 'Gaussian Cluster'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, mu_, sigma_, amp_, c_, c_combo_: c_combo_ * func1(x) * base_gaus_pdf(x, mu_, sigma_, amp_, c_)
    func3_args = (*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1)))
    func3_name = 'Combination'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pp, ax_pp = plt.subplots(dpi=144, figsize=(6, 3))
    fig_norm, ax_norm = plt.subplots(dpi=144, figsize=(7, 3))
    fig_norm_convo, ax_norm_convo = plt.subplots(dpi=144, figsize=(7, 3))
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax_norm.axhline(0, color='black')
    ax_norm_convo.axhline(0, color='black')
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_minus_p2_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        pp_terms, pp_minus_p2_terms, pp_minus_p2_norm_terms = [], [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width)
            p = width / (2 * np.pi)
            pp_terms.append(pp)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_minus_p2_norm_terms.append(pp_minus_p2 / (p * (1 - p)))
        ax_norm.plot(widths, pp_minus_p2_norm_terms, label=func_name)
        if func_name == 'Combination':
            ax_norm_convo.plot(widths, pp_minus_p2_norm_terms, color='green', label=func_name)
        ax.plot(widths, pp_minus_p2_terms, label=func_name)
        ax_pp.plot(widths, pp_terms, label=func_name)
        pp_minus_p2_dict.update({func_name: pp_minus_p2_terms})
        # pp_minus_p2_dict.update({func_name: pp_terms})
    ax_norm.set_xlabel('Partition Width (w)')
    ax_norm.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    ax_pp.legend()
    fig_norm.tight_layout()
    fig_norm_convo.tight_layout()
    fig.tight_layout()
    fig_pp.tight_layout()

    fig_v2_div_eff, ax_v2_div_eff = plt.subplots(dpi=144, figsize=(7, 3))
    ax_v2_div_eff.axhline(0, color='black')
    for func_name, y in pp_minus_p2_dict.items():
        ax_v2_div_eff.plot(widths, y, label=func_name)
    combo_div_eff = np.array(pp_minus_p2_dict['Combination']) - np.array(pp_minus_p2_dict['Efficiency'])
    ax_v2_div_eff.plot(widths, combo_div_eff, label='Combo - Eff')
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_v2_div_eff.legend()
    fig_v2_div_eff.tight_layout()

    fig_ft, ax_ft = plt.subplots(dpi=144)
    yf_dict, yf_plt_dict = {}, {}
    for func_name, y in pp_minus_p2_dict.items():
        # print(f'{func_name}:\n{y}')
        y = np.nan_to_num(y)
        # print(f'{func_name}:\n{y}')
        n = len(y)
        yf = fft(y)
        xf = fftfreq(n, 2 * np.pi / n)[:n//2]
        ax_ft.plot(xf, 2.0 / n * np.abs(yf[0:n//2]), label=func_name)
        yf_dict.update({func_name: yf})
        yf_plt_dict.update({func_name: 2.0 / n * np.abs(yf[0:n//2])})
    # ax_ft.plot(xf, 2.0 / n * np.abs((yf_dict['Convolution'] / yf_dict['Elliptic Flow'])[0:n//2]),
    #            label='Convolution / Flow')
    # ax_ft.plot(xf, yf_plt_dict['Convolution'] / yf_plt_dict['Elliptic Flow'], label='Convolution / Flow')
    ax_ft.legend()
    ax_ft.grid()
    ax_ft.set_xlabel('Frequency')
    ax_ft.set_ylabel('Power')
    fig_ft.tight_layout()

    plt.show()


def eff_gaus_combo2():
    func1_pre = get_efficiency_pdf(62)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x: func1_pre(x) / func1_norm
    func1_args = ()
    func1_name = 'Efficiency'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    mu = np.pi
    sigma = 0.8
    amp = 0.1
    # func2 = base_gaus_pdf
    func2 = base_gaus_pdf_wrap
    # func2_args = (mu, sigma, amp, 1. / get_norm(func2, (mu, sigma, amp, 1)))
    func2_args = (mu, sigma, amp)
    func2_name = 'Gaussian Cluster'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    # func3 = lambda x, mu_, sigma_, amp_, c_, c_combo_: c_combo_ * func1(x) * base_gaus_pdf(x, mu_, sigma_, amp_, c_)
    func3 = lambda x, mu_, sigma_, amp_, c_combo_: c_combo_ * func1(x) * base_gaus_pdf_wrap(x, mu_, sigma_, amp_)
    func3_args = [*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1))]
    func3_name = 'Combination'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pdfs, ax_pdfs = plt.subplots(dpi=144, figsize=(8, 4))
    xs = np.linspace(0, 2 * np.pi, 1000)
    ax_pdfs.plot(xs, func1(xs, *func1_args), label=func1_name)
    ax_pdfs.plot(xs, func2(xs, *func2_args), label=func2_name)
    # ax_pdfs.plot(xs, func3(xs, *func3_args), label=func3_name)
    ax_pdfs.set_xlabel(r'$\phi$')
    ax_pdfs.set_ylabel('Probability')
    ax_pdfs.set_ylim(bottom=0)
    ax_pdfs.legend()
    fig_pdfs.tight_layout()

    fig_pp, ax_pp = plt.subplots(dpi=144, figsize=(6, 3))
    fig_norm, ax_norm = plt.subplots(dpi=144, figsize=(7, 3))
    fig_norm_convo, ax_norm_convo = plt.subplots(dpi=144, figsize=(7, 3))
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax_norm.axhline(0, color='black')
    ax_norm_convo.axhline(0, color='black')
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_minus_p2_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        pp_terms, pp_minus_p2_terms, pp_minus_p2_norm_terms = [], [], []
        for width in widths:
            if func_name == 'Combination':
                pp_minus_p2_psi, pp_psi = [], []
                for psi in np.linspace(0, 2 * np.pi, 100):
                    func_args[0] = psi
                    pp_minus_p2, pp = get_partition_variance(func, func_args, width)
                    pp_minus_p2_psi.append(pp_minus_p2)
                    pp_psi.append(pp)
                pp_minus_p2 = np.mean(pp_minus_p2_psi)
                pp = np.mean(pp)
            else:
                pp_minus_p2, pp = get_partition_variance(func, func_args, width)
            p = width / (2 * np.pi)
            pp_terms.append(pp)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_minus_p2_norm_terms.append(pp_minus_p2 / (p * (1 - p)))
        ax_norm.plot(widths, pp_minus_p2_norm_terms, label=func_name)
        if func_name == 'Combination':
            ax_norm_convo.plot(widths, pp_minus_p2_norm_terms, color='green', label=func_name)
        ax.plot(widths, pp_minus_p2_terms, label=func_name)
        ax_pp.plot(widths, pp_terms, label=func_name)
        pp_minus_p2_dict.update({func_name: pp_minus_p2_terms})
        # pp_minus_p2_dict.update({func_name: pp_terms})
    ax_norm.set_xlabel('Partition Width (w)')
    ax_norm.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    ax_pp.legend()
    fig_norm.tight_layout()
    fig_norm_convo.tight_layout()
    fig.tight_layout()
    fig_pp.tight_layout()

    fig_v2_div_eff, ax_v2_div_eff = plt.subplots(dpi=144, figsize=(7, 3))
    ax_v2_div_eff.axhline(0, color='black')
    for func_name, y in pp_minus_p2_dict.items():
        ax_v2_div_eff.plot(widths, y, label=func_name)
    combo_div_eff = np.array(pp_minus_p2_dict['Combination']) - np.array(pp_minus_p2_dict['Efficiency'])
    ax_v2_div_eff.plot(widths, combo_div_eff, ls='--', label='Combo - Efficiency')
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax_v2_div_eff.legend(loc='upper left')
    fig_v2_div_eff.tight_layout()

    fig_ft, ax_ft = plt.subplots(dpi=144)
    yf_dict, yf_plt_dict = {}, {}
    for func_name, y in pp_minus_p2_dict.items():
        # print(f'{func_name}:\n{y}')
        y = np.nan_to_num(y)
        # print(f'{func_name}:\n{y}')
        n = len(y)
        yf = fft(y)
        xf = fftfreq(n, 2 * np.pi / n)[:n//2]
        ax_ft.plot(xf, 2.0 / n * np.abs(yf[0:n//2]), label=func_name)
        yf_dict.update({func_name: yf})
        yf_plt_dict.update({func_name: 2.0 / n * np.abs(yf[0:n//2])})
    # ax_ft.plot(xf, 2.0 / n * np.abs((yf_dict['Convolution'] / yf_dict['Elliptic Flow'])[0:n//2]),
    #            label='Convolution / Flow')
    # ax_ft.plot(xf, yf_plt_dict['Convolution'] / yf_plt_dict['Elliptic Flow'], label='Convolution / Flow')
    ax_ft.legend()
    ax_ft.grid()
    ax_ft.set_xlabel('Frequency')
    ax_ft.set_ylabel('Power')
    fig_ft.tight_layout()

    plt.show()


def gaus_v2_combo():
    mu = np.pi
    sigma = 0.8
    amp = 0.4
    func1 = base_gaus_pdf_wrap
    func1_args = [mu, sigma, amp, 1. / get_norm_scipy(func1, (mu, sigma, amp, 1))]
    func1_name = 'Gaussian Cluster'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = [psi, v2, n]
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, psi_1_, psi_2_, sigma_, amp_, c_, v2_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf_wrap(x, psi_1_, sigma_, amp_, c_) * vn_pdf(x, psi_2_, v2_, n_)
    func3_args = [mu, psi, sigma, amp, func1_args[-1], v2, n, 1]
    func3_args[-1] = 1. / get_norm_scipy(func3, func3_args)
    func3_name = 'Combination'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig_pdfs, ax_pdfs = plt.subplots(dpi=144, figsize=(8, 4))
    xs = np.linspace(0, 2 * np.pi, 1000)
    ax_pdfs.plot(xs, func1(xs, *func1_args), label=func1_name)
    ax_pdfs.plot(xs, func2(xs, *func2_args), label=func2_name)
    # ax_pdfs.plot(xs, func3(xs, *func3_args), label=func3_name)
    ax_pdfs.set_xlabel(r'$\phi$')
    ax_pdfs.set_ylabel('Probability')
    ax_pdfs.set_ylim(bottom=0)
    ax_pdfs.legend()
    fig_pdfs.tight_layout()

    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    dsigma_dict = {}
    for func, func_args, func_name in \
            [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]:
        dsigma_terms = []
        for width in widths:
            print(f'{func_name}: width = {width}')
            if func_name == 'Combination':
                pass
                psi_bounds = [[0, 2 * np.pi], [0, 2 * np.pi]]
                dsigma = nquad(integrate_partition, psi_bounds, args=((func, width, func_args),))[0]
            else:
                dsigma = quad(integrate_partition, 0, 2 * np.pi, args=((func, width, func_args),))[0]
            dsigma_terms.append(dsigma)
        ax.plot(widths, dsigma_terms, label=func_name)
        dsigma_dict.update({func_name: dsigma_terms})
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax.legend()
    fig.tight_layout()

    fig_v2_div_eff, ax_v2_div_eff = plt.subplots(dpi=144, figsize=(8, 3))
    ax_v2_div_eff.axhline(0, color='black')
    for func_name, y in dsigma_dict.items():
        ax_v2_div_eff.plot(widths, y, label=func_name)
    combo_div_eff = np.array(dsigma_dict['Combination']) - np.array(dsigma_dict[func2_name])
    # combo_div_err = np.array(dsigma_dict['Combination'])**1.5 + np.mean(dsigma_dict[func2_name])**1.5
    combo_div_err = (np.array(dsigma_dict['Combination']) * np.mean(dsigma_dict[func2_name]))**0.75
    # combo_div_err_theory = np.array(dsigma_dict[func1_name])**1.5 + np.array(dsigma_dict[func2_name])**1.5
    combo_div_err_theory = combo_div_err.copy()
    combo_div_err *= 2
    combo_div_err_theory *= 0
    ax_v2_div_eff.plot(widths, combo_div_eff, ls='--', label='Combo - Elliptic Flow', color='red')
    ax_v2_div_eff.fill_between(widths, combo_div_eff - combo_div_err, combo_div_eff + combo_div_err, color='red',
                               alpha=0.3)
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax_v2_div_eff.legend(loc='upper left')
    fig_v2_div_eff.tight_layout()

    fig_corr_dev, ax_corr_dev = plt.subplots(dpi=144, figsize=(8, 3))
    ax_corr_dev.axhline(0, color='black')
    ax_corr_dev.plot(widths, np.array(dsigma_dict[func1_name]) - combo_div_eff)
    cor_diff = np.array(dsigma_dict[func1_name]) - combo_div_eff
    ax_corr_dev.fill_between(widths, cor_diff - combo_div_err, cor_diff + combo_div_err, color='blue', alpha=0.3)
    ax_corr_dev.fill_between(widths, cor_diff - combo_div_err_theory, cor_diff + combo_div_err_theory, color='green',
                             alpha=0.3)
    ax_corr_dev.set_title('Deviation of Correction from True')
    ax_corr_dev.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax_corr_dev.set_xlabel('Partition Width (w)')
    fig_corr_dev.tight_layout()

    plt.show()


def eff_plotting():
    func1_pre = get_efficiency_pdf(62)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x: func1_pre(x) / func1_norm
    func1_args = ()
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    # ax.set_ylim(bottom=0)
    # fig.tight_layout()

    plot_az_bin_example(func1, func1_args, np.pi / 4, np.pi / 4 + np.deg2rad(120))
    ax.set_title('Efficiency Probability Distribution')
    fig.tight_layout()

    plt.show()


def v2_plotting():
    func2 = vn_pdf
    n = 2
    v2 = 0.1
    psi = 0
    func2_args = (psi, v2, n)
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))

    plot_az_bin_example(func2, func2_args, np.pi / 4, np.pi / 4 + np.deg2rad(100))
    ax.set_title('Elliptic Flow Azimuthal Probability Distribution')
    fig.tight_layout()

    plt.show()


def vn_analytic_plotting():
    v_mag = 0.1
    psi = np.pi / 3
    ns = [2, 3, 4, 5]

    fig_prob, ax_prob = plt.subplots(dpi=144, figsize=(8, 3))
    fig_dsig, ax_dsig = plt.subplots(dpi=144, figsize=(8, 3))
    fig_dsig_v2_comp, ax_dsig_v2_comp = plt.subplots(dpi=144, figsize=(8, 3))
    w = np.linspace(0, 2 * np.pi, 1000)

    for n in ns:
        ax_prob.plot(w, vn_pdf(w, psi, v_mag, n), label=fr'$v_{n}$')  # w in place of phi, same range though different
        ax_dsig.plot(w, vn_divs(w, v_mag, n=n), label=fr'$v_{n}$')

    dsigs = []
    for w_i in w:
        ep_diff, ep2 = get_partition_variance(vn_pdf, (psi, v_mag, 2), w_i)
        dsigs.append(ep_diff)
    ax_dsig_v2_comp.plot(w, dsigs, label='Numerical')
    ax_dsig_v2_comp.plot(w, vn_divs(w, v_mag, n=2), ls='--', label='Analytic')

    for ax in [ax_prob, ax_dsig, ax_dsig_v2_comp]:
        ax.legend()
        ax.set_ylim(bottom=0)

    ax_prob.set_title('Probability Densities')
    ax_dsig.set_title(r'$\Delta\sigma^2$')
    ax_dsig_v2_comp.set_title(r'$\Delta\sigma^2$ Analytic vs Numerical')
    ax_prob.set_ylabel('Probability')
    ax_dsig.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax_dsig_v2_comp.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2$')
    ax_prob.set_xlabel(r'$\phi$')
    ax_dsig.set_xlabel(r'Azimuthal Partition Width $w$')
    ax_dsig_v2_comp.set_xlabel(r'Azimuthal Partition Width $w$')
    fig_prob.tight_layout()
    fig_dsig.tight_layout()
    fig_dsig_v2_comp.tight_layout()

    plt.show()


def plot_pdf(func, pars, xs=None):
    if xs is None:
        xs = np.linspace(0, 2 * np.pi, 1000)
    plt.plot(xs, func(xs, *pars))
    plt.xlabel(r'$\phi$')
    plt.ylabel('Probability')


def plot_az_bin_example(func, pars, bin_low, bin_high):
    xs = np.linspace(0, 2 * np.pi, 1000)
    # plt.axhline(0, color='black')
    plot_pdf(func, pars, xs=xs)
    xs_bin = np.linspace(bin_low, bin_high, 1000)
    plt.fill_between(xs_bin, func(xs_bin, *pars), color='gray')
    func_max = max(func(xs_bin, *pars))
    plt.text(bin_low, -0.1 * func_max, r'$\psi$', ha='center')
    plt.text(bin_high, -0.1 * func_max, r'$\psi+w$', ha='center')
    plt.xlim(0, 2 * np.pi)
    plt.ylim(bottom=0)
    plt.tight_layout()


def plot_variance(func, func_args):
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func, func_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    plot_az_bin_example(func, func_args, 2, 2 + np.pi / 3)
    widths = np.linspace(0, 2 * np.pi, 100)
    pp_terms, pp_minus_p2_terms = [], []
    for width in widths:
        pp_minus_p2, pp = get_partition_variance(func, func_args, width)
        p = width / (2 * np.pi)
        pp_terms.append(pp)
        pp_minus_p2_terms.append(pp_minus_p2 / (p * (1 - p)))
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    ax.plot(widths, pp_terms)
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3.5))
    ax.axhline(0, color='black')
    ax.plot(widths, pp_minus_p2_terms)
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    fig.tight_layout()

    # var = get_partition_variance(func, func_args, width)
    plt.show()


# def get_partition_variance_scipy(width, p, p_args):
#     # norm = quad(p, 0, 2 * np.pi, args=p_args)
#     # Need to deal with wrapping non-periodic functions!
#     # def p2(x,  vn, psi, n):
#     #     pass
#     # p2 = lambda x, vn, psi, n: p(x, vn, psi, n)**2
#     # pq = lambda x, vn, psi, n: p(x, vn, psi, n) * (1 - p(x, vn, psi, n))
#     ep = quad(int_pdf, 0, 2 * np.pi, args=(width, p, p_args))[0] / (2 * np.pi)
#     ep2 = quad(int_pdf2, 0, 2 * np.pi, args=(width, p, p_args))[0] / (2 * np.pi)
#     # epq = quad(int_pdf, 0, 2 * np.pi, args=(width, pq, p_args))
#
#     return ep2
#     # return ep2[0] - ep[0]**2
#     # print(quad(int_pdf, 0, 2 * np.pi, args=(width, p, p_args)))
#
#
# def int_pdf(low_bound, width, func, func_args):
#     # print(func_args)
#     return quad(func, low_bound, low_bound + width, args=func_args)[0]
#
#
# def int_pdf2(low_bound, width, func, func_args):
#     # print(func_args)
#     # print(f'low= {low_bound}, width= {width}, func: {func}, func_args: {func_args}')
#     return quad(func, low_bound, low_bound + width, args=func_args)[0] ** 2


def get_partition_variance(func, func_args, width, bounds=(0, 2 * np.pi), points=1000):
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx
    pdf /= pdf_norm
    pdf_wrap = np.append(pdf, np.append(pdf[1:], pdf[1]))  # To represent a periodic boundary glue duplicate array to end

    width_points = round(width / dx)
    probs = [np.sum(pdf_wrap[:width_points])]
    mids = [width / 2]
    for index in range(points):
        probs.append(probs[-1] - pdf_wrap[index] + pdf_wrap[index + width_points])
        mids.append(mids[-1] + dx)
    probs = np.array(probs) * width / width_points
    probs2 = probs**2
    ep = np.mean(probs)
    ep2 = np.mean(probs2)

    return ep2 - ep**2, ep2


def get_partition_variance_scipy(func, func_args, width, bounds=(0, 2 * np.pi)):
    """
    Assumes a normalized input pdf func. Numerically integrates to calculate delta_sigma^2.
    :param func: PDF function, must be normalized
    :param func_args: Arguments to be passed to func
    :param width: Width of azimuthal partition in radians
    :param bounds: Bounds of the psi integration
    :return:
    """
    def func_square(psi, func, width, func_args):
        return quad(func, psi, psi + width, args=func_args)[0] ** 2

    e_p2 = quad(func_square, *bounds, args=(func, width, func_args))[0] / (2 * np.pi)
    e_p = width / (2 * np.pi)

    return e_p2 - e_p**2, e_p2


def integrate_partition(*args):
    phis, func, width, func_args = args[:-1], *args[-1]  # Split input into n phi integral values and function arguments
    func_args[:len(phis)] = phis  # Replace phis in func_args with current integration values
    if len(phis) > 1:
        func_args[-1] = 1  # Reset the function normalization to 1
        func_args[-1] = 1. / get_norm_scipy(func, func_args)  # Recalculate proper normalization with current phis
    p2 = (width / (2 * np.pi))**2
    return (quad(func, 0, width, args=tuple(func_args))[0] ** 2 - p2) / (2 * np.pi) ** len(phis)


def get_norm(func, func_args):
    points = 1000
    bounds = (0, 2 * np.pi)
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    # pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx

    return pdf_norm


def get_norm_scipy(func, func_args, bounds=(0, 2 * np.pi)):
    return quad(func, *bounds, args=tuple(func_args))[0]


def get_partitions_covariance(func, func_pars, width, sep):
    points = 1000
    bounds = (0, 2 * np.pi)
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_pars)
    norm = np.sum(func(xs, *func_pars)) * dx
    pdf /= norm
    pdf_wrap = np.append([pdf, pdf, pdf])  # To represent a periodic boundary glue duplicate array to end

    width_points = round(width / bound_range * points)
    sep_points = round(sep / bound_range * points)
    probs_a = [np.sum(pdf_wrap[:width_points])]
    probs_b = [np.sum(pdf_wrap[sep_points:sep_points + width_points])]
    mids = [width / 2]
    for index in range(points):
        probs_a.append(probs_a[-1] - pdf_wrap[index] + pdf_wrap[index + width_points])
        probs_b.append(probs_b[-1] - pdf_wrap[index + sep_points] + pdf_wrap[index + sep_points + width_points])
        mids.append(mids[-1] + bound_range / points)
    probs_a = np.array(probs_a) * width / width_points
    probs_b = np.array(probs_b) * width / width_points
    # probs_a2 = probs_a**2
    # probs_b2 = probs_b**2
    np.mean((probs_a - probs_b)**2)
    # pqs = probs * (1 - probs)
    # epq = np.mean(pqs)
    # ep = np.mean(probs)
    # ep2 = np.mean(probs2)

    # fig, ax = plt.subplots()
    # ax.plot(np.linspace(bounds[0], 2 * bounds[1], 2 * points), pdf_wrap, color='blue', label='pdf')
    # ax.plot(mids, probs, color='green', label='bin probs')
    # ax.plot(mids, probs2, color='red', label='bin probs squared')
    # ax.axhline(ep, color='green', ls='--', label='expectation of p')
    # ax.axhline(ep**2, color='olive', label='(expectation of p)^2')
    # ax.axhline(ep2, color='red', ls='--', label='expectation of p^2')
    # ax.legend()
    # print(ep2, ep**2)

    # return ep2 - ep**2, epq


def get_efficiency_pdf(energy=62, plot=False):
    set_name = 'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0'
    qa_path = f'F:/Research/Data/default/{set_name}/{energy}GeV/QA_{energy}GeV.root'
    with uproot.open(qa_path) as file:
        hist_name = f'post_phi_{set_name}_{energy}'
        hist = file[hist_name]
        hist_centers, hist_vals = hist.axis().centers(), hist.values()
    interp = get_periodic_interp(hist_centers, hist_vals)
    if plot:
        plt.scatter(hist_centers, hist_vals, color='blue', marker='.')
        xs = np.linspace(0, 2 * np.pi, 10000)
        plt.plot(xs, interp(xs), color='red', alpha=0.6)
        plt.show()
    return interp


def get_periodic_interp(x, y):
    x_step = x[-1] - x[-2]
    x = np.append(x, [x[-1] + x_step])
    y = np.append(y, [y[0]])
    interp = CubicSpline(x, y, bc_type='periodic')
    return interp


def vn_pdf(phi, psi, vn, n=2):
    return (1 + 2 * vn * np.cos(n * (phi - psi))) / (2 * np.pi)


def gaus_pdf(phi, mu, sigma):
    return np.exp(-0.5 * ((phi - mu) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))


def base_gaus_pdf(phi_in, mu, sigma, amp, normalization):
    phi = phi_in % (2 * np.pi)
    return normalization * (1 + amp * np.exp(-0.5 * ((phi - mu) / sigma)**2))


def base_gaus_pdf_wrap(phi_in, mu, sigma, amp, normalization, wrap=3):
    phi = phi_in % (2 * np.pi)
    gaus_pdf = 1 + amp * np.exp(-0.5 * ((phi - mu) / sigma)**2)
    for i in range(1, wrap + 1):
        phi_right = phi + i * 2 * np.pi
        gaus_pdf += amp * np.exp(-0.5 * ((phi_right - mu) / sigma)**2)
        phi_left = phi - i * 2 * np.pi
        gaus_pdf += amp * np.exp(-0.5 * ((phi_left - mu) / sigma)**2)

    # return gaus_pdf / get_norm(base_gaus_pdf, (mu, sigma, amp, 1))
    return normalization * gaus_pdf


# def quad_180_sin(x, a, c, )
def sin(x, a, f):
    return a * (1 + np.sin(2 * np.pi * f * x))


def cos(x, a, f):
    return a * (1 + np.cos(2 * np.pi * f * x))


def cos_pi(x, a, f):
    return a * (1 + np.cos(2 * np.pi * f * x + np.pi))


def cos_pi_fixed(x, a):
    return a * (1 - np.cos(x))


def quad_cos(x, a1, c, a2, sigma):
    # gaus_mod = np.exp(((x - np.pi) / sigma)**2) / np.sqrt(2 * np.pi * sigma)
    gaus_mod = np.exp(-((x - np.pi) / sigma)**2)
    quad_mod = quad_180_rad(x, a1, c) * gaus_mod
    # quad_mod = 0
    cos_mod = cos_pi_fixed(x, a2) * (1 - gaus_mod)
    # cos_mod = 0

    return quad_mod + cos_mod


def quad_gaus(x, a, c, c2, sigma):
    return quad_180(x, a, c) * np.exp(((x - np.pi) / sigma)**2) + c2


def cos_gaus(x, a, sigma):
    return a * (1 - np.cos(x)) / np.exp(-((x - np.pi) / sigma)**2)


def cos_sin_gaus(x, a, b, sigma):
    return a * (1 - np.cos(x) + b * np.sin(3.5 * x)) / np.exp(-((x - np.pi) / sigma)**2)


def cos_sin(x, a, b):
    return a * (1 - np.cos(x) + b * np.sin(1.5 * x))


def cos_sin4(x, a, b, c, d):
    return a * (1 - np.cos(x) + b * np.sin(1.5 * x) + c * np.sin(2.5 * x) + d * np.sin(3.5 * x))


if __name__ == '__main__':
    main()
