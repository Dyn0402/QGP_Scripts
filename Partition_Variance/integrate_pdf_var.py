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
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mtick
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
    # eff_v2_combo()
    eff_gaus_combo()
    # gaus_v2_combo()
    # eff_plotting()
    # v2_plotting()
    # vn_analytic_plotting()
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
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')

    width_fit_low, width_fit_high = np.deg2rad([60, 300])
    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
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
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
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
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
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
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
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

    ax.set_xlabel('Azimuthal Partition Width (w)')
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
        ax.set_xlabel('Azimuthal Partition Width (w)')
        ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
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


def eff_v2_combo():
    func1_pre, func1_points = get_efficiency_pdf(62, return_centers=True)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x, psi: func1_pre(x - psi) / func1_norm
    psi_eff = 0
    func1_args = [psi_eff, ]
    func1_name = 'Efficiency'
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

    func3 = lambda x, psi_eff_, psi_, v2_, n_, c_combo_: c_combo_ * func1(x, psi_eff_) * vn_pdf(x, psi_, v2_, n_)
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

    func_list = [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]
    func_points = {func1_name: func1_points, func3_name: func1_points}
    widths, dsigma_dict = integrate_funcs(func_list, func_points=func_points)
    subtract_funcs(widths, dsigma_dict, func3_name, func2_name, func1_name)
    plt.show()


def eff_gaus_combo():
    func1_pre, func1_points = get_efficiency_pdf(62, return_centers=True)
    func1_norm = get_norm(func1_pre, ())
    func1 = lambda x, psi: func1_pre(x - psi) / func1_norm
    psi_eff = 0
    func1_args = [psi_eff, ]
    func1_name = 'Efficiency'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    mu = np.pi
    sigma = 0.8
    p = 0.3
    amp = np.sqrt(2 * np.pi) / sigma * p / (1 - p)
    func2 = base_gaus_pdf_wrap
    func2_args = [mu, sigma, amp, 1. / get_norm_scipy(func2, (mu, sigma, amp, 1))]
    func2_name = 'Gaussian Cluster'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, psi_, mu_, sigma_, amp_, c_, c_combo_: c_combo_ * func1(x, psi_) * \
                                                         base_gaus_pdf_wrap(x, mu_, sigma_, amp_, c_)
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
    ax_pdfs.plot(xs, func3(xs, *func3_args), label=func3_name)
    ax_pdfs.set_xlabel(r'$\phi$')
    ax_pdfs.set_ylabel('Probability')
    ax_pdfs.set_ylim(bottom=0)
    ax_pdfs.legend()
    fig_pdfs.tight_layout()

    func_list = [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]
    func_points = {func1_name: func1_points, func3_name: func1_points}
    widths, dsigma_dict = integrate_funcs(func_list, func_points=func_points)
    subtract_funcs(widths, dsigma_dict, func3_name, func1_name, func2_name)
    plt.show()


def gaus_v2_combo():
    mu = np.pi
    sigma = 0.8
    p = 0.3
    amp = np.sqrt(2 * np.pi) / sigma * p / (1 - p)
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
    ax_pdfs.plot(xs, func3(xs, *func3_args), label=func3_name)
    ax_pdfs.set_xlabel(r'$\phi$')
    ax_pdfs.set_ylabel('Probability')
    ax_pdfs.set_ylim(bottom=0)
    ax_pdfs.legend()
    fig_pdfs.tight_layout()

    func_list = [(func1, func1_args, func1_name), (func2, func2_args, func2_name), (func3, func3_args, func3_name)]
    widths, dsigma_dict = integrate_funcs(func_list)
    subtract_funcs(widths, dsigma_dict, func3_name, func2_name, func1_name)
    plt.show()


def integrate_funcs(func_list, func_points=None):
    fig, ax = plt.subplots(dpi=144, figsize=(7, 3))
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    dsigma_dict = {}
    for func, func_args, func_name in func_list:
        dsigma_terms = []
        for width in widths:
            print(f'{func_name}: width = {width}')

            if func_name == 'Combination':
                psi_bounds = [[0, 2 * np.pi], [0, 2 * np.pi]]
            else:
                psi_bounds = [[0, 2 * np.pi]]

            points = func_points[func_name] if func_points is not None and func_name in func_points else None
            dsigma_terms.append(calc_dsigma(func, func_args, width, psi_bounds, points))

        ax.plot(widths, dsigma_terms, label=func_name)
        dsigma_dict.update({func_name: dsigma_terms})
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax.legend()
    fig.tight_layout()

    return widths, dsigma_dict


def subtract_funcs(widths, dsigma_dict, combo_name='Combination', name_to_subtract=None, base_name=None):
    names = [x for x in dsigma_dict.keys() if x != combo_name]
    if name_to_subtract is None:
        name_to_subtract = names.pop(0)
    if base_name is None:
        base_name = names.pop(0)

    fig = plt.figure(dpi=144, figsize=(8, 5))
    gs = GridSpec(2, 1, height_ratios=[0.7, 0.3])

    # fig_v2_div_eff, ax_integral = plt.subplots(dpi=144, figsize=(8, 3))
    ax_integral = fig.add_subplot(gs[0])
    ax_integral.axhline(0, color='black')
    for func_name, y in dsigma_dict.items():
        ax_integral.plot(widths, y, label=func_name)
    combo_div_eff = np.array(dsigma_dict[combo_name]) - np.array(dsigma_dict[name_to_subtract])
    # combo_div_err = (np.array(dsigma_dict[combo_name]) * np.mean(dsigma_dict[name_to_subtract]))**0.75
    # combo_div_err = np.array(dsigma_dict[combo_name]) * np.sqrt(np.nanmean(dsigma_dict[name_to_subtract]))
    combo_div_err = (np.array(dsigma_dict[base_name]) * np.nanmean(dsigma_dict[name_to_subtract]))**0.75
    # combo_div_err = np.array(dsigma_dict[combo_name]) * np.sqrt(dsigma_dict[name_to_subtract])
    # combo_div_err /= 2
    diff_name = f'{combo_name} - {name_to_subtract}'
    ax_integral.plot(widths, combo_div_eff, ls='--', label=diff_name, color='red')
    ax_integral.fill_between(widths, combo_div_eff - combo_div_err, combo_div_eff + combo_div_err, color='red',
                               alpha=0.3)
    # ax_integral.set_xlabel('Azimuthal Partition Width (w)')
    ax_integral.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax_integral.legend(loc='upper left')
    ax_integral.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    # fig_v2_div_eff.tight_layout()

    # fig_corr_dev, ax_corr_dev = plt.subplots(dpi=144, figsize=(8, 3))
    ax_corr_dev = fig.add_subplot(gs[1], sharex=ax_integral)
    ax_corr_dev.axhline(0, color='black')
    ax_corr_dev.plot(widths, np.array(dsigma_dict[base_name]) - combo_div_eff, color='red',
                     label=f'{base_name} - ({diff_name})')
    cor_diff = np.array(dsigma_dict[base_name]) - combo_div_eff
    ax_corr_dev.fill_between(widths, cor_diff - combo_div_err, cor_diff + combo_div_err, color='red', alpha=0.3)
    # ax_corr_dev.set_title('Deviation of Correction from True')
    ax_corr_dev.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax_corr_dev.set_xlabel('Azimuthal Partition Width (w)')
    ax_corr_dev.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax_corr_dev.legend(loc='lower left')
    # fig_corr_dev.tight_layout()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.025)


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
    fig = plt.figure(dpi=144, figsize=(8, 5))
    gs = GridSpec(2, 1, height_ratios=[0.7, 0.3])
    fig_dsig, ax_dsig = plt.subplots(dpi=144, figsize=(8, 3))
    # fig_dsig_v2_comp, ax_dsig_v2_comp = plt.subplots(dpi=144, figsize=(8, 3))
    ax_dsig_v2_comp = fig.add_subplot(gs[0])
    ax_dsig_v2_comp_dev = fig.add_subplot(gs[1], sharex=ax_dsig_v2_comp)
    ax_dsig_v2_comp_dev.axhline(0, color='black', lw=0.5)
    w = np.linspace(0, 2 * np.pi, 1000)

    for n in ns:
        ax_prob.plot(w, vn_pdf(w, psi, v_mag, n), label=fr'$v_{n}$')  # w in place of phi, same range though different
        ax_dsig.plot(w, vn_divs(w, v_mag, n=n), label=fr'$v_{n}$')

    dsigs = []
    for w_i in w:
        dsig = nquad(integrate_partition, [[0, 2 * np.pi]], args=((vn_pdf, w_i, [psi, v_mag, 2], None),))[0]
        dsigs.append(dsig)
    ax_dsig_v2_comp.plot(w, dsigs, label='Numerical')
    ax_dsig_v2_comp.plot(w, vn_divs(w, v_mag, n=2), ls='--', label='Analytic')
    ax_dsig_v2_comp_dev.plot(w, dsigs - vn_divs(w, v_mag, n=2), label='Numerical - Analytic')

    for ax in [ax_prob, ax_dsig, ax_dsig_v2_comp]:
        ax.legend()
        ax.set_ylim(bottom=0)
    ax_dsig_v2_comp_dev.legend()

    ax_prob.set_title('Probability Densities')
    ax_dsig.set_title(r'$\Delta\sigma^2$')
    ax_dsig_v2_comp.set_title(r'$\Delta\sigma^2$ Analytic vs Numerical')
    ax_prob.set_ylabel('Probability')
    ax_dsig.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax_dsig_v2_comp.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax_dsig_v2_comp_dev.set_ylabel(r'$\frac{1}{2\pi}\int_{0}^{2\pi} \delta p(\psi)^2 \,d\psi$')
    ax_prob.set_xlabel(r'$\phi$')
    ax_dsig.set_xlabel(r'Azimuthal Partition Width $w$')
    # ax_dsig_v2_comp.set_xlabel(r'Azimuthal Partition Width $w$')
    ax_dsig_v2_comp_dev.set_xlabel(r'Azimuthal Partition Width $w$')
    ax_dsig_v2_comp.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax_dsig_v2_comp.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_dsig_v2_comp_dev.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    fig_prob.tight_layout()
    fig_dsig.tight_layout()
    # fig_dsig_v2_comp.tight_layout()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)

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
    ax.set_xlabel('Azimuthal Partition Width (w)')
    ax.set_ylabel(r'$\left[\frac{1}{2\pi}\int_{0}^{2\pi}p(\psi)^2 \,d\psi - p^2\right] / \left[p (1-p)\right]$')
    fig.tight_layout()

    # var = get_partition_variance(func, func_args, width)
    plt.show()


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


def calc_dsigma(func, func_args, width, psi_bounds, points):
    if points is None:
        dsigma = nquad(integrate_partition, psi_bounds, args=((func, width, func_args, points),))[0]
    else:
        if len(psi_bounds) == 1:
            dsigma = integrate_partition_delta_prob(func, func_args, width, points)
        else:
            psi_vals = [np.linspace(x[0], x[1], 1000) for x in psi_bounds[:-1]]  # Do last psi average in integral bound
            psi_vecs = np.stack([x.ravel() for x in np.meshgrid(*psi_vals, indexing='ij')], axis=-1)  # All psi combos
            dsigs = []
            for psi_vec in psi_vecs:
                func_args[:len(psi_vec)] = psi_vec
                dsigs.append(integrate_partition_delta_prob(func, func_args, width, points))
            dsigma = np.mean(dsigs)
    return dsigma


def integrate_partition_delta_prob(func, func_args, width, xs):
    dx = xs[1] - xs[0]
    pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx
    pdf /= pdf_norm
    pdf_wrap = np.append(pdf, np.append(pdf, pdf))  # To represent a periodic boundary glue duplicate array to end
    delta_pdf_wrap = pdf_wrap - 1. / (2 * np.pi)

    width_points = round(width / dx)
    delta_probs = [np.sum(delta_pdf_wrap[:width_points])]
    for index in range(xs.size):
        delta_probs.append(delta_probs[-1] - delta_pdf_wrap[index] + delta_pdf_wrap[index + width_points])
    delta_probs = np.array(delta_probs) * width / width_points
    delta_probs2 = delta_probs**2
    edp2 = np.mean(delta_probs2)

    return edp2


def integrate_partition(*args):
    phis, func, width, func_args, pnts = args[:-1], *args[-1]  # Split input into n phi integral values and other

    func_args[:len(phis)] = phis  # Replace phis in func_args with current integration values

    if len(phis) > 1:
        func_args[-1] = 1  # Reset the function normalization to 1
        func_args[-1] = 1. / get_norm_scipy(func, func_args)  # Recalculate proper normalization with current phis
    integral = quad(func, 0, width, args=tuple(func_args))[0]

    p2 = (width / (2 * np.pi)) ** 2
    return (integral ** 2 - p2) / (2 * np.pi) ** len(phis)


def get_norm(func, func_args, bounds=(0, 2 * np.pi), xs=None):
    if xs is None:
        points = 1000
        xs = np.linspace(*bounds, points)
    else:
        points = xs.size

    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
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


def get_efficiency_pdf(energy=62, plot=False, return_centers=False):
    set_name = 'rapid05_resample_norotate_seed_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_calcv2_0'
    qa_path = f'F:/Research/Data/default/{set_name}/{energy}GeV/QA_{energy}GeV.root'
    with uproot.open(qa_path) as file:
        hist_name = f'post_phi_{set_name}_{energy}'
        hist = file[hist_name]
        hist_centers, hist_vals = hist.axis().centers(), hist.values()
    interpolation = get_periodic_interp(hist_centers, hist_vals)
    if plot:
        plt.scatter(hist_centers, hist_vals, color='blue', marker='.')
        xs = np.linspace(0, 2 * np.pi, 10000)
        plt.plot(xs, interpolation(xs), color='red', alpha=0.6)
        plt.show()
    if return_centers:
        return interpolation, hist_centers
    return interpolation


def get_periodic_interp(x, y):
    x_step = x[-1] - x[-2]
    x = np.append(x, [x[-1] + x_step])
    y = np.append(y, [y[0]])
    # interpolation = CubicSpline(x, y, bc_type='periodic', extrapolate='periodic')
    interpolation = lambda xx: np.interp(xx, x, y, period=2 * np.pi)
    return interpolation


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
