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
from scipy.integrate import quad
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
    gaus_fit_dependence()
    # convo_test()
    # eff_v2_combo()
    # eff_v2_combo2()
    # eff_gaus_combo()
    # eff_gaus_combo2()
    # gaus_v2_combo2()
    print('donzo')


def vn_test():
    func = vn_pdf
    n = 2
    v = 0.07
    psi = np.pi / 3
    func_args = (v, psi, n)

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
    sigma = 0.8
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
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')

    width_fit_low, width_fit_high = np.deg2rad([60, 300])
    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
    ax.legend()
    fig.tight_layout()

    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, pp_minus_p2_terms)
    ax.set_ylabel('Linear Terms')
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
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
    ax.legend()
    fig.tight_layout()

    plt.show()


def gaus_fit_dependence():
    mu = np.pi
    sigmas = np.linspace(0.3, 1.0, 10)
    amp = 0.5
    bounds = (0, 2 * np.pi)
    func = base_gaus_pdf
    widths = np.linspace(*bounds, 100)
    xs_fit_plot = np.linspace(min(widths), max(widths), 1000)

    quad_mags, quad_consts, cos_mags, gaus_sigs = [], [], [], []
    quad_mags_err, quad_consts_err, cos_mags_err, gaus_sigs_err = [], [], [], []
    for sigma in sigmas:
        func_args = (mu, sigma, amp, 1. / get_norm(func, (mu, sigma, amp, 1)))

        pp_minus_p2_terms, pp_terms = [], []
        for width in widths:
            pp_minus_p2, pp = get_partition_variance(func, func_args, width, bounds)
            pp_minus_p2_terms.append(pp_minus_p2)
            pp_terms.append(pp)

        p0 = (-3e-4, 0.002, 0.0004, 1)
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
        ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
        ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)
        ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *p0), color='olive', alpha=0.7, ls='--', label='Combo Initial Guess')
        ax.plot(xs_fit_plot, quad_cos(xs_fit_plot, *popt), color='olive', label='Combo')
        ax.set_title(f'Fit for sigma={sigma:.2f}')
        ax.legend()
        fig.tight_layout()

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
    func2_args = (v2, psi, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, mu_, sigma_, amp_, c_, v2_, psi_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf(x, mu_, sigma_, amp_, c_) * vn_pdf(x, v2_, psi_, n_)
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    func2_args = (v2, psi, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, v2_, psi_, n_, c_combo_: c_combo_ * func1(x) * vn_pdf(x, v2_, psi_, n_)
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    ax_v2_div_eff.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
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
    func2_args = (v2, psi, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, v2_, psi_, n_, c_combo_: c_combo_ * func1(x) * vn_pdf(x, v2_, psi_, n_)
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    ax_v2_div_eff.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    ax_v2_div_eff.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    ax_v2_div_eff.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
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


def gaus_v2_combo2():
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
    func2_args = (v2, psi, n)
    func2_name = 'Elliptic Flow'
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, mu_, sigma_, amp_, c_, v2_, psi_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf(x, mu_, sigma_, amp_, c_) * vn_pdf(x, v2_, psi_, n_)
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
                    func_args[5] = psi
                    # print(func_args)
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
    ax_norm.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm.legend()
    ax_norm_convo.set_xlabel('Partition Width (w)')
    ax_norm_convo.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    ax_norm_convo.legend()
    ax.set_xlabel('Partition Width (w)')
    ax.set_ylabel(r'$\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2$')
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
    combo_div_eff = np.array(pp_minus_p2_dict['Combination']) - np.array(pp_minus_p2_dict['Elliptic Flow'])
    ax_v2_div_eff.plot(widths, combo_div_eff, label='Combo - Ellpitic Flow')
    ax_v2_div_eff.set_xlabel('Partition Width (w)')
    ax_v2_div_eff.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
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


def plot_pdf(func, pars, xs=None):
    if xs is None:
        xs = np.linspace(0, 2 * np.pi, 1000)
    plt.plot(xs, func(xs, *pars))
    plt.xlabel('phi')
    plt.ylabel('Probability')


def plot_az_bin_example(func, pars, bin_low, bin_high):
    xs = np.linspace(0, 2 * np.pi, 1000)
    # plt.axhline(0, color='black')
    plot_pdf(func, pars, xs=xs)
    xs_bin = np.linspace(bin_low, bin_high, 1000)
    plt.fill_between(xs_bin, func(xs_bin, *pars), color='gray')
    plt.xlim(0, 2 * np.pi)
    plt.ylim(bottom=0)


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
    ax.set_ylabel(r'$\left[\int_{0}^{2\pi}p(\psi)^2 \,d\phi - p^2\right] / \left[p (1-p)\right]$')
    fig.tight_layout()

    # var = get_partition_variance(func, func_args, width)
    plt.show()


def get_partition_variance_scipy(width, p, p_args):
    # norm = quad(p, 0, 2 * np.pi, args=p_args)
    # Need to deal with wrapping non-periodic functions!
    # def p2(x,  vn, psi, n):
    #     pass
    # p2 = lambda x, vn, psi, n: p(x, vn, psi, n)**2
    # pq = lambda x, vn, psi, n: p(x, vn, psi, n) * (1 - p(x, vn, psi, n))
    ep = quad(int_pdf, 0, 2 * np.pi, args=(width, p, p_args))[0] / (2 * np.pi)
    ep2 = quad(int_pdf2, 0, 2 * np.pi, args=(width, p, p_args))[0] / (2 * np.pi)
    # epq = quad(int_pdf, 0, 2 * np.pi, args=(width, pq, p_args))

    return ep2
    # return ep2[0] - ep[0]**2
    # print(quad(int_pdf, 0, 2 * np.pi, args=(width, p, p_args)))


def int_pdf(low_bound, width, func, func_args):
    # print(func_args)
    return quad(func, low_bound, low_bound + width, args=func_args)[0]


def int_pdf2(low_bound, width, func, func_args):
    # print(func_args)
    # print(f'low= {low_bound}, width= {width}, func: {func}, func_args: {func_args}')
    return quad(func, low_bound, low_bound + width, args=func_args)[0]**2


def get_partition_variance(func, func_args, width, bounds=(0, 2 * np.pi)):
    points = 1000
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx
    pdf /= pdf_norm
    pdf_wrap = np.append(pdf, pdf)  # To represent a periodic boundary glue duplicate array to end

    width_points = round(width / bound_range * points)
    probs = [np.sum(pdf_wrap[:width_points])]
    mids = [width / 2]
    for index in range(points):
        probs.append(probs[-1] - pdf_wrap[index] + pdf_wrap[index + width_points])
        mids.append(mids[-1] + bound_range / points)
    probs = np.array(probs) * width / width_points
    probs2 = probs**2
    ep = np.mean(probs)
    ep2 = np.mean(probs2)

    # fig, ax = plt.subplots()
    # ax.plot(np.linspace(bounds[0], 2 * bounds[1], 2 * points), pdf_wrap, color='blue', label='pdf')
    # ax.plot(mids, probs, color='green', label='bin probs')
    # ax.plot(mids, probs2, color='red', label='bin probs squared')
    # ax.axhline(ep, color='green', ls='--', label='expectation of p')
    # ax.axhline(ep**2, color='olive', label='(expectation of p)^2')
    # ax.axhline(ep2, color='red', ls='--', label='expectation of p^2')
    # ax.legend()
    # print(ep2, ep**2)

    return ep2 - ep**2, ep2


def get_partition_variance(func, func_args, width, bounds=(0, 2 * np.pi)):
    points = 1000
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx
    pdf /= pdf_norm
    pdf_wrap = np.append(pdf, pdf)  # To represent a periodic boundary glue duplicate array to end

    width_points = round(width / bound_range * points)
    probs = [np.sum(pdf_wrap[:width_points])]
    mids = [width / 2]
    for index in range(points):
        probs.append(probs[-1] - pdf_wrap[index] + pdf_wrap[index + width_points])
        mids.append(mids[-1] + bound_range / points)
    probs = np.array(probs) * width / width_points
    probs2 = probs**2
    ep = np.mean(probs)
    ep2 = np.mean(probs2)

    # fig, ax = plt.subplots()
    # ax.plot(np.linspace(bounds[0], 2 * bounds[1], 2 * points), pdf_wrap, color='blue', label='pdf')
    # ax.plot(mids, probs, color='green', label='bin probs')
    # ax.plot(mids, probs2, color='red', label='bin probs squared')
    # ax.axhline(ep, color='green', ls='--', label='expectation of p')
    # ax.axhline(ep**2, color='olive', label='(expectation of p)^2')
    # ax.axhline(ep2, color='red', ls='--', label='expectation of p^2')
    # ax.legend()
    # print(ep2, ep**2)

    return ep2 - ep**2, ep2


def get_norm(func, func_args):
    points = 1000
    bounds = (0, 2 * np.pi)
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    # pdf = func(xs, *func_args)
    pdf_norm = np.sum(func(xs, *func_args)) * dx

    return pdf_norm


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
    set_name = 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_0'
    qa_path = f'F:/Research/Data/default_resample/{set_name}/{energy}GeV/QA_{energy}GeV.root'
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


def vn_pdf(phi, vn, psi, n=2):
    return (1 + 2 * vn * np.cos(n * (phi - psi))) / (2 * np.pi)


def gaus_pdf(phi, mu, sigma):
    return np.exp(-0.5 * ((phi - mu) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))


def base_gaus_pdf(phi, mu, sigma, amp, normalization):
    return normalization * (1 + amp * np.exp(-0.5 * ((phi - mu) / sigma)**2))


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


if __name__ == '__main__':
    main()
