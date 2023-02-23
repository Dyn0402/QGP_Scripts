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

from Binom_Slices.analyze_binom_slices import quad_180


def main():
    # vn_test()
    # gaus_test()
    convo_test()
    print('donzo')


def vn_test():
    func = vn_pdf
    n = 3
    v2 = 0.1
    psi = np.pi / 3
    func_args = (v2, psi, n)
    width = np.pi / 3

    # plot_pdf(func, func_args)
    plot_az_bin_example(func, func_args, 2, 2 + np.pi / 3)
    plt.show()
    widths = np.linspace(0, 2 * np.pi, 100)
    lin_terms, const_terms = [], []
    for width in widths:
        # vars_a.append(get_partition_variance_scipy(width, func, func_args))
        lin, const = get_partition_variance(func, func_args, width)
        lin_terms.append(lin)
        const_terms.append(const)
    fig, ax = plt.subplots()
    ax.plot(widths, lin_terms)
    fig, ax = plt.subplots()
    ax.plot(widths, const_terms)
    # var = get_partition_variance(func, func_args, width)
    plt.show()


def gaus_test():
    mu = 3
    sigma = 1.5
    # func = norm(mu, sigma).pdf
    # func_args = ()
    func = base_gaus_pdf
    func_args = (mu, sigma, )
    bounds = (0, 2 * np.pi)

    plot_pdf(func, func_args)
    widths = np.linspace(*bounds, 100)
    lin_terms, const_terms = [], []
    for width in widths:
        lin, const = get_partition_variance(func, func_args, width, bounds)
        lin_terms.append(lin)
        const_terms.append(const)
    fig, ax = plt.subplots()
    ax.axhline(0, color='black')
    ax.plot(widths, lin_terms)
    ax.set_ylabel('Linear Terms')

    width_fit_low, width_fit_high = 60, 300
    fig, ax = plt.subplots(dpi=144)
    ax.axhline(0, color='black')
    ax.scatter(widths, lin_terms)
    ax.set_ylabel('Linear Terms')
    width_filter = (np.deg2rad(width_fit_low) < widths) & (widths < np.deg2rad(width_fit_high))
    # print(widths[width_filter], np.array(lin_terms)[width_filter])
    popt, pcov = cf(quad_180, np.rad2deg(widths[width_filter]), np.array(lin_terms)[width_filter])
    xs_fit_plot = np.linspace(min(widths), max(widths), 1000)
    ax.plot(xs_fit_plot, quad_180(np.rad2deg(xs_fit_plot), *popt), color='red', label='Quadratic')
    ax.axvline(np.deg2rad(width_fit_low), ls='--', color='orange', label='Fit Range')
    ax.axvline(np.deg2rad(width_fit_high), ls='--', color='orange')
    ax.set_ylim(bottom=-ax.get_ylim()[-1] * 0.05)
    fig.tight_layout()

    p0 = (0.075, 0.2)
    popt, pcov = cf(cos_pi, widths[width_filter], np.array(lin_terms)[width_filter], p0=p0)
    # ax.plot(xs_fit_plot, cos_pi(xs_fit_plot, *p0), color='purple', alpha=0.3, label='Sine p0')
    ax.plot(xs_fit_plot, cos_pi(xs_fit_plot, *popt), color='purple', label='Sine')
    # ax.plot(xs_fit_plot, 5 * cos_pi(xs_fit_plot, *popt)**2, color='green', label='Sine**2')

    ax.legend()

    # fig, ax = plt.subplots()
    # ax.axhline(0, color='black')
    # ax.plot(widths, const_terms)
    # ax.set_ylabel('Constant Terms')
    # var = get_partition_variance(func, func_args, width)
    plt.show()


def convo_test():
    mu = np.pi
    sigma = 0.1
    amp = 0.5
    func1 = base_gaus_pdf
    func1_args = (mu, sigma, amp, 1. / get_norm(func1, (mu, sigma, amp, 1)))
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func1, func1_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func2 = vn_pdf
    n = 2
    v2 = 0.07
    psi = np.pi / 3
    func2_args = (v2, psi, n)
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func2, func2_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    func3 = lambda x, mu_, sigma_, base_, c_, v2_, psi_, n_, c_combo_: \
        c_combo_ * base_gaus_pdf(x, mu_, sigma_, base_, c_) * vn_pdf(x, v2_, psi_, n_)
    func3_args = (*func1_args, *func2_args, 1. / get_norm(func3, (*func1_args, *func2_args, 1)))
    fig, ax = plt.subplots(dpi=144, figsize=(6, 3))
    plot_pdf(func3, func3_args)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    ax.axhline(0, color='black')
    widths = np.linspace(0, 2 * np.pi, 100)
    for func, func_args in [(func1, func1_args), (func2, func2_args), (func3, func3_args)]:
        lin_terms = []
        for width in widths:
            full_lin, p2 = get_partition_variance(func, func_args, width)
            lin_terms.append(full_lin)
        ax.plot(widths, lin_terms, label=func.__name__)
    ax.legend()
    fig.tight_layout()
    fig2.tight_layout()
    # var = get_partition_variance(func, func_args, width)
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


def get_norm(func, func_args):
    points = 1000
    bounds = (0, 2 * np.pi)
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_args)
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


if __name__ == '__main__':
    main()
