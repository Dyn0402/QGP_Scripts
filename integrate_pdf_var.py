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


def main():
    # vn_test()
    gaus_test()
    print('donzo')


def vn_test():
    func = vn_pdf
    n = 3
    v2 = 0.02
    psi = np.pi / 3
    func_args = (v2, psi, n)
    width = np.pi / 3

    plot_pdf(func, func_args)
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
    sigma = 0.5
    func = norm(mu, sigma).pdf
    func_args = ()

    plot_pdf(func, func_args)
    widths = np.linspace(0, 2 * np.pi, 100)
    lin_terms, const_terms = [], []
    for width in widths:
        lin, const = get_partition_variance(func, func_args, width)
        lin_terms.append(lin)
        const_terms.append(const)
    fig, ax = plt.subplots()
    ax.plot(widths, lin_terms)
    fig, ax = plt.subplots()
    ax.plot(widths, const_terms)
    # var = get_partition_variance(func, func_args, width)
    plt.show()


def plot_pdf(func, pars):
    xs = np.linspace(0, 2 * np.pi, 1000)
    plt.plot(xs, func(xs, *pars))
    plt.xlabel('phi')
    plt.ylabel('Probability')


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


def get_partition_variance(func, func_pars, width):
    points = 1000
    bounds = (0, 2 * np.pi)
    xs = np.linspace(*bounds, points)
    bound_range = bounds[-1] - bounds[0]
    dx = bound_range / points
    pdf = func(xs, *func_pars)
    norm = np.sum(func(xs, *func_pars)) * dx
    pdf /= norm
    pdf_wrap = np.append(pdf, pdf)  # To represent a periodic boundary glue duplicate array to end

    width_points = round(width / bound_range * points)
    probs = [np.sum(pdf_wrap[:width_points])]
    mids = [width / 2]
    for index in range(points):
        probs.append(probs[-1] - pdf_wrap[index] + pdf_wrap[index + width_points])
        mids.append(mids[-1] + bound_range / points)
    probs = np.array(probs) * width / width_points
    probs2 = probs**2
    pqs = probs * (1 - probs)
    epq = np.mean(pqs)
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

    return ep2 - ep**2, epq


def vn_pdf(phi, vn, psi, n=2):
    return (1 + 2 * vn * np.cos(n * (phi - psi))) / (2 * np.pi)


def gaus(x, a, b, c):
    pass


if __name__ == '__main__':
    main()
