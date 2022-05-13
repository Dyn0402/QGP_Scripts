#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 11 4:29 PM 2022
Created in PyCharm
Created as QGP_Scripts/sim_comp_fits.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def main():
    spreads = [0.01, 0.1, 0.2, 0.5, 0.6, 0.65, 0.7000000000000001, 0.75, 0.8, 0.8500000000000001, 0.8999999999999999,
               1.0, 1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0]
    amps = [0.14372594165489835, 0.04485752036521288, 0.02735948930506084, 0.015407327801912606, 0.014521546458387337,
            0.013938072798862717, 0.013890054595783833, 0.013863623528661457, 0.013871419116892365,
            0.013845253000711478, 0.014134780681031902, 0.014034857347270549, 0.019120916782309448, 0.0340362628330822,
            0.05001593431142087, 0.08107800985334815, 0.13093587396755724, 0.21827866260369758, 0.33779853405711235,
            0.45057315850250823, 0.5202582884414335, 0.523722398499884, 0.4939609892987397, 0.5167445769330699]
    amps_err = [0.005700400117642654, 0.0015041008040337304, 0.000709766154299732, 0.0001847216516244048,
                0.00020317758108751076, 0.00019914731210134878, 0.0002154410003142111, 0.0002494976402945341,
                0.0002536845555085603, 0.00016804748922418083, 0.0001785812154133484, 0.0002599658163093812,
                0.00033576628874251886, 0.00129648761478393, 0.0010276192095396985, 0.0006352971475365707,
                0.0010299309760495286, 0.0029162449957741395, 0.0022397992830927695, 0.0013355169473621767,
                0.0010299164319222801, 0.002249342194583998, 0.001448085275625134, 0.0001841953180415634]

    x_low_left, x_high_left = -0.15, 1.1
    p0_left = (0.05, 7, 1, 0.013)
    spreads_fit, amps_fit, amps_err_fit = trim(spreads, x_low_left, x_high_left, amps, amps_err)
    popt_left, pcov_left = curve_fit(exp_zero, spreads_fit, amps_fit, p0=p0_left, sigma=amps_err_fit, absolute_sigma=True)
    print(popt_left)

    fig_left, ax_left = plt.subplots()
    fit_xs = np.linspace(x_low_left, x_high_left, 1000)
    ax_left.plot(fit_xs, exp_zero(fit_xs, *p0_left), ls='--', color='gray', label='Initial Guess')
    ax_left.plot(fit_xs, exp_zero(fit_xs, *popt_left), color='red', label='Fit')
    ax_left.errorbar(spreads_fit, amps_fit, yerr=amps_err_fit, marker='o', ls='none')
    ax_left.set_ylim((0, max(amps_fit) * 1.1))
    ax_left.grid()
    ax_left.set_xlabel('spread')
    ax_left.set_ylabel('min amp')
    ax_left.legend()
    fig_left.tight_layout()

    # x_low_left, x_high_left = 0.15, 1.1
    # p0_left = (0.05, 7, 0.013)
    # spreads_fit, amps_fit, amps_err_fit = trim(spreads, x_low_left, x_high_left, amps, amps_err)
    # popt_left, pcov_left = curve_fit(exp, spreads_fit, amps_fit, p0=p0_left, sigma=amps_err_fit, absolute_sigma=True)
    # print(popt_left)
    #
    # fig_left, ax_left = plt.subplots()
    # fit_xs = np.linspace(x_low_left, x_high_left, 1000)
    # ax_left.plot(fit_xs, exp(fit_xs, *p0_left), ls='--', color='gray', label='Initial Guess')
    # ax_left.plot(fit_xs, exp(fit_xs, *popt_left), color='red', label='Fit')
    # ax_left.errorbar(spreads_fit, amps_fit, yerr=amps_err_fit, marker='o', ls='none')
    # ax_left.set_ylim((0, max(amps_fit) * 1.1))
    # ax_left.grid()
    # ax_left.set_xlabel('spread')
    # ax_left.set_ylabel('min amp')
    # ax_left.legend()
    # fig_left.tight_layout()

    x_low_right, x_high_right = 1.1, 3.2
    p0_right = (0.01, -0.77, 0.013)
    spreads_fit, amps_fit, amps_err_fit = trim(spreads, x_low_right, x_high_right, amps, amps_err)
    popt_right, pcov_right = curve_fit(exp, spreads_fit, amps_fit, p0=p0_right, sigma=amps_err_fit, absolute_sigma=True)
    print(popt_right)

    fig_right, ax_right = plt.subplots()
    fit_xs = np.linspace(x_low_right, x_high_right, 1000)
    ax_right.plot(fit_xs, exp(fit_xs, *p0_right), ls='--', color='gray', label='Initial Guess')
    ax_right.plot(fit_xs, exp(fit_xs, *popt_right), color='red', label='Fit')
    ax_right.errorbar(spreads_fit, amps_fit, yerr=amps_err_fit, marker='o', ls='none')
    ax_right.set_ylim((0, max(amps_fit) * 1.1))
    ax_right.grid()
    ax_right.set_xlabel('spread')
    ax_right.set_ylabel('min amp')
    ax_right.legend()
    fig_right.tight_layout()

    x_low_full, x_high_full = -0.1, 3.1
    p0_full = (*popt_left, *popt_right[:-1])
    spreads_fit, amps_fit, amps_err_fit = trim(spreads, x_low_full, x_high_full, amps, amps_err)
    popt_full, pcov_full = curve_fit(exp_up_down, spreads_fit, amps_fit, p0=p0_full, sigma=amps_err_fit,
                                     absolute_sigma=True)
    print(popt_full)

    fig_full, ax_full = plt.subplots()
    fit_xs = np.linspace(x_low_full, x_high_full, 1000)
    ax_full.plot(fit_xs, exp_up_down(fit_xs, *p0_full), ls='--', color='gray', label='Initial Guess')
    ax_full.plot(fit_xs, exp_up_down(fit_xs, *popt_full), color='red', label='Fit')
    ax_full.errorbar(spreads_fit, amps_fit, yerr=amps_err_fit, marker='o', ls='none')
    ax_full.set_ylim((0, max(amps_fit) * 1.1))
    ax_full.grid()
    ax_full.set_xlabel('spread')
    ax_full.set_ylabel('min amp')
    ax_full.legend()
    fig_full.tight_layout()

    p0_full = (*popt_left, *popt_right[:-1])
    spreads_fit, amps_fit, amps_err_fit = trim(spreads, x_low_full, x_high_full, amps, amps_err)
    popt_full, pcov_full = curve_fit(exp_up_down, spreads_fit, amps_fit, p0=p0_full, sigma=amps_err_fit,
                                     absolute_sigma=True)
    print(popt_full)

    fig_full_res, ax_full_res = plt.subplots()
    ax_full_res.axhline(0, ls='--', color='black')
    ax_full_res.errorbar(spreads_fit, amps_fit - exp_up_down(spreads_fit, *popt_full), yerr=amps_err_fit, marker='o',
                         ls='none')
    ax_full_res.grid()
    ax_full_res.set_xlabel('spread')
    ax_full_res.set_ylabel('min amp')
    fig_full_res.tight_layout()

    plt.show()
    print('donzo')


def trim(xs, x_low, x_high, *ys):
    return (np.array(array) for array in zip(*[(xi, *yi) for xi, *yi in zip(xs, *ys) if x_low < xi < x_high]))


def exp_zero(x, a, b, c, d):
    return a * np.exp(-b * x) / (1 - np.exp(-c * x)) + d


def exp(x, a, b, c):
    return a * np.exp(-b * x) + c


def exp_up_down(x, a, b, c, d, e, f):
    return a * np.exp(-b * x) / (1 - np.exp(-c * x)) + d + e * np.exp(-f * x)


if __name__ == '__main__':
    main()
