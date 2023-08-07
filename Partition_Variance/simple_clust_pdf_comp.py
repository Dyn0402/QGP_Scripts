#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 06 6:01 PM 2023
Created in PyCharm
Created as QGP_Scripts/simple_clust_pdf_comp.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from integrate_pdf_var import base_gaus_pdf_wrap, get_norm_scipy


def main():
    mu = 1.2
    sigma = 0.1
    p = 0.5
    # amp = p / (1 - p) * (sigma * np.sqrt(2 * np.pi))
    amp = p / (1 - p) * np.sqrt(2 * np.pi) / sigma
    # amp = np.pi
    print(amp)
    widths_edges = np.linspace(0, 2 * np.pi, 1000)
    widths = (widths_edges[:-1] + widths_edges[1:]) / 2
    analytic_pdf = base_gaus_pdf_wrap(widths, mu, sigma, amp, 1, 8)
    analytic_pdf /= get_norm_scipy(base_gaus_pdf_wrap, [mu, sigma, amp, 1, 8])
    plt.plot(widths, analytic_pdf, label='Analytic')
    plt.plot(widths, gen_simple_clust_mc_pdf(mu, sigma, p, widths_edges, 10000000), label='Simulation')
    plt.ylim(bottom=0)
    plt.legend()
    plt.show()
    print('donzo')


def gen_simple_clust_mc_pdf(mu, sigma, p, width_edges, n=1000):
    gaus = np.random.normal(mu, sigma, int(p * n)) % (2 * np.pi)
    flat = np.random.uniform(0, 2 * np.pi, int((1 - p) * n))
    hist = np.histogram(gaus, width_edges)[0]
    hist += np.histogram(flat, width_edges)[0]
    hist = hist / (2 * np.pi * np.sum(hist) / len(hist))

    return hist


if __name__ == '__main__':
    main()
