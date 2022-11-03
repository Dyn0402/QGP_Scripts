#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 9:10 PM 2022
Created in PyCharm
Created as QGP_Scripts/ClustDist.py

@author: Dylan Neff, Dylan
"""

import numpy as np
from scipy.stats import rv_continuous
from scipy.stats import norm
from scipy.interpolate import interp1d


class ClustDist(rv_continuous):
    def __init__(self, means, sd, amp, a, b, base=1, wrap_num=1):
        self.means = means
        self.sd = sd
        self.base = base
        self.amp = amp  # Positive for clustering, negative for anti-clustering
        self.amp_norm = 1 / (np.sqrt(2 * np.pi) * self.sd)
        self.dists = []
        for mean in means:
            self.dists.append([norm(mean, sd)])
            for wrap_i in range(1, wrap_num + 1):
                self.dists[-1].extend([norm(mean - 2 * wrap_i * np.pi, sd), norm(mean + 2 * wrap_i * np.pi, sd)])
        # self.dists = [[norm(mean, sd), norm(mean - 2 * np.pi, sd), norm(mean + 2 * np.pi, sd)] for mean in means]
        self.n_points = 1000
        self.x = np.linspace(a, b, self.n_points)
        self.prob = np.ones(self.n_points)

        for dist in self.dists:
            clust_pdf = self.base
            for dist_sub in dist:
                clust_pdf += self.amp * dist_sub.pdf(self.x) / self.amp_norm
            self.prob *= clust_pdf
        self.prob /= np.sum(self.prob) / self.n_points * (b - a)

        # self.prob_interp = CubicSpline(self.x, self.prob)  # , bc_type='periodic')
        self.prob_interp = interp1d(self.x, self.prob)

        rv_continuous.__init__(self, a=a, b=b)

    def _pdf(self, x):
        return self.prob_interp(x)


class ClustDistPlusMinus(rv_continuous):
    def __init__(self, means, sd_plus, sd_minus, amp_plus, amp_minus, a, b, base=1, wrap_num=1):
        self.means = means
        self.base = base
        self.sd_plus = sd_plus
        self.sd_minus = sd_minus
        self.amp_plus = amp_plus  # Positive for clustering, negative for anti-clustering
        self.amp_minus = amp_minus  # Positive for clustering, negative for anti-clustering
        self.amp_norm_plus = 1 / (np.sqrt(2 * np.pi) * self.sd_plus)
        self.amp_norm_minus = 1 / (np.sqrt(2 * np.pi) * self.sd_minus)
        self.dists_plus, self.dists_minus = [], []

        for mean in means:
            self.dists_plus.append([norm(mean, sd_plus)])
            self.dists_minus.append([norm(mean, sd_minus)])
            for wrap_i in range(1, wrap_num + 1):
                self.dists_plus[-1].extend([norm(mean - 2 * wrap_i * np.pi, self.sd_plus),
                                            norm(mean + 2 * wrap_i * np.pi, self.sd_plus)])
                self.dists_minus[-1].extend([norm(mean - 2 * wrap_i * np.pi, self.sd_minus),
                                            norm(mean + 2 * wrap_i * np.pi, self.sd_minus)])
        # self.dists = [[norm(mean, sd), norm(mean - 2 * np.pi, sd), norm(mean + 2 * np.pi, sd)] for mean in means]
        self.n_points = 1000
        self.x = np.linspace(a, b, self.n_points)
        self.prob = np.ones(self.n_points)

        for dist_plus, dist_minus in zip(self.dists_plus, self.dists_minus):
            clust_pdf = self.base
            for dist_plus_sub in dist_plus:
                clust_pdf += self.amp_plus * dist_plus_sub.pdf(self.x) / self.amp_norm_plus
            for dist_minus_sub in dist_minus:
                clust_pdf += self.amp_minus * dist_minus_sub.pdf(self.x) / self.amp_norm_minus
            self.prob *= clust_pdf
        self.prob /= np.sum(self.prob) / self.n_points * (b - a)

        # self.prob_interp = CubicSpline(self.x, self.prob)  # , bc_type='periodic')
        self.prob_interp = interp1d(self.x, self.prob)

        rv_continuous.__init__(self, a=a, b=b)

    def _pdf(self, x):
        return self.prob_interp(x)


class ClustDist_slow(rv_continuous):
    def __init__(self, means, sd, amp, a, b, base=1, anti=True):
        self.means = means
        self.sd = sd
        self.base = base
        self.amp = amp
        self.amp_norm = 1 / (np.sqrt(2 * np.pi) * self.sd)
        self.dists = [[norm(mean, sd), norm(mean - 2 * np.pi, sd), norm(mean + 2 * np.pi, sd)] for mean in means]
        if anti:
            self.sign = -1
        else:
            self.sign = 1

        rv_continuous.__init__(self, a=a, b=b)

    def _pdf(self, x):
        try:
            prob_clust = np.ones(size=len(x))
        except TypeError:
            prob_clust = 1

        for dist in self.dists:
            clust_pdf = self.base
            for dist_sub in dist:
                clust_pdf += self.sign * self.amp * dist_sub.pdf(x) / self.amp_norm
            prob_clust *= clust_pdf

        return prob_clust


class Testrv(rv_continuous):
    def __init__(self, a, b):
        a = 0

        rv_continuous.__init__(self, a=a, b=b)

    def _pdf(self, x):
        norm = self.b**2 / 2

        return x / norm
