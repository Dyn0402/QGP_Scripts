#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 2:55 PM 2020
Created in PyCharm
Created as QGP_Scripts/DistStats.py

@author: Dylan Neff, dylan
"""


from Measure import Measure
from scipy.special import binom


class DistStats:
    def __init__(self, dist=None):
        """
        Initiate with 1D distribution
        :param dist: Distribution to calculate stats for. dist ~ dictionary{key=x_value, value=counts}
        """
        self.dist = dist
        self.raw_moments = {}
        self.cent_moments = {}
        self.m = {}
        self.total_counts = None

    def calc_total_counts(self):
        self.calc_raw_moments(1, 0)

    def calc_raw_moments(self, n_min=1, n_max=1):
        """
        Calculate raw moments of self.dist from n_min to n_max order
        :param n_min: Lowest order raw moment to calculate
        :param n_max: Highest order raw moment to calculate
        :return:
        """
        if len([n for n in range(n_min, n_max+1) if n not in self.raw_moments]) == 0:
            return

        if n_min <= 0:
            n_min = 1
        self.raw_moments.update({0: 1})

        for ni in range(n_min, n_max+1):
            self.raw_moments.update({ni: 0})

        self.total_counts = 0
        for x, counts in self.dist.items():
            self.total_counts += counts
            for ni in range(n_min, n_max+1):
                self.raw_moments[ni] += x**ni * counts

        for ni in range(n_min, n_max+1):
            self.raw_moments[ni] /= self.total_counts

    def calc_cent_moments(self, n_min=1, n_max=1):
        """
        Calculate central moments of self.dist from n_min to n_max order from self.raw_moments
        :param n_min: Lowest order central moment to calculate
        :param n_max: Highest order central moment to calculate
        :return:
        """
        for ni in range(n_min, n_max+1):
            if ni not in self.raw_moments:
                self.calc_raw_moments(n_min, n_max)
                break

        for ni in range(n_min, n_max+1):
            self.cent_moments.update({ni: 0})
            for nj in range(0, ni+1):
                self.cent_moments[ni] += \
                    binom(ni, nj) * (-1)**(ni - nj) * self.raw_moments[nj] * self.raw_moments[1]**(ni - nj)

        for ni in range(2, n_max+1):
            self.m.update({ni: self.cent_moments[ni] / self.cent_moments[2]**(0.5 * ni)})

    def get_mean(self):
        """
        Get mean of distribution with error from delta theorem
        :return: Mean as a Measure
        """
        self.calc_cent_moments(1, 2)
        val = self.raw_moments[1]
        err = (self.cent_moments[2] / self.total_counts)**0.5

        return Measure(val, err)

    def get_sd(self):
        """
        Get standard deviation of distribution with error from delta theorem
        :return: Standard deviation as a Measure
        """
        self.calc_cent_moments(1, 4)
        val = self.cent_moments[2]**0.5
        err = ((self.m[4] - 1) * self.cent_moments[2] / (4 * self.total_counts))**0.5

        return Measure(val, err)

    def get_skewness(self):
        """
        Get skewness of distribution with error from delta theorem
        :return: Skewness as a Measure
        """
        self.calc_cent_moments(1, 6)
        m = self.m
        val = self.cent_moments[3] / self.cent_moments[2]**1.5
        err = ((9 - 6 * m[4] + m[3]**2 * (35 + 9 * m[4]) / 4 - 3 * m[3] * m[5] + m[6])
               / self.total_counts) ** 0.5

        return Measure(val, err)

    def get_kurtosis(self):
        """
        Get kurtosis of distribution with error from delta theorem
        :return: Kurtosis as a Measure
        """
        self.calc_cent_moments(1, 8)
        m = self.m
        val = self.cent_moments[4] / self.cent_moments[2]**2 - 3
        # for i, rawi in self.raw_moments.items():
        #     print(f'kurt raw{i}: {rawi}')
        # for i, centi in self.cent_moments.items():
        #     print(f'kurt cent{i}: {centi}')
        # for i, mi in m.items():
        #     print(f'kurt m{i}: {mi}')
        err = ((-m[4]**2 + 4 * m[4]**3 + 16 * m[3]**2 * (1 + m[4]) - 8 * m[3] * m[5] - 4 * m[4] * m[6] + m[8])
               / self.total_counts) ** 0.5

        return Measure(val, err)

    def get_kurt_var(self):
        """
        Get kurtosis*variance of distribution with error from delta theorem
        :return: Kurtosis*variance as a Measure
        """
        self.calc_cent_moments(1, 8)
        m = self.m
        val = self.cent_moments[4] / self.cent_moments[2] - 3 * self.cent_moments[2]
        err = ((-9 + 6 * m[4]**2 + m[4]**3 + 8 * m[3]**2 * (5 + m[4]) - 8 * m[3] * m[5] + m[4] * (9 - 2 * m[6])
                - 6 * m[6] + m[8]) * self.cent_moments[2]**2 / self.total_counts) ** 0.5

        return Measure(val, err)
