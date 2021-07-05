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
import warnings


class DistStats:
    def __init__(self, dist=None, debug=True):
        """
        Initiate with 1D distribution
        :param dist: Distribution to calculate stats for. dist ~ dictionary{key=x_value, value=counts}
        """
        self.dist = dist
        self.raw_moments = {}
        self.cent_moments = {}
        self.m = {}
        self.total_counts = None
        if debug:
            warnings.filterwarnings('error')

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
            try:
                self.raw_moments[ni] /= self.total_counts
            except RuntimeWarning:
                print(f'Warning treated as error for debug:\n'
                      f'raw_moment {ni}: {self.raw_moments[ni]}, total_counts: {self.total_counts}')

    def calc_cent_moments(self, n_min=1, n_max=1):
        """
        Calculate central moments of self.dist from n_min to n_max order from self.raw_moments
        :param n_min: Lowest order central moment to calculate
        :param n_max: Highest order central moment to calculate
        :return:
        """
        for ni in range(n_min, n_max+1):
            if ni not in self.raw_moments:
                self.calc_raw_moments(n_min, n_max+1)
                break

        for ni in range(n_min, n_max+1):
            self.cent_moments.update({ni: 0})
            for nj in range(0, ni+1):
                self.cent_moments[ni] += \
                    binom(ni, nj) * (-1)**(ni - nj) * self.raw_moments[nj] * self.raw_moments[1]**(ni - nj)

        for ni in range(2, n_max+1):
            if self.cent_moments[2] == 0:
                self.m.update({ni: float('nan')})
            else:
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
        err = (self.m[4] - 1) * self.cent_moments[2]
        if err >= 0:
            err = (err / (4 * self.total_counts))**0.5
        else:
            err = float('nan')

        return Measure(val, err)

    def get_skewness(self):
        """
        Get skewness of distribution with error from delta theorem
        :return: Skewness as a Measure
        """
        self.calc_cent_moments(1, 6)
        m = self.m
        if self.cent_moments[2] == 0:
            val = float('nan')
            err = float('nan')
        else:
            val = self.cent_moments[3] / self.cent_moments[2]**1.5
            err = 9 - 6 * m[4] + m[3]**2 * (35 + 9 * m[4]) / 4 - 3 * m[3] * m[5] + m[6]
            if err >= 0:
                err = (err / self.total_counts) ** 0.5
            else:
                err = float('nan')

        return Measure(val, err)

    def get_kurtosis(self):
        """
        Get kurtosis of distribution with error from delta theorem
        :return: Kurtosis as a Measure
        """
        self.calc_cent_moments(1, 8)
        m = self.m
        if self.cent_moments[2] == 0:
            val = float('nan')
            err = float('nan')
        else:
            val = self.cent_moments[4] / self.cent_moments[2]**2 - 3
            err = -m[4]**2 + 4 * m[4]**3 + 16 * m[3]**2 * (1 + m[4]) - 8 * m[3] * m[5] - 4 * m[4] * m[6] + m[8]
            if err >= 0:
                err = (err / self.total_counts) ** 0.5
            else:
                err = float('nan')

        return Measure(val, err)

    def get_kurt_var(self):
        """
        Get kurtosis*variance of distribution with error from delta theorem
        :return: Kurtosis*variance as a Measure
        """
        self.calc_cent_moments(1, 8)
        m = self.m
        if self.cent_moments[2] == 0:
            val = float('nan')
            err = float('nan')
        else:
            val = self.cent_moments[4] / self.cent_moments[2] - 3 * self.cent_moments[2]
            err = -9 + 6 * m[4]**2 + m[4]**3 + 8 * m[3]**2 * (5 + m[4]) - 8 * m[3] * m[5] + m[4] * (9 - 2 * m[6]) \
                  - 6 * m[6] + m[8]
            if err >= 0:
                err = (err * self.cent_moments[2]**2 / self.total_counts) ** 0.5
            else:
                err = float('nan')

        return Measure(val, err)

    def get_cumulant(self, order):
        """
        Calculate cumulant of distribution for given order from central moments with error from delta theorem
        !!!! DELTA THEOREM ERRORS NOT YET IMPLEMENTED !!!!!
        :param order: Order of cumulant to calculate
        :return: Cumulant value and error as Measure object
        """
        self.calc_cent_moments(1, order)
        # Error catching for delta theorem error calc
        err = 1
        # Maybe write in bell polynomials later
        if order == 1:
            val = self.raw_moments[1]
        elif order == 2:
            val = self.cent_moments[2]
        elif order == 3:
            val = self.cent_moments[3]
        elif order == 4:
            val = self.cent_moments[4] - 3 * self.cent_moments[2]**2
        elif order == 5:
            val = self.cent_moments[5] - 10 * self.cent_moments[2] * self.cent_moments[3]

        return Measure(val, err)

    def get_k_stat(self, order):
        """
        Calculate k statistic of distribution for given order from central moments
        !!!! ERRORS NOT YET IMPLEMENTED !!!!!
        :param order: Order of k statistic to colculate
        :return: K statistic value and error as Measure object
        """
        self.calc_raw_moments(1, order)
        # Placeholder for actual error calculation
        err = 1
        # Maybe try to find analytical formula later
        n = self.total_counts
        if n < order:  # Will get /0 errors
            print(f'Too few entries ({n}) for {order} order k-statistic!')
            return Measure(float('NaN'), float('NaN'))
        if order == 1:
            val = self.raw_moments[1]
        elif order == 2:
            val = n / (n - 1) * self.cent_moments[2]
        elif order == 3:
            val = n**2 / ((n - 1) * (n - 2)) * self.cent_moments[3]
        elif order == 4:
            val = n**2 * ((n + 1) * self.cent_moments[4] - 3 * (n - 1) * self.cent_moments[2] ** 2) \
                  / ((n - 1) * (n - 2) * (n - 3))

        return Measure(val, err)
