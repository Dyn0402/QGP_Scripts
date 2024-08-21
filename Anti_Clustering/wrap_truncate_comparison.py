#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 21 06:55 2024
Created in PyCharm
Created as QGP_Scripts/wrap_truncate_comparison

@author: Dylan Neff, dn277127
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    """
    Compare Gaussian distributions wrapped around the azimuth to those truncated.
    :return:
    """
    a, x0, sigma = 1, 0, 1.5
    n_points = 1000
    n_wrap = 5

    az_range_points = np.linspace(-np.pi, np.pi, n_points)
    bin_width = az_range_points[1] - az_range_points[0]

    truncated_gaus = gaus(az_range_points, a, x0, sigma)
    truncated_gaus_integral = np.sum(truncated_gaus) * bin_width
    truncated_gaus /= truncated_gaus_integral

    wrapped_gaus = gaus(az_range_points, a, x0, sigma)
    for wrap_i in range(n_wrap):
        neg_wrap_gaus = gaus(az_range_points - 2 * np.pi * (wrap_i + 1), a, x0, sigma)
        pos_wrap_gaus = gaus(az_range_points + 2 * np.pi * (wrap_i + 1), a, x0, sigma)
        wrapped_gaus += neg_wrap_gaus + pos_wrap_gaus
    wrapped_gaus_integral = np.sum(wrapped_gaus) * bin_width
    wrapped_gaus /= wrapped_gaus_integral

    fig, (ax, ax2) = plt.subplots(1, 2, sharey='all', figsize=(10, 4))
    ax.plot(az_range_points, truncated_gaus, label='Truncated Gaus')
    ax.plot(az_range_points, wrapped_gaus, ls='--', label='Wrapped Gaus')
    ax.set_title('Azimuth Centered')
    ax.legend()

    rotate_azimuth_angle = np.pi / 4
    truncated_gaus_rotated = shift_array(truncated_gaus, int(rotate_azimuth_angle / (2 * np.pi) * n_points))
    wrapped_gaus_rotated = shift_array(wrapped_gaus, int(rotate_azimuth_angle / (2 * np.pi) * n_points))
    # fig, ax = plt.subplots()
    ax2.plot(az_range_points, truncated_gaus_rotated, label='Truncated Gaus')
    ax2.plot(az_range_points, wrapped_gaus_rotated, ls='--', label='Wrapped Gaus')
    ax2.set_title('Azimuth Rotated')
    ax2.legend()
    fig.tight_layout()

    plt.show()

    print('donzo')


def shift_array(arr, n):
    return np.concatenate((arr[n:], arr[:n]))


def gaus(x, a, x0, sigma):
    return a * np.exp(-0.5 * (x - x0)**2 / sigma**2)


if __name__ == '__main__':
    main()
