#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 30 2:53 PM 2020
Created in PyCharm
Created as QGP_Scripts/exp_ws_func_test.py

@author: Dylan Neff, dylan
"""


import numpy as np
import matplotlib.pyplot as plt


def main():
    amp = 2
    width = 10
    center = 5
    decays = np.linspace(0, 1, 20)
    x = np.linspace(-100, 100, 1000)
    for decay in decays:
        y = exp_ws(x, amp, center, width, decay)
        plt.plot(x, y, label=f'decay = {decay}')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Exponential * Woods Saxon | amp={amp} width={width} center={center}')
    plt.grid()
    plt.show()

    print('donzo')


def exp_ws(x, amp, x_bar, width, decay):
    """
    Calculate exponential * Woods Saxon
    :param x: x value
    :param amp: amplitude
    :param x_bar: Center of both exponential and Woods Saxon functions
    :param width: Width of Woods Saxon function
    :param decay: Decay parameter for exponential
    :return: Result
    """
    return amp * np.exp(-decay / (2 * width) * (x - x_bar)) / (1 + np.exp(-(x - x_bar) / width))


if __name__ == '__main__':
    main()
