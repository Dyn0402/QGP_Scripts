#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 07 9:40 PM 2022
Created in PyCharm
Created as QGP_Scripts/interp_test

@author: Dylan Neff, Dyn04
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp


def main():
    scipy_doc_test_edit()

    print('donzo')


def gaus2d_test():
    x = np.linspace(0, 10, 50)
    y = np.linspace(-6, -1, 10)
    xv, yv = np.meshgrid(x, y)
    xx = np.ravel(xv)
    yy = np.ravel(yv)

    mu = np.array([2, -2])
    sigma = np.array([[0.5, 0], [0, 1]])

    # print(np.column_stack((xx, yy)))

    z = np.array([gaus(xyi, mu, sigma) for xyi in np.column_stack((xx, yy))])
    # print(z)
    print(z.shape)
    print(z.reshape((len(x), len(y))))

    plt.figure()
    plt.title('Original')
    plt.imshow(z.reshape((len(y), len(x))), extent=[min(x), max(x), min(y), max(y)], aspect='auto', origin='lower',
               cmap='jet')

    print(gaus(np.array([2, -3]), mu, sigma))

    tcks = interp.bisplrep(xx, yy, z, kx=1, ky=1)

    # f = interp2d(df_chi['amp'], df_chi['spread'], df_chi['chi2_sum'], kind='cubic')
    x1 = np.linspace(min(x), max(x), 10)
    y1 = np.linspace(min(y), max(y), 5)
    xx1, yy1 = np.meshgrid(x1, y1)
    # z = f(x, y)
    z1 = interp.bisplev(x1, y1, tcks)
    print(z1)
    print(z1.T)
    print(z1.shape)
    print(z1.T.shape)

    plt.figure()
    plt.title('Interpolated')
    plt.imshow(z1.T, extent=[min(x1), max(x1), min(y1), max(y1)], aspect='auto', origin='lower', cmap='jet')

    plt.show()


def gaus(x, mu, sigma):
    return np.exp(-0.5 * np.dot(np.dot((x - mu), np.linalg.inv(sigma)), (x - mu).T))


def scipy_doc_test_edit():
    x = np.arange(0, 4, 0.25)
    y = np.array([0, 2, 3, 3.5, 3.75, 3.875, 3.9375, 4])
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.pi * X / 2) * np.exp(Y / 2)

    x2 = np.arange(0, 4, 0.05)
    y2 = np.arange(0, 4, 0.05)
    f = interp.interp2d(x, y, Z, kind='cubic')
    Z2 = f(x2, y2)

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].pcolormesh(X, Y, Z, cmap='jet')

    X2, Y2 = np.meshgrid(x2, y2)
    ax[1].pcolormesh(X2, Y2, Z2, cmap='jet')

    fig2, ax2 = plt.subplots()
    x2_index = dict((val, key) for key, val in dict(enumerate(x2)).items())
    ax2.grid()
    for x_index, x_val in enumerate(x):
        ax2.scatter(y, Z.T[x_index], label=f'x: {x_val}')
    ax2.legend()
    ax2.set_xlabel('y')
    print(x)
    print(x2)

    plt.show()


def scipy_doc_test():
    x = np.linspace(0, 4, 13)
    y = np.array([0, 2, 3, 3.5, 3.75, 3.875, 3.9375, 4])
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.pi * X / 2) * np.exp(Y / 2)

    x2 = np.linspace(0, 4, 65)
    y2 = np.linspace(0, 4, 65)
    f = interp.interp2d(x, y, Z, kind='cubic')
    Z2 = f(x2, y2)

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].pcolormesh(X, Y, Z)

    X2, Y2 = np.meshgrid(x2, y2)
    ax[1].pcolormesh(X2, Y2, Z2)

    plt.show()

if __name__ == '__main__':
    main()
