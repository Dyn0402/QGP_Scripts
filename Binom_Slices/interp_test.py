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
from matplotlib.pyplot import cm
import scipy.interpolate as interp
import pandas as pd


def main():
    # scipy_doc_test_edit()
    chi2_data_test()

    print('donzo')


def chi2_data_test():
    base_path = 'D:/Research/Results/Azimuth_Analysis/'
    chi_df_name = 'chi_df.csv'
    df = pd.read_csv(base_path + chi_df_name)
    print(df)
    print(df.columns)
    for spread in pd.unique(df['spread']):
        print(f'spread {spread}, amps: {list(df[df["spread"] == spread]["amp"])}')
    return
    df = df[(df['spread'] != 0) & (df['spread'] != 4)]  # & (df['spread'] != 0.2) & (df['spread'] != 0.5) &
            # (df['spread'] != 1)]
    df = df.sort_values(by=['spread', 'amp'])

    print(df['amp'], df['spread'], np.log10(df['chi2_sum']))
    x, y, z = [np.asarray(d) for d in [df['amp'], df['spread'], np.log10(df['chi2_sum'])]]
    print(x, y, z)
    # f = interp.interp2d(x, y, z, kind='cubic')

    fig1 = plt.figure()
    x_unq, y_unq = np.unique(x), np.unique(y)
    X, Y = np.meshgrid(np.unique(x), np.unique(y))
    Z = []
    for spread in np.unique(df['spread']):
        Z.append(list(np.log10(df[df['spread'] == spread]['chi2_sum'])))
    pcm1 = plt.pcolormesh(X, Y, Z, cmap='jet')
    plt.xlabel('amp')
    plt.ylabel('spread')
    fig1.colorbar(pcm1)
    fig1.tight_layout()

    print(np.asarray(Z).shape, np.unique(x).shape, np.unique(y).shape)

    rbs = interp.RectBivariateSpline(sorted(np.unique(x)), sorted(np.unique(y)), np.asarray(Z).T, kx=2, ky=2)
    f = lambda x, y: rbs.ev(x, y)

    print(f(0.1, 1))

    fig2 = plt.figure()
    x_intp = np.linspace(0, .2, 100)
    y_intp = np.linspace(0.2, 2.5, 100)
    X_intp, Y_intp = np.meshgrid(x_intp, y_intp)
    Z_intp = []
    for yi in y_intp:
        Zi_intp = []
        for xi in x_intp:
            Zi_intp.append(f(xi, yi))
        Z_intp.append(Zi_intp)
    pcm2 = plt.pcolormesh(X_intp, Y_intp, Z_intp, cmap='jet')
    plt.xlabel('amp')
    plt.ylabel('spread')
    fig2.colorbar(pcm2)
    fig2.tight_layout()

    fig3 = plt.figure()
    spreads = np.unique(df['spread'])
    color = iter(cm.rainbow(np.linspace(0, 1, len(spreads))))
    for spread in spreads:
        c = next(color)
        df_s = df[df['spread'] == spread]
        print(f'spread {spread}, num: {len(df_s)}')
        plt.scatter(df_s['amp'], np.log10(df_s['chi2_sum']), color=c, label=f'spread {spread}')
        zs = []
        amps = np.linspace(min(df_s['amp']), max(df_s['amp']), 1000)
        for amp in amps:
            zs.append(f(amp, spread))
        plt.plot(amps, zs, color=c)
    plt.xlabel('amp')
    plt.ylabel('chi2_sum')
    plt.legend()
    fig3.tight_layout()

    fig4 = plt.figure()
    amps = np.unique(df['amp'])
    color = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
    for amp in amps:
        c = next(color)
        df_a = df[df['amp'] == amp]
        print(f'amp {amp}, num: {len(df_a)}')
        plt.scatter(df_a['spread'], np.log10(df_a['chi2_sum']), color=c, label=f'amp {amp}')
        zs = []
        spreads = np.linspace(min(df_a['spread']), max(df_a['spread']), 1000)
        for spread in spreads:
            zs.append(f(amp, spread))
        plt.plot(spreads, zs, color=c)
    plt.xlabel('spread')
    plt.ylabel('chi2_sum')
    plt.legend(bbox_to_anchor=(1.05, 1))
    fig4.tight_layout()

    plt.show()


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
    print('x: ', x)
    print('y: ', y)
    print('Z: ', Z)
    f = interp.interp2d(x, y, Z, kind='linear')
    Z2 = f(x2, y2)

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].pcolormesh(X, Y, Z, cmap='jet')

    X2, Y2 = np.meshgrid(x2, y2)
    ax[1].pcolormesh(X2, Y2, Z2, cmap='jet')

    fig2, ax2 = plt.subplots()
    x2_index = dict((val, key) for key, val in dict(enumerate(x2)).items())
    ax2.grid()
    color = iter(cm.rainbow(np.linspace(0, 1, len(x))))
    for x_index, x_val in enumerate(x):
        c = next(color)
        ax2.scatter(y, Z.T[x_index], label=f'x: {x_val}', color=c)
        z_xi = []
        for yi in y:
            # print(xi, y_val, f(xi, y_val)[0])
            z_xi.append(f(x_val, yi)[0])
        # print(f(x, np.full(len(x), y_val)))
        ax2.plot(y, z_xi, color=c)
    ax2.legend()
    ax2.set_xlabel('y')

    fig3, ax3 = plt.subplots()
    x2_index = dict((val, key) for key, val in dict(enumerate(x2)).items())
    ax3.grid()
    color = iter(cm.rainbow(np.linspace(0, 1, len(y))))
    for y_index, y_val in enumerate(y):
        c = next(color)
        ax3.scatter(x, Z[y_index], label=f'y: {y_val}', color=c)
        print(f'y={y_val}')
        z_yi = []
        for xi in x:
            print(xi, y_val, f(xi, y_val)[0])
            z_yi.append(f(xi, y_val)[0])
        # print(f(x, np.full(len(x), y_val)))
        ax3.plot(x, z_yi, color=c)
    ax3.legend()
    ax3.set_xlabel('y')
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
