#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 28 2:46 PM 2020
Created in PyCharm
Created as QGP_Scripts/rapidity_vs_pseudo.py

@author: Dylan Neff, dylan
"""


import numpy as np
import matplotlib.pyplot as plt


def main():
    pts = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]  # GeV
    etas = np.linspace(-1, 1, 1000)
    m = 0.93827  # GeV, proton mass
    for pt in pts:
        rap = rapidity(etas, pt, m)
        plt.plot(etas, rap, label=f'pt = {pt}')
    plt.plot(etas, etas, label='y = eta', color='black', linestyle='--')
    plt.legend()
    plt.xlabel('Pseudo-rapidity (eta)')
    plt.ylabel('Rapidity (y)')
    plt.title(f'Rapidity vs Pseudo-rapidity for m={m}GeV')
    plt.grid()
    plt.show()

    for pt in pts:
        rap = rapidity(etas, pt, m)
        ratio = [x / y for x, y in zip(rap, etas)]
        plt.plot(etas, rap/etas, label=f'pt = {pt}')
    plt.plot(etas, [1]*len(etas), label='y = eta', color='black', linestyle='--')
    plt.legend()
    plt.xlabel('Pseudo-rapidity (eta)')
    plt.ylabel('Rapidity / Pseudo-rapidity')
    plt.title(f'Rapidity / Pseudo-rapidity vs Pseudo-rapidity for m={m}GeV')
    plt.grid()
    plt.show()

    print('donzo')


def rapidity(eta, pt, m):
    """
    Calculate rapidity as a function of pesudo-rapidity, transverse momentum and mass
    :param eta: Pseudo-rapidity of particle trajectory
    :param pt: Transverse momentum of particle (units of energy matching m)
    :param m: Mass of particle (units of energy matching pt)
    :return: Rapidity of particle
    """
    return np.log((np.sqrt(m**2 + pt**2 * np.cosh(eta)**2) + pt * np.sinh(eta)) / np.sqrt(m**2 + pt**2))


if __name__ == '__main__':
    main()
