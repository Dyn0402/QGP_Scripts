#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 06 12:52 AM 2020
Created in PyCharm
Created as QGP_Scripts/ref3_corr_62GeV_vz.py

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    print(vz_correction(0.845083, '11GeV_set1'))
    print(vz_correction(0.0880397, '11GeV_set1'))
    plot_corr_11()
    print('donzo')


def plot_corr_62():
    vz = np.linspace(-30, 30, 5000)
    corr1 = [vz_correction(z, 'set1') for z in vz]
    corr2 = [vz_correction(z, 'set2') for z in vz]
    corr2_corrected = [vz_correction(z, 'set2_corrected') for z in vz]
    plt.plot(vz, corr1, color='blue', label='set1')
    plt.plot(vz, corr2, color='red', label='set2')
    plt.semilogy()
    plt.title('Refmult3 Vz Correction Factor (hovno) vs vz')
    plt.xlabel('z vertex position (cm)')
    plt.ylabel('refmult3 correction factor "hovno"')
    plt.legend()
    plt.show()
    plt.plot(vz, corr1, color='blue', label='set1')
    plt.plot(vz, corr2, color='red', label='set2')
    plt.plot(vz, corr2_corrected, color='green', label='set2_corrected')
    plt.title('Refmult3 Vz Correction Factor (hovno) vs vz')
    plt.xlabel('z vertex position (cm)')
    plt.ylabel('refmult3 correction factor "hovno"')
    plt.legend()
    plt.show()


def plot_corr_11():
    vz = np.linspace(-30, 30, 5000)
    corr1 = [vz_correction(z, '11GeV_set1') for z in vz]
    corr2 = [vz_correction(z, '11GeV_set2') for z in vz]
    plt.plot(vz, corr1, color='blue', label='set1')
    plt.plot(vz, corr2, color='red', label='set2')
    plt.title('Refmult3 Vz Correction Factor (hovno) vs vz')
    plt.xlabel('z vertex position (cm)')
    plt.ylabel('refmult3 correction factor "hovno"')
    plt.legend()
    plt.show()


def vz_correction(z, s='set1'):
    pars = {'set1': [663.421, 0.661982, 0.0220411, -0.00173208, -0.000123877, 1.71494E-06, 1.04242E-07, 8.911],
            'set2': [672.332, 0.354958, -0.0108052, 0.000160221, 3.39671E-06, -1.46475, -2.06264E-09, 0],
            'set2_corrected': [672.332, 0.354958, -0.0108052, 0.000160221, 3.39671E-06, -1.46475e-07, -2.06264E-09, 0],
            '11GeV_set2': [411.646,0.036548,-0.00333334,1.68893e-05,-1.0805e-06,-1.95689e-09,1.89723e-10,0],
            '11GeV_set1': [417.985,0.0264489,0.000980243,0.000123403,-4.58899e-06,-4.17197e-08,9.44972e-10,-6.339]}

    hovno = 1.0
    refmult_ref = pars[s][0]
    refmult_z = 0
    for i in range(7):
        refmult_z += pars[s][i] * z**i

    if refmult_z > 0:
        hovno = (refmult_ref + pars[s][7]) / refmult_z

    return hovno


if __name__ == '__main__':
    main()
