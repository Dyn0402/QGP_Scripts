#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 05 5:04 PM 2023
Created in PyCharm
Created as QGP_Scripts/gaus_fits

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf

import uproot
import awkward as ak
import vector

from Measure import Measure


def main():
    energy = 7
    set_name = 'rapid05_resample_norotate_dca1_nsprx1_m2r6_m2s0_nhfit20_epbins1_0'
    path = f'F:/Research/Data/default_resample_epbins1/{set_name}/{energy}GeV/QA_{energy}GeV.root'

    cut_variables = init_params()

    with uproot.open(path) as file:
        for var_name, var_pars in cut_variables.items():
            hist = file[f'{var_pars["hist_name"]}{set_name}_{energy}']
            xs = hist.axis().centers()
            ys = hist.values()
            y_errs = hist.errors()
            sd = np.sqrt(np.average((xs - np.average(xs, weights=ys))**2, weights=ys))
            fig, ax = plt.subplots(dpi=144)
            ax.axhline(0, color='black')
            ax.set_title(f'{energy} GeV {var_name}')
            fig.canvas.manager.set_window_title(f'{energy} GeV {var_name}')
            plt.errorbar(xs, ys, y_errs, ls='none', marker='.')
            fit_range = var_pars['fit_range'][energy]
            xs_fit, ys_fit, y_errs_fit = zip(*[(x, y, y_err) for x, y, y_err in zip(xs, ys, y_errs)
                                               if fit_range[0] < x < fit_range[1]])
            try:
                popt, pcov = cf(gaus, xs_fit, ys_fit, sigma=y_errs_fit, absolute_sigma=True, p0=var_pars['p0'])
                popt_meas = [Measure(val, err) for val, err in zip(popt, np.sqrt(np.diag(pcov)))]
                print(popt_meas)
                ax.annotate(rf'$\mu$: {popt_meas[1]}' '\n' rf'$\sigma$: {popt_meas[2]}' f'\nsd: {sd:.4f}',
                            xy=(0.7, 0.95), xycoords='axes fraction',
                            bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
                            verticalalignment='top', horizontalalignment='center')
                xs_fit_plt = np.linspace(*fit_range, 100)
                plt.plot(xs_fit_plt, gaus(xs_fit_plt, *popt), color='red')
            except RuntimeError:
                ax.annotate(f'sd: {sd:.4f}',
                            xy=(0.7, 0.95), xycoords='axes fraction',
                            bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
                            verticalalignment='top', horizontalalignment='center')
    plt.show()
    print('donzo')


def init_params():
    cut_variables = {
        'dca': {
            'hist_name': 'pre_dca_',
            'fit_range': {
                7: [0.05, 0.3],
                11: [0.05, 0.3],
                19: [0.05, 0.3],
                27: [0.05, 0.3],
                39: [0.05, 0.3],
                62: [0.05, 0.3],
            },
            'p0': [1e5, 0.5, 0.2]
        },
        'nsigma_proton': {
            'hist_name': 'pre_nsigma_',
            'fit_range': {
                7: [-0.5, 0.5],
                11: [-0.5, 0.5],
                19: [-0.5, 0.5],
                27: [-0.5, 0.5],
                39: [-0.5, 0.5],
                62: [-0.5, 0.5],
            },
            'p0': [1e5, 0.0, 0.5]
        },
        'm2': {
            'hist_name': 'pre_m2_',
            'fit_range': {
                7: [0.8, 0.95],
                11: [0.8, 0.95],
                19: [0.8, 0.95],
                27: [0.8, 0.95],
                39: [0.8, 0.95],
                62: [0.8, 0.95],
            },
            'p0': [1e5, 1.0, 0.2]
        },
        'n_hits_fit': {
            'hist_name': 'pre_nhits_fit_',
            'fit_range': {
                7: [30, 45],
                11: [30, 44],
                19: [30, 44],
                27: [30, 44],
                39: [30, 44],
                62: [30, 44],
            },
            'p0': [1e6, 40, 5]
        },
        'pt': {
            'hist_name': 'pre_pt_',
            'fit_range': {
                7: [0.6, 1.0],
                11: [0.6, 1.0],
                19: [0.6, 1.0],
                27: [0.6, 1.0],
                39: [0.6, 1.0],
                62: [0.6, 1.0],
            },
            'p0': [1e5, 1.0, 0.5]
        },
        'p': {
            'hist_name': 'pre_p_',
            'fit_range': {
                7: [0.5, 1.6],
                11: [0.5, 1.6],
                19: [0.5, 1.6],
                27: [0.5, 1.6],
                39: [0.5, 1.6],
                62: [0.5, 1.6],
            },
            'p0': [1e5, 1.4, 0.5]
        },
    }

    return cut_variables


def gaus(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


if __name__ == '__main__':
    main()
