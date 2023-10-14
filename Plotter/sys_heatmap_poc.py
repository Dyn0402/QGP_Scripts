#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 13 4:58 PM 2023
Created in PyCharm
Created as QGP_Scripts/sys_heatmap_poc

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

# Define your systematic variables and values for each energy and centrality
systematic_variables = ['default', 'vz', 'm2r', 'nhfit', 'Efficiency', 'dca', 'nsprx', 'dcxyqa', 'pileupqa',
                        'sysrefshift']
energies = [7.7, 11.5, 19.6, 27, 39, 62.4]
# Determine the color limits based on the minimum and maximum values across all heatmaps
vmin = 0
vmax = 100

cent_sets = [[8, 7, 6, 5, 4], [3, 2, 1, 0]]

for centralities in cent_sets:
    errors = np.random.rand(len(systematic_variables), len(centralities),
                            len(energies))  # Replace with your actual error values

    # Create subplots for each centrality
    fig, axes = plt.subplots(1, len(centralities), figsize=(14, 5), dpi=144, sharey=True)

    # Iterate through centralities and create a heatmap for each
    for i, centrality in enumerate(centralities):
        ax = axes[i]
        ax.set_title(f'Centrality {centrality}')
        hm = sns.heatmap(errors[:, i, :] * 100, cmap='bwr', annot=True, fmt='.0f', xticklabels=energies,
                         yticklabels=systematic_variables, ax=ax, vmin=vmin, vmax=vmax, cbar=False)
        for t in ax.texts:
            t.set_text(t.get_text() + "%")
        ax.set_title(f'Centrality {centrality}')  # Keep subfigure titles
        # Remove tick marks on the outside of the axes
        ax.tick_params(left=False, right=False, top=False, bottom=False)


        def format_percent(x, pos):
            return f'{x:.1f}%'


        ax.yaxis.set_major_formatter(mtick.FuncFormatter(format_percent))
        ax.axhline(y=1, color='white', linewidth=10)
        # ax.axhline(y=1, color='black', linewidth=2)

    axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
    plt.tight_layout()
plt.show()
