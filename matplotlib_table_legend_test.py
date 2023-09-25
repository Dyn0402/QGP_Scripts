#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 24 7:46 PM 2023
Created in PyCharm
Created as QGP_Scripts/matplotlib_table_legend_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    two_leg_test()
    print('donzo')


def two_leg_test():
    # Create some sample data and plots
    x = [1, 2, 3, 4, 5]
    y1 = [1, 2, 3, 4, 5]
    y2 = [5, 4, 3, 2, 1]

    fig, ax = plt.subplots()

    # Plot the data and set labels
    line1, = ax.plot(x, y1, label='Data_Set_1')
    line2, = ax.plot(x, y2, label='Data_Set_2')

    # Get the handles and labels of all lines
    handles, labels = ax.get_legend_handles_labels()

    # Create the first legend with half of the entries
    legend1 = ax.legend(handles[:len(handles) // 2], labels[:len(labels) // 2], loc='upper left')

    # Add the second legend with the other half of the entries
    ax.add_artist(legend1)
    legend2 = ax.legend(handles[len(handles) // 2:], labels[len(labels) // 2:], loc='upper right')

    plt.show()

def table_test():
    # Define your categories and energies
    categories = ["STAR", "AMPT"]
    energies = [7, 11, 19, 27, 39, 62]

    # Create an empty legend dictionary to store handles and labels
    legend_dict = {}

    # Define a list of marker colors to use
    marker_colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']

    # Your existing loop for plotting error bars
    for i, data_set in enumerate(categories):
        for j, energy in enumerate(energies):
            # Simulate error bar data (replace with your actual data)
            x = np.arange(10)
            y = np.random.rand(10)
            yerr = np.random.rand(10)

            # Plot the error bar with markers and store the handle
            handle = plt.errorbar(x, y, yerr=yerr, fmt='o', label=f"{data_set} - {energy}", color=marker_colors[j],
                                  markersize=10)

            # Store the marker color and handle in the legend dictionary
            legend_dict[(i, j)] = (marker_colors[j], handle)

    # Create a list to hold the legend handles and labels
    legend_handles = []
    legend_labels = []

    # Add the first row with energy labels
    for energy in energies:
        legend_labels.append(f"{energy}")
        legend_handles.append(plt.Line2D([0], [0], marker='', color='w', label=f"{energy}"))

    # Add the second row with "STAR" label in the first column
    legend_labels.append("STAR")
    legend_handles.append(plt.Line2D([0], [0], marker='', color='w', label="STAR"))
    for j, energy in enumerate(energies):
        legend_labels.append("")
        marker_color, handle = legend_dict[(0, j)]
        legend_handles.append(handle)

    # Add the third row with "AMPT" label in the first column
    legend_labels.append("AMPT")
    legend_handles.append(plt.Line2D([0], [0], marker='', color='w', label="AMPT"))
    for j, energy in enumerate(energies):
        legend_labels.append("")
        marker_color, handle = legend_dict[(1, j)]
        legend_handles.append(handle)

    # Create the legend with a 3x6 grid layout
    fig, ax = plt.subplots()
    ax.legend(legend_handles, legend_labels, loc='upper center', ncol=6, bbox_to_anchor=(0.5, 1.1))

    # Hide the axis to display only the legend
    ax.axis('off')

    plt.show()


if __name__ == '__main__':
    main()
