#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 16 7:49 PM 2023
Created in PyCharm
Created as QGP_Scripts/c++_timing_plots.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt


def main():
    vec_sizes = [100, 250]
    for vec_size in vec_sizes:
        file_path = f'F:/Research/Results/c++_timing_tests/timing_results_vecsize{vec_size}.txt'

        # Load data from output text file
        with open(file_path, "r") as f:
            data = f.readlines()

        # Extract x and y values from data
        x = []
        y1 = []
        y2 = []
        for line in data:
            values = line.strip().split()
            x.append(int(values[0]))
            y1.append(float(values[1]) / 1000)
            y2.append(float(values[2]) / 1000)

        # Create plot
        fig, ax = plt.subplots()
        ax.plot(x, y1, label="Shuffle")
        ax.plot(x, y2, label="Erase")

        # Set plot title and axes labels
        ax.set_title(f"Timing Comparison {vec_size} Vector Size")
        ax.set_xlabel("Number of Samples")
        ax.set_ylabel("Elapsed Time (ms)")

        # Add legend
        ax.legend()

    # Show plot
    plt.show()

    print('donzo')


if __name__ == '__main__':
    main()
