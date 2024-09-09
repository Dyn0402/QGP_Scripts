#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 09 2:00 PM 2024
Created in PyCharm
Created as QGP_Scripts/one_shot_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
num_vectors = 2  # Number of random momentum vectors
magnitude_range = (1, 10)  # Range for the magnitude of the momentum vectors

# Generate random 2D momentum vectors
angles = np.random.uniform(0, 2 * np.pi, num_vectors)
magnitudes = np.random.uniform(magnitude_range[0], magnitude_range[1], num_vectors)
momentum_vectors = np.vstack((magnitudes * np.cos(angles), magnitudes * np.sin(angles))).T

# Calculate total momentum
p_total = np.sum(momentum_vectors, axis=0)

# Function to rotate vector by angle theta
def rotate_vector(vector, theta):
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                                [np.sin(theta),  np.cos(theta)]])
    return np.dot(rotation_matrix, vector)

# Rotate each vector by -p_total/p_i
rotated_momentum_vectors = np.zeros_like(momentum_vectors)
for i, vector in enumerate(momentum_vectors):
    # Angle to rotate each vector
    theta = -np.arctan2(p_total[1], p_total[0]) / np.linalg.norm(vector)
    rotated_momentum_vectors[i] = rotate_vector(vector, theta)

# Recalculate total momentum after rotation
p_total_rotated = np.sum(rotated_momentum_vectors, axis=0)

# Print results
print(f"Initial Total Momentum: {p_total}")
print(f"Total Momentum after Rotation: {p_total_rotated}")

# Visualization
plt.figure(figsize=(8, 8))
plt.quiver([0]*num_vectors, [0]*num_vectors, momentum_vectors[:, 0], momentum_vectors[:, 1], color='r', angles='xy', scale_units='xy', scale=1, label='Original Vectors')
plt.quiver([0]*num_vectors, [0]*num_vectors, rotated_momentum_vectors[:, 0], rotated_momentum_vectors[:, 1], color='b', angles='xy', scale_units='xy', scale=1, label='Rotated Vectors')
plt.quiver(0, 0, p_total[0], p_total[1], color='g', scale_units='xy', scale=1, label='Initial Total Momentum', width=0.01)
plt.quiver(0, 0, p_total_rotated[0], p_total_rotated[1], color='m', scale_units='xy', scale=1, label='Rotated Total Momentum', width=0.01)
plt.xlim(-np.max(magnitudes), np.max(magnitudes))
plt.ylim(-np.max(magnitudes), np.max(magnitudes))
plt.grid()
plt.legend()
plt.title('2D Momentum Vectors and Total Momentum')
plt.xlabel('x-component of Momentum')
plt.ylabel('y-component of Momentum')
plt.show()
