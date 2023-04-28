#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 27 11:18 AM 2023
Created in PyCharm
Created as QGP_Scripts/momentum_conservation_model

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def main():
    # chat_gpt_test()
    momentum_test()
    print('donzo')


def momentum_test():
    n_particles = 50
    angle_frac = 0.8 / np.sqrt(n_particles)
    iterations = 4
    momenta = np.random.uniform(-5, 5, size=(n_particles, 3))

    net_momentum_vec = np.sum(momenta, axis=0)
    init_net_momentum = np.linalg.norm(net_momentum_vec)
    print(f'Initial: Net-Momentum={init_net_momentum}')
    plot_momenta(momenta, net_momentum_vec)
    net_momentum = init_net_momentum

    net_momentum_list = [net_momentum]
    for i in range(iterations):
        for j in range(len(momenta)):
            momenta[j] = rotate_vector(momenta[j], -net_momentum_vec, angle_frac * net_momentum / init_net_momentum)

        net_momentum_vec = np.sum(momenta, axis=0)
        net_momentum = np.linalg.norm(net_momentum_vec)
        print(f'Iteration {i}: Net-Momentum={net_momentum}')
        plot_momenta(momenta, net_momentum_vec)
        net_momentum_list.append(net_momentum)

    fig, ax = plt.subplots(dpi=144)
    ax.plot(range(len(net_momentum_list)), net_momentum_list)
    ax.axhline(0, color='black')
    ax.set_xlabel('Number of Iterations')
    ax.set_ylabel('Net Momentum Magnitude')
    fig.tight_layout()
    plt.show()


def chat_gpt_test():
    # Define the vectors
    vec1 = np.array([1, 1, 0.5])  # Initial unit vector to be rotated
    # vec1 /= np.linalg.norm(vec1)
    vec2 = np.array([0.1, 1, 0.2])  # Target unit vector
    # vec2 /= np.linalg.norm(vec2)

    rotated_vec1 = rotate_vector(vec1, vec2)

    # Create a 3D figure and add axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set limits for the axes
    ax.set_xlim([0, np.linalg.norm(vec2)])
    ax.set_ylim([0, np.linalg.norm(vec2)])
    ax.set_zlim([0, np.linalg.norm(vec2)])

    # Set the origin point
    ax.quiver(0, 0, 0, 0, 0, 0, color='black')

    # Plot the vectors
    ax.quiver(0, 0, 0, *(vec1 / np.linalg.norm(vec1)), color='red', length=np.linalg.norm(vec1), label='vec1')
    ax.quiver(0, 0, 0, *(vec2 / np.linalg.norm(vec2)), color='blue', length=np.linalg.norm(vec2), label='vec2')
    ax.quiver(0, 0, 0, *(rotated_vec1 / np.linalg.norm(rotated_vec1)), color='green',
              length=np.linalg.norm(rotated_vec1), label='rotated vec1')
    # ax.quiver(0, 0, 0, vec1[0], vec1[1], vec1[2], color='red', length=np.linalg.norm(vec1), label='vec1')
    # ax.quiver(0, 0, 0, vec2[0], vec2[1], vec2[2], color='blue', length=np.linalg.norm(vec2), label='vec2')
    # ax.quiver(0, 0, 0, rotated_vec1[0], rotated_vec1[1], rotated_vec1[2], color='green',
    #           length=np.linalg.norm(rotated_vec1), label='rotated vec1')
    print(vec1)
    print(vec2)
    print(rotated_vec1)

    # Add labels for the vectors
    ax.text(vec1[0], vec1[1], vec1[2], "vec1", color='red')
    ax.text(vec2[0], vec2[1], vec2[2], "vec2", color='blue')

    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.legend()

    # Show the plot
    plt.show()


def rotate_vector(vec, vec_target, angle_fraction=None):
    vec_unit = vec / np.linalg.norm(vec)
    vec_target_unit = vec_target / np.linalg.norm(vec_target)

    # Find the angle and axis of rotation
    angle = np.arccos(np.dot(vec_unit, vec_target_unit))
    angle = angle * angle_fraction if angle_fraction is not None else angle
    axis = np.cross(vec_unit, vec_target_unit)
    axis = axis / np.linalg.norm(axis)

    # Create the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    rot_matrix = np.array(
        [[axis[0] ** 2 * C + c, axis[0] * axis[1] * C - axis[2] * s, axis[0] * axis[2] * C + axis[1] * s],
         [axis[1] * axis[0] * C + axis[2] * s, axis[1] ** 2 * C + c, axis[1] * axis[2] * C - axis[0] * s],
         [axis[2] * axis[0] * C - axis[1] * s, axis[2] * axis[1] * C + axis[0] * s, axis[2] ** 2 * C + c]])

    # Apply the rotation matrix to vec1
    rotated_vec1 = np.dot(rot_matrix, vec)

    return rotated_vec1


def plot_vectors(vectors):
    # Create a 3D figure and add axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set the origin point
    ax.quiver(0, 0, 0, 0, 0, 0, color='black')

    # Plot the vectors
    max_axis = 0
    for vec in vectors:
        vec_len = np.linalg.norm(vec)
        ax.quiver(0, 0, 0, *(vec / vec_len), length=vec_len)
        for vec_i in vec:
            max_axis = vec_i if vec_i > max_axis else max_axis

    ax.set_xlim([-max_axis, max_axis])
    ax.set_ylim([-max_axis, max_axis])
    ax.set_zlim([-max_axis, max_axis])

    # Set labels for the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def plot_momenta(momenta, net_momentum_vec):
    plot_vectors(momenta)
    net_momentum = np.linalg.norm(net_momentum_vec)
    plt.quiver(0, 0, 0, *(net_momentum_vec / net_momentum), length=net_momentum, color='red')
    plt.quiver(0, 0, 0, *(-net_momentum_vec / net_momentum), length=net_momentum, color='orange')
    plt.tight_layout()


if __name__ == '__main__':
    main()
