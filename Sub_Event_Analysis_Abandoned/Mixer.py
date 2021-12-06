#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 08 1:40 PM 2021
Created in PyCharm
Created as QGP_Scripts/Mixer

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt


class Mixer:
    """
    Collects events for specific dataset (energy) in classes and mixes to produce mixed set when enough available
    """
    def __init__(self):
        cent_bins = np.arange(0, 9)
        ep_bins = np.linspace(0, np.pi, 10)
        vz_bins = np.linspace(-50, 50, 10)
        classes = np.zeros((len(cent_bins), len(ep_bins), len(vz_bins)))

    def add_events(self, tracks):
        cents = tracks['cent']

    def add_event(self, ):
