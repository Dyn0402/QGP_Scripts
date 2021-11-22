#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 08 2:08 PM 2021
Created in PyCharm
Created as QGP_Scripts/RefMultCorr

@author: Dylan Neff, dylan
"""

import numpy as np
import matplotlib.pyplot as plt


class RefMultCorr:
    def __init__(self, energy, ref=3):
        param_path = '/home/dylan/git/Research/QGP_Fluctuations/Tree_Reader/StRefMultCorr/Param.h'

