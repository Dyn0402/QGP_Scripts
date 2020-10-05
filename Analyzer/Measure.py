#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 2:58 PM 2020
Created in PyCharm
Created as QGP_Scripts/Measure.py

@author: Dylan Neff, dylan
"""


class Measure:
    def __init__(self, val=None, err=None):
        self.val = val
        self.err = err

    def __str__(self):
        return f'{self.val} ± {self.err}'
