#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 14 10:55 AM 2022
Created in PyCharm
Created as QGP_Scripts/dist_stat_tests.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt
from DistStats import DistStats


def main():
    dist = {0: 553, 1: 141, 2: 1030, 3: 3132, 4: 5493, 5: 8353, 6: 5644, 7: 2670, 8: 1134, 9: 276, 10: 142, 11: 32,
            12: 18, 13: 26, 14: 110, 15: 46}
    stats = DistStats(dist)
    print(stats.get_cumulant(4))
    print('donzo')


if __name__ == '__main__':
    main()
