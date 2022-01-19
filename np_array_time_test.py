#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on January 17 4:18 PM 2022
Created in PyCharm
Created as QGP_Scripts/np_array_time_test.py

@author: Dylan Neff, Dylan
"""

import numpy as np
import matplotlib.pyplot as plt

from math import log10, floor
import timeit
import scipy.stats as sps

from DistStats import DistStats


def main():
    print('Test 1')
    dist_stats_rand_tests()
    print('Test 2')
    dist_stats_rand_tests2()
    return

    # compare_results3()
    # return

    nums = range(10)
    mults = range(len(nums))
    pows = list(range(6))
    sig_figs = 3
    num_times = 100
    print_res = True

    time_test1(nums, pows, num_times, sig_figs, print_res)
    time_test2(nums, mults, pows, num_times, sig_figs, print_res)

    nums_lists = [range(x) for x in range(1, 100)]
    mults_lists = [range(len(x)) for x in nums_lists]
    plot1(nums_lists, pows, num_times)
    plot2(nums_lists, mults_lists, pows, num_times)
    plot3(nums_lists, mults_lists, pows, num_times)
    plot4(nums_lists, mults_lists, pows, num_times)

    plt.legend()
    plt.show()

    print('donzo')


def plot1(nums_list, pows, num_times):
    nums_lens = []
    np_over_loops = []
    for nums in nums_list:
        loop_time, np_time = time_test1(nums, pows, num_times, print_res=False)
        nums_lens.append(len(nums))
        np_over_loops.append(np_time / loop_time * 100)

    plt.grid()
    plt.axhline(100, ls='--', color='black')
    plt.plot(nums_lens, np_over_loops, label='num ** pow')
    plt.xlabel('Number of elements in list to operate on')
    plt.ylabel('Percentage of time it takes numpy compared to loop')


def plot2(nums_list, mults_list, pows, num_times):
    nums_lens = []
    np_over_loops = []
    for nums, mults in zip(nums_list, mults_list):
        loop_time, np_time = time_test2(nums, mults, pows, num_times, print_res=False)
        nums_lens.append(len(nums))
        np_over_loops.append(np_time / loop_time * 100)

    plt.grid()
    plt.axhline(100, ls='--', color='black')
    plt.plot(nums_lens, np_over_loops, label='num ** pow * mult')
    plt.xlabel('Number of elements in list to operate on')
    plt.ylabel('Percentage of time it takes numpy compared to loop')


def plot3(nums_list, mults_list, pows, num_times):
    nums_lens = []
    np_over_loops = []
    for nums, mults in zip(nums_list, mults_list):
        loop_time, np_time = time_test3(nums, mults, pows, num_times, print_res=False)
        nums_lens.append(len(nums))
        np_over_loops.append(np_time / loop_time * 100)

    plt.grid()
    plt.axhline(100, ls='--', color='black')
    plt.plot(nums_lens, np_over_loops, label='sum(num ** pow * mult)')
    plt.xlabel('Number of elements in list to operate on')
    plt.ylabel('Percentage of time it takes numpy compared to loop')


def plot4(nums_list, mults_list, pows, num_times):
    nums_lens = []
    np_over_loops = []
    for nums, mults in zip(nums_list, mults_list):
        loop_time, np_time = time_test4(nums, mults, pows, num_times, print_res=False)
        nums_lens.append(len(nums))
        np_over_loops.append(np_time / loop_time * 100)

    plt.grid()
    plt.axhline(100, ls='--', color='black')
    plt.plot(nums_lens, np_over_loops, label='np 1Ds sum(num ** pow * mult)')
    plt.xlabel('Number of elements in list to operate on')
    plt.ylabel('Percentage of time it takes numpy compared to loop')


def time_test1(nums, pows, num_times, sig_figs=3, print_res=True):
    loop_time = timeit.timeit(lambda: loop_based(nums, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs loop based: {round_sig(loop_time, sig_figs)}s')

    np_time = timeit.timeit(lambda: np_based(nums, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs numpy based: {round_sig(np_time, sig_figs)}s')

        print(f'Numpy ran in {round_sig(np_time / loop_time * 100, sig_figs)}% of the time for {len(nums)} elements '
              f'and powers: {pows}')

    return loop_time, np_time


def time_test2(nums, mults, pows, num_times, sig_figs=3, print_res=True):
    loop_time = timeit.timeit(lambda: loop_based2(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs loop based: {round_sig(loop_time, sig_figs)}s')

    np_time = timeit.timeit(lambda: np_based2(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs numpy based: {round_sig(np_time, sig_figs)}s')

        print(f'Numpy ran in {round_sig(np_time / loop_time * 100, sig_figs)}% of the time for {len(nums)} elements '
              f'and powers: {pows}')

    return loop_time, np_time


def time_test3(nums, mults, pows, num_times, sig_figs=3, print_res=True):
    loop_time = timeit.timeit(lambda: loop_based3(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs loop based: {round_sig(loop_time, sig_figs)}s')

    np_time = timeit.timeit(lambda: np_based3(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs numpy based: {round_sig(np_time, sig_figs)}s')

        print(f'Numpy ran in {round_sig(np_time / loop_time * 100, sig_figs)}% of the time for {len(nums)} elements '
              f'and powers: {pows}')

    return loop_time, np_time


def time_test4(nums, mults, pows, num_times, sig_figs=3, print_res=True):
    loop_time = timeit.timeit(lambda: loop_based3(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs loop based: {round_sig(loop_time, sig_figs)}s')

    np_time = timeit.timeit(lambda: np_based4(nums, mults, pows), number=num_times)
    if print_res:
        print(f'Total time for {num_times} runs numpy based: {round_sig(np_time, sig_figs)}s')

        print(f'Numpy ran in {round_sig(np_time / loop_time * 100, sig_figs)}% of the time for {len(nums)} elements '
              f'and powers: {pows}')

    return loop_time, np_time


def compare_results():
    nums = [1, 2, 3]
    pows = [0, 1]

    loop_res = []
    for num in nums:
        loop_res.append([])
        for power in pows:
            loop_res[-1].append(num ** power)
    print(np.array(loop_res))

    nums = np.tile(nums, (len(pows), 1)).T
    print(nums ** pows)


def compare_results2():
    nums = [1, 2, 3]
    mults = [0, 1, 2]
    pows = [0, 1]

    loop_res = []
    for num, mult in zip(nums, mults):
        loop_res.append([])
        for power in pows:
            loop_res[-1].append(num ** power * mult)
    print(np.array(loop_res))

    nums = np.tile(nums, (len(pows), 1))
    print(((nums.T ** pows).T * mults).T)


def compare_results3():
    nums = [1, 2, 3]
    mults = [0, 1, 2]
    pows = [0, 1, 2, 3, 4]

    print(loop_based3(nums, mults, pows))
    print(np_based3(nums, mults, pows))


def round_sig(num, sig_figs=1):
    return round(num, -int(floor(log10(abs(num)))) + (sig_figs - 1))


def setup():
    pass


def loop_based(nums, pows):
    for num in nums:
        for power in pows:
            num ** power


def np_based(nums, pows):
    nums = np.tile(nums, (len(pows), 1)).T
    nums ** pows


def loop_based2(nums, mults, pows):
    for num, mult in zip(nums, mults):
        for power in pows:
            num ** power * mult


def np_based2(nums, mults, pows):
    nums = np.tile(nums, (len(pows), 1))
    (nums.T ** pows).T * mults


def loop_based3(nums, mults, pows):
    pow_list = {p: 0 for p in pows}
    for num, mult in zip(nums, mults):
        for power in pows:
            pow_list[power] += num ** power * mult

    return pow_list


def np_based3(nums, mults, pows):
    nums = np.tile(nums, (len(pows), 1))
    res = np.sum((nums.T ** pows).T * mults, axis=1)
    pow_list = {p: x for p, x in zip(pows, res)}

    return pow_list


# def loop_based3(nums, mults, pows):
#     pow_list = {p: 0 for p in pows}
#     for num, mult in zip(nums, mults):
#         for power in pows:
#             pow_list[power] += num ** power * mult
#
#     return pow_list


def np_based4(nums, mults, pows):
    pow_list = {p: np.sum(np.array(nums) ** p * mults) for p in pows}

    return pow_list


def dist_stats_test():
    test_dist = {x: y for x, y in zip(range(10), np.ones(10))}
    d = DistStats(test_dist)
    print(len(test_dist), d.get_raw_moment(5))

    print(list(test_dist.keys()))
    print(d.get_central_moment(4), sps.moment(list(test_dist.keys()), moment=4))

    test_dist.update({11: 0})
    d2 = DistStats(test_dist)
    print(len(test_dist), d2.get_raw_moment(5))

    dist = [8.77122009, 6.05084343, 5.1532753, 2.7667714, 3.81792581, 0.02699633, 3.08537567, 8.75912662, 0.99800769,
            2.4536896]
    ds = DistStats(dist, unbinned=True)
    print(ds.get_central_moment(2))
    ds2 = DistStats(dict(zip(dist, np.ones(len(dist), dtype=int))))
    print(ds2.get_central_moment(2))
    return


def dist_stats_rand_tests():
    num_tests = 1000
    num_vals = 100
    dist = sps.uniform(0, 10)
    central_moments = range(1, 10)
    max_deviation = 0.0000001

    for i in range(num_tests):
        vals = dist.rvs(num_vals)
        ds = DistStats(vals, unbinned=True)
        for cm in reversed(central_moments):
            ds_moment = ds.get_central_moment(cm)
            scipy_moment = sps.moment(vals, moment=cm)
            if ds_moment + scipy_moment == 0:
                if ds_moment == scipy_moment:
                    frac_deviation = 0
                else:
                    print('Div zero')
                    frac_deviation = 1
            else:
                frac_deviation = (ds_moment - scipy_moment) / ((ds_moment + scipy_moment) / 2)
            if frac_deviation > max_deviation:
                print(f'\nTest #{i}  Central Moment Order: {cm}\n'
                      f'ds = {ds_moment},  scipy = {scipy_moment},  frac_deviation = {frac_deviation}\n'
                      f'Data: {vals}')


def dist_stats_rand_tests2():
    num_tests = 100
    num_vals = 70
    dist = sps.uniform(0, 10)
    mult_dist = sps.uniform(5, 15)
    central_moments = range(1, 3)
    max_deviation = 0.0000001

    for i in range(num_tests):
        vals = dist.rvs(num_vals)
        mults = np.array(mult_dist.rvs(num_vals), dtype=int)
        ds = DistStats(dict(zip(vals, mults)))
        for cm in reversed(central_moments):
            ds_moment = ds.get_central_moment(cm)
            scipy_moment = sps.moment([val for val, mult in zip(vals, mults) for j in range(mult)], moment=cm)
            if ds_moment + scipy_moment == 0:
                if ds_moment == scipy_moment:
                    frac_deviation = 0
                else:
                    print('Div zero')
                    frac_deviation = 1
            else:
                frac_deviation = (ds_moment - scipy_moment) / ((ds_moment + scipy_moment) / 2)
            if frac_deviation > max_deviation:
                print(f'\nTest #{i}  Central Moment Order: {cm}\n'
                      f'ds = {ds_moment},  scipy = {scipy_moment},  frac_deviation = {frac_deviation}\n'
                      f'Data: {vals}\n Mults: {mults}')


if __name__ == '__main__':
    main()
