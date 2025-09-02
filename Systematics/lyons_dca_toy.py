#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on August 27 13:49 2025
Created in PyCharm
Created as QGP_Scripts/lyons_dca_toy

@author: Dylan Neff, dn277127
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from math import sqrt


def main():
    # --- Simulation parameters (you can change these) ---
    np.random.seed(42)  # reproducible
    true_value = 1.0  # the "true" underlying value the measurement approaches
    beta = 0.3  # bias coefficient (measurement = true_value + beta * dca)
    k = 3.0  # controls how fast number of events rises with dca (saturation)
    n_per_dca = 500  # maximum number of events at large dca cut
    sigma_per_event = 1.0  # per-event spread (will give stat uncertainty ~ sigma/sqrt(N))

    # DCA cuts to simulate (6 values, including a small one near 0)
    # dca_cuts = np.array([0.01, 0.03, 0.07, 0.15, 0.35, 0.8])
    dca_cuts = np.array([0.01, 1.0, 1.5, 2.0, 3.0, 4.0])

    # choose a default dca cut to highlight (user asked for default + list)
    default_dca = 2.0

    # Helper functions
    def mu_of_dca(dca):
        # measurement bias increases with dca (simple linear bias here)
        return true_value + beta * dca

    def f_dca(d):
        # event density (PDF)
        # return n_per_dca * k * np.exp(-k * d)
        return n_per_dca + (d * 0)  # flat distribution for simplicity

    def n_of_dca(dca_cut):
        # cumulative number of events up to dca_cut, integrate f_dca from 0 to dca_cut
        num, _ = quad(lambda d: f_dca(d), 0, dca_cut)
        return num

    def weighted_mu_avg(dca_cut):
        # numerator = ∫ mu(d) f(d) dd
        num, _ = quad(lambda d: mu_of_dca(d) * f_dca(d), 0, dca_cut)
        # denominator = ∫ f(d) dd = N_expected
        den = n_of_dca(dca_cut)
        return num / den

    dca_plt = np.linspace(np.min(dca_cuts), np.max(dca_cuts), 100)
    # n_dca_plt = np.vectorize(n_of_dca)(dca_plt)
    # f_dca_plt = np.vectorize(weighted_mu_avg)(dca_plt)
    n_dca_plt = np.array([n_of_dca(d) for d in dca_plt])
    f_dca_plt = np.array([f_dca(d) for d in dca_plt])
    fig_stats, ax_stats = plt.subplots()
    ax_cum = ax_stats.twinx()
    ax_stats.plot(dca_plt, f_dca_plt, color='blue', label='f(d)')
    ax_cum.plot(dca_plt, n_dca_plt, color='orange', label='N(d)')
    ax_cum.set_ylabel('Cumulative events N(d)', color='orange')
    ax_stats.set_xlabel('DCA')
    ax_stats.set_ylabel('Event density f(d)')
    ax_stats.set_title('Event density vs DCA')
    ax_stats.legend(loc='upper left')
    ax_cum.legend(loc='upper right')
    fig_stats.tight_layout()

    fig_mu, ax_mu = plt.subplots()
    mu_plt = np.array([mu_of_dca(d) for d in dca_plt])
    wmu_plt = np.array([weighted_mu_avg(d) for d in dca_plt])
    ax_mu.plot(dca_plt, mu_plt, label='mu(d)')
    ax_mu.plot(dca_plt, wmu_plt, label='weighted mu avg')
    ax_mu.set_xlabel('DCA')
    ax_mu.set_ylabel('Measurement value')
    ax_mu.set_title('Measurement bias vs DCA')
    ax_mu.axhline(true_value, color='black', linestyle='--', label='true value')
    ax_mu.legend()
    fig_mu.tight_layout()

    rows = []
    for d in dca_cuts:
        n = n_of_dca(d)
        mu_avg = weighted_mu_avg(d)
        stat_uncertainty = sigma_per_event / sqrt(n)
        observed = np.random.normal(loc=mu_avg, scale=stat_uncertainty)
        rows.append({
            "dca_cut": d,
            "n": n,
            "mu_true": mu_avg,  # now the integrated weighted mean
            "observed": observed,
            "stat_uncertainty": stat_uncertainty
        })

    #
    # # Simulate measurement for each cut
    # rows = []
    # for d in dca_cuts:
    #     N_expected = N_of_dca(d)
    #     # ensure at least 1 event
    #     N = max(1, int(round(N_expected)))
    #     mu = mu_of_dca(d)
    #     stat_uncertainty = sigma_per_event / sqrt(N)  # sigma / sqrt(N)
    #     observed = np.random.normal(loc=mu, scale=stat_uncertainty)
    #     rows.append({
    #         "dca_cut": d,
    #         "N_expected": N_expected,
    #         "N_used": N,
    #         "mu_true": mu,
    #         "observed": observed,
    #         "stat_uncertainty": stat_uncertainty
    #     })

    # Also compute the "default" measurement (sampled again at default_dca for demonstration)
    n_def = n_of_dca(default_dca)
    mu_def = weighted_mu_avg(default_dca)
    stat_unc_def = sigma_per_event / sqrt(n_def)
    obs_def = np.random.normal(loc=mu_def, scale=stat_unc_def)

    df = pd.DataFrame(rows)

    # Display table to user with nicely formatted numbers
    df_display = df.copy()
    df_display["dca_cut"] = df_display["dca_cut"].map(lambda x: f"{x:.3f}")
    df_display["n"] = df_display["n"].astype(int)
    df_display["mu_true"] = df_display["mu_true"].map(lambda x: f"{x:.4f}")
    df_display["observed"] = df_display["observed"].map(lambda x: f"{x:.4f}")
    df_display["stat_uncertainty"] = df_display["stat_uncertainty"].map(lambda x: f"{x:.5f}")

    # Plot 1: observed measurement with stat error bars vs dca_cut
    plt.figure(figsize=(6, 4))
    d = df["dca_cut"].values
    obs = df["observed"].values
    err = df["stat_uncertainty"].values
    plt.errorbar(d, obs, yerr=err, marker='o', linestyle='-')
    plt.axhline(true_value, color='black', linestyle='--', label='true value')
    plt.axvline(default_dca, color='red', linestyle=':', label='default dca cut')
    plt.legend()
    plt.xlabel("DCA cut")
    plt.ylabel("Observed measurement")
    plt.title("Observed measurement ± statistical uncertainty vs DCA cut")
    plt.grid(True)
    plt.tight_layout()

    # Plot 2: N_expected and N_used vs dca_cut
    plt.figure(figsize=(6, 4))
    plt.plot(d, df["n"].values, marker='o', linestyle='-')
    plt.xlabel("DCA cut")
    plt.ylabel("Number of events (expected and used)")
    plt.title("Event counts vs DCA cut")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Print default measurement
    print(
        f"Default dca = {default_dca:.3f} -> N_used = {n_def}, mu_true = {mu_def:.4f}, observed = {obs_def:.4f}, stat_unc = {stat_unc_def:.5f}")

    # Save CSV for later use
    # csv_path = "/mnt/data/simulated_dca_measurements.csv"
    # df.to_csv(csv_path, index=False)
    # print(f"\nSaved full numeric results to: {csv_path}")

    print('donzo')


if __name__ == '__main__':
    main()
