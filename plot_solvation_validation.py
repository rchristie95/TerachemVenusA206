#!/usr/bin/env python3
"""Plot the tandem-Venus solvation/dephasing numerical validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


LIGHT_CM_PER_S = 2.99792458e10


def correlation(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float) - np.mean(values)
    n = len(values)
    return np.correlate(values, values, mode="full")[n - 1 :] / np.arange(n, 0, -1)


def cumulative_trapezoid(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros(len(values), dtype=float)
    out[1:] = np.cumsum(0.5 * (values[1:] + values[:-1]) * np.diff(x))
    return out


def curves(a: np.ndarray, b: np.ndarray, dt_fs: float):
    difference = np.asarray(a) - np.asarray(b)
    corr = correlation(difference)
    normalized = corr / corr[0]
    time_fs = np.arange(len(corr)) * dt_fs
    time_s = time_fs * 1.0e-15
    c_omega = (2.0 * np.pi * LIGHT_CM_PER_S) ** 2 * corr
    g = cumulative_trapezoid(cumulative_trapezoid(c_omega, time_s), time_s)
    return time_fs, normalized, np.exp(-np.clip(g, 0.0, 700.0))


def first_crossing(time_fs: np.ndarray, values: np.ndarray, threshold: float) -> float:
    index = int(np.flatnonzero(values <= threshold)[0])
    if index == 0:
        return float(time_fs[0])
    x0, x1 = time_fs[index - 1 : index + 1]
    y0, y1 = values[index - 1 : index + 1]
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap-inputs", type=Path, nargs=2, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--pme", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    parts = [np.load(path) for path in args.gap_inputs]
    energy = np.concatenate([part["energies_cm"] for part in parts])
    dt_fs = float(np.median(np.diff(parts[0]["time_fs"])))
    total = energy.sum(axis=2)
    protein = energy[:, :, 0]
    water = energy[:, :, 1]
    time_fs = np.arange(len(energy)) * dt_fs
    lag, total_corr, total_coh = curves(total[:, 0], total[:, 1], dt_fs)
    _, protein_corr, _ = curves(protein[:, 0], protein[:, 1], dt_fs)
    _, water_corr, water_coh = curves(water[:, 0], water[:, 1], dt_fs)

    pme_data = np.load(args.pme)
    pme_energy = pme_data["pme_cm"]
    pme_dt = float(np.median(np.diff(pme_data["time_fs"])))
    pme_lag, pme_corr, pme_coh = curves(pme_energy[:, 0], pme_energy[:, 1], pme_dt)
    validation = json.loads(args.validation.read_text(encoding="utf-8"))

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5))
    ax = axes[0, 0]
    differential = total[:, 0] - total[:, 1]
    ax.plot(time_fs / 1000.0, differential - differential.mean(), lw=0.8, color="0.15")
    ax.set(xlabel="trajectory time (ps)", ylabel=r"$\delta(E_A-E_B)$ (cm$^{-1}$)")
    ax.set_title("A  Differential electrostatic gap")

    ax = axes[0, 1]
    ax.plot(lag, total_corr, label="all MM, 7-pair probe", lw=1.8)
    ax.plot(lag, protein_corr, label="protein only", lw=1.2)
    ax.plot(lag, water_corr, label="water only", lw=1.2)
    ax.plot(pme_lag, pme_corr, label="full PME, CR2 probe", lw=1.5, ls="--")
    ax.axhline(1 / np.e, color="0.55", ls=":", lw=1)
    ax.axhline(0, color="0.75", lw=0.7)
    ax.set(xlim=(0, 500), ylim=(-0.35, 1.03), xlabel="lag (fs)", ylabel="normalized correlation")
    ax.set_title("B  Bath memory is multicomponent")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 0]
    ax.plot(lag, total_coh, label="7-pair minimum image", lw=2)
    ax.plot(pme_lag, pme_coh, label="CR2 full PME", lw=1.7, ls="--")
    ax.plot(lag, water_coh, label="water only", lw=1.3)
    ax.plot(lag, np.exp(-lag / 60.0), label="assumed 60 fs exponential", lw=1.2, color="0.45", ls=":")
    ax.axhline(1 / np.e, color="tab:red", ls=":", lw=1)
    ax.set(xlim=(0, 180), ylim=(0, 1.02), xlabel="time (fs)", ylabel=r"coherence $|L(t)|$")
    ax.set_title("C  Solvation memory is not $T_2^*$")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 1]
    block = validation["block_metrics"]
    block_t2 = np.asarray([row["total_with_ions"]["classical_cumulant_T2_1e_fs"] for row in block])
    labels = ["seg. 1", "seg. 2", "8 ps", "PME", "water", "assumed"]
    values = [
        validation["segments"][0]["total_with_ions"]["classical_cumulant_T2_1e_fs"],
        validation["segments"][1]["total_with_ions"]["classical_cumulant_T2_1e_fs"],
        validation["full_trajectory"]["total_with_ions"]["classical_cumulant_T2_1e_fs"],
        first_crossing(pme_lag, pme_coh, 1 / np.e),
        validation["full_trajectory"]["water"]["classical_cumulant_T2_1e_fs"],
        60.0,
    ]
    colors = ["tab:blue", "tab:blue", "tab:blue", "tab:green", "tab:cyan", "0.55"]
    ax.bar(np.arange(len(values)), values, color=colors, alpha=0.85)
    ax.errorbar(
        2,
        values[2],
        yerr=np.std(block_t2, ddof=1),
        fmt="none",
        ecolor="black",
        capsize=4,
        lw=1,
        label="1 ps block SD",
    )
    ax.set_xticks(np.arange(len(labels)), labels, rotation=25, ha="right")
    ax.set(ylabel=r"$T_2^*$ at $|L|=e^{-1}$ (fs)", ylim=(0, 85))
    ax.set_title("D  Dephasing estimate and checks")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Tandem Venus: electrostatic solvation and decoherence test", fontsize=13)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
