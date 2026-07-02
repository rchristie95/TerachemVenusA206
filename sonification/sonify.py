#!/usr/bin/env python3
r"""
sonify.py -- turn the paper's physical results into short sonified videos for
the *openquantumsonification* channel.

Each clip pairs a data-driven animation with an audio track synthesised so that
the sound is *physically* controlled by the same quantities on screen:

  1. decoherence.mp4    -- open-quantum-system dynamics.  A single stochastic
     (SSE) trajectory: the excitation audibly pans between the two chromophore
     "speakers" (site populations -> L/R gain) while the tone's shimmer tracks
     the coherence |rho_12(t)| and fades as the state localises (sub-100 fs).
  2. davydov_chord.mp4  -- the two exciton eigenstates split by 2|J|.  Rendered
     as two tones whose detuning IS the coupling; as the solvent relaxes
     (optical -> static dielectric over the Debye time) the coupling is screened
     and the interval collapses from audible beating into a unison.
  3. nvt_breathing.mp4  -- the real restrained-NVT dimer movie, sonified: the
     carrier pitch follows the excitonic coupling J(t) over the trajectory and
     broadband noise amplitude follows the per-frame potential energy, so you
     hear the protein "breathing" and modulating the coupling.
  4. parameter_sweep.mp4 -- the dephasing-time sweep.  As T2* scans 200 -> 20 fs
     the pitch glides up and the tremolo speeds up: faster dephasing = a shriller,
     more agitated tone, over the animated coherence-decay curve.

Pure NumPy/Matplotlib/SciPy audio synthesis (no GPU); video muxed with ffmpeg.

    python sonify.py            # build all four into this folder
    python sonify.py 1 3        # build only clips 1 and 3

Data sources (relative to the repo root, one level up):
  open_quantum_dynamics.py                          (imported: OQS solvers)
  videos/nvt_restrained.mp4                          (NVT dimer movie, 200f/10s)
  nvt_restrained.log                                 (200 per-frame NVT energies)
  coupling_paper_steom_thermal/coupling_samples.csv  (200 per-frame STEOM J)
"""

import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.io import wavfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from open_quantum_dynamics import make_params, solve_me, solve_sse, J_of_t  # noqa: E402

SR = 44100          # audio sample rate
FPS = 25            # video frame rate
SEED = 20260618

# --------------------------------------------------------------------------- #
# audio helpers
# --------------------------------------------------------------------------- #
def _resample(ctrl, n):
    """Stretch a control array to n samples by linear interpolation."""
    ctrl = np.asarray(ctrl, float)
    x = np.linspace(0.0, 1.0, len(ctrl))
    return np.interp(np.linspace(0.0, 1.0, n), x, ctrl)


def _norm01(a):
    a = np.asarray(a, float)
    lo, hi = np.min(a), np.max(a)
    return (a - lo) / (hi - lo) if hi > lo else np.zeros_like(a)


def _osc(freq_per_sample, sr=SR):
    """Sine oscillator for a per-sample instantaneous frequency array."""
    phase = 2.0 * np.pi * np.cumsum(freq_per_sample) / sr
    return np.sin(phase)


def _fade(sig, sr=SR, sec=0.06):
    n = int(sr * sec)
    if n > 0 and len(sig) > 2 * n:
        env = np.ones(len(sig))
        env[:n] = np.linspace(0, 1, n)
        env[-n:] = np.linspace(1, 0, n)
        sig = sig * env
    return sig


def _write_wav(path, left, right, sr=SR):
    L = _fade(np.asarray(left, float))
    R = _fade(np.asarray(right, float))
    peak = max(np.max(np.abs(L)), np.max(np.abs(R)), 1e-9)
    scale = 0.9 / peak
    stereo = np.stack([L * scale, R * scale], axis=1)
    wavfile.write(path, sr, (stereo * 32767).astype(np.int16))


def _save_anim(anim, silent_mp4, fps=FPS, dpi=100):
    writer = animation.FFMpegWriter(fps=fps, bitrate=2400,
                                    extra_args=["-pix_fmt", "yuv420p"])
    anim.save(silent_mp4, writer=writer, dpi=dpi)
    plt.close("all")


def _mux(silent_mp4, wav, out_mp4, reencode_video=False, scale=None):
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-i", silent_mp4, "-i", wav]
    if reencode_video:
        vf = ["-vf", f"scale={scale}"] if scale else []
        cmd += vf + ["-c:v", "libx264", "-pix_fmt", "yuv420p"]
    else:
        cmd += ["-c:v", "copy"]
    cmd += ["-c:a", "aac", "-b:a", "192k", "-shortest", out_mp4]
    subprocess.run(cmd, check=True)


# --------------------------------------------------------------------------- #
# 1. decoherence
# --------------------------------------------------------------------------- #
def make_decoherence(tmp, dur=12.0):
    p = make_params()
    me = solve_me(p)
    sse = solve_sse(p, seed=SEED)
    n = int(dur * SR)

    P1 = _resample(sse["P1"], n)
    P2 = _resample(sse["P2"], n)
    coh_s = _resample(sse["coh"], n)
    coh_m = _resample(me["coh"], n)
    coh_s = _norm01(coh_s) if coh_s.max() > 0 else coh_s

    fL, fR = 196.0, 294.0                     # G3 & D4 -> the two sites
    toneL = _osc(np.full(n, fL))
    toneR = _osc(np.full(n, fR))
    shimmerL = coh_s * _osc(np.full(n, fL * 3))
    shimmerR = coh_s * _osc(np.full(n, fR * 3))
    drone = 0.25 * coh_m * _osc(np.full(n, 98.0))   # ensemble coherence reference
    left = P1 * (toneL + 0.6 * shimmerL) + drone
    right = P2 * (toneR + 0.6 * shimmerR) + drone
    wav = os.path.join(tmp, "decoherence.wav")
    _write_wav(wav, left, right)

    # animation
    t = sse["t"]
    nfr = int(dur * FPS)
    fig, ax = plt.subplots(3, 1, figsize=(8, 6.2), sharex=True)
    fig.suptitle("Hearing decoherence: one stochastic quantum trajectory",
                 fontsize=12, weight="bold")
    ax[0].plot(t, me["coh"], color="#1f77b4", lw=2, label=r"$|\rho_{12}|$ ensemble")
    ax[0].plot(t, sse["coh"], color="#d62728", lw=0.8, alpha=0.8, label=r"$|\rho_{12}|$ single")
    ax[0].set_ylabel("coherence"); ax[0].legend(loc="upper right", fontsize=8)
    ax[1].plot(t, sse["P1"], color="#2ca02c", lw=1, label=r"$P_L$ (left)")
    ax[1].plot(t, sse["P2"], color="#9467bd", lw=1, label=r"$P_R$ (right)")
    ax[1].set_ylabel("site pop."); ax[1].legend(loc="upper right", fontsize=8)
    ax[2].plot(t, me["PB"], color="#ff7f0e", lw=2, label="bright")
    ax[2].plot(t, me["PD"], color="#8c564b", lw=2, label="dark")
    ax[2].set_ylabel("adiabatic"); ax[2].set_xlabel("time (ps)")
    ax[2].legend(loc="right", fontsize=8)
    cursors = [a.axvline(0, color="k", lw=1.2) for a in ax]
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    def upd(i):
        tc = (i / max(nfr - 1, 1)) * t[-1]
        for c in cursors:
            c.set_xdata([tc, tc])
        return cursors

    anim = animation.FuncAnimation(fig, upd, frames=nfr, blit=False)
    silent = os.path.join(tmp, "decoherence_silent.mp4")
    _save_anim(anim, silent)
    _mux(silent, wav, os.path.join(HERE, "decoherence.mp4"))


# --------------------------------------------------------------------------- #
# 2. davydov chord
# --------------------------------------------------------------------------- #
def make_davydov(tmp, dur=12.0, t_end_ps=40.0):
    p = make_params()
    n = int(dur * SR)
    tt = np.linspace(0.0, t_end_ps, n)
    J = J_of_t(tt, p)                       # cm^-1, relaxes 74.4 -> ~1.7
    f0, K = 330.0, 0.15
    g = K * J                                # Hz detuning (audible beating)
    f1 = _osc(f0 + g)
    f2 = _osc(f0 - g)
    mono = 0.5 * (f1 + f2)
    wav = os.path.join(tmp, "davydov.wav")
    _write_wav(wav, mono, mono)

    nfr = int(dur * FPS)
    tt_fr = np.linspace(0.0, t_end_ps, nfr)
    J_fr = J_of_t(tt_fr, p)
    E0 = 0.0
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.2))
    fig.suptitle(r"The Davydov 'chord': solvent screening collapses $2|J|$",
                 fontsize=12, weight="bold")
    ax[0].plot(tt_fr, E0 + J_fr, color="#d62728", lw=2, label=r"$E+J(t)$ (bright)")
    ax[0].plot(tt_fr, E0 - J_fr, color="#1f77b4", lw=2, label=r"$E-J(t)$ (dark)")
    ax[0].set_ylabel(r"exciton energy (cm$^{-1}$)")
    ax[0].set_xlabel("time (ps)"); ax[0].legend(loc="upper right", fontsize=8)
    lev_dots = ax[0].plot([0, 0], [J_fr[0], -J_fr[0]], "ko", ms=6)[0]

    x = np.linspace(-90, 90, 600)

    def lorentz(x, x0, w=6.0):
        return (w**2) / ((x - x0)**2 + w**2)

    (specline,) = ax[1].plot(x, lorentz(x, J_fr[0]) + lorentz(x, -J_fr[0]),
                             color="#6a3d9a", lw=2)
    ax[1].set_ylabel("absorption (arb.)"); ax[1].set_xlabel(r"detuning (cm$^{-1}$)")
    ax[1].set_ylim(0, 2.2)
    txt = ax[1].text(0.02, 0.9, "", transform=ax[1].transAxes, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    def upd(i):
        j = J_fr[i]
        lev_dots.set_data([tt_fr[i], tt_fr[i]], [j, -j])
        specline.set_ydata(lorentz(x, j) + lorentz(x, -j))
        txt.set_text(rf"$2|J| = {2*j:5.1f}\ \mathrm{{cm}}^{{-1}}$")
        return lev_dots, specline, txt

    anim = animation.FuncAnimation(fig, upd, frames=nfr, blit=False)
    silent = os.path.join(tmp, "davydov_silent.mp4")
    _save_anim(anim, silent)
    _mux(silent, wav, os.path.join(HERE, "davydov_chord.mp4"))


# --------------------------------------------------------------------------- #
# 3. NVT breathing (real dimer movie + energy/coupling-driven audio)
# --------------------------------------------------------------------------- #
def _load_nvt_energy():
    path = os.path.join(ROOT, "nvt_restrained.log")
    E = []
    with open(path) as fh:
        for line in fh:
            if "NVT] frame=" in line and "energy=" in line:
                E.append(float(line.split("energy=")[1].split()[0]))
    return np.array(E)


def _load_coupling_J():
    path = os.path.join(ROOT, "coupling_paper_steom_thermal", "coupling_samples.csv")
    J = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            J.append(float(row["J_cm"]))
    return np.array(J)


def make_nvt_breathing(tmp, dur=10.0):
    movie = os.path.join(ROOT, "videos", "nvt_restrained.mp4")
    if not os.path.exists(movie):
        print("  [skip] videos/nvt_restrained.mp4 not found")
        return
    E = _load_nvt_energy()
    J = _load_coupling_J()
    n = int(dur * SR)

    Jn = _norm01(_resample(J, n))
    En = _norm01(_resample(E, n))
    freq = 190.0 + Jn * (380.0 - 190.0)      # pitch tracks the coupling J(t)
    carrier = 0.5 * _osc(freq)
    rng = np.random.default_rng(SEED)
    noise = rng.standard_normal(n)
    # gentle low-pass on the noise so it reads as "thermal rustle", not hiss
    k = 40
    noise = np.convolve(noise, np.ones(k) / k, mode="same")
    noise_amp = 0.08 + 0.55 * En            # amplitude tracks potential energy
    mix = carrier + noise_amp * noise
    wav = os.path.join(tmp, "nvt.wav")
    _write_wav(wav, mix, mix)
    # mux straight onto the real dimer movie (downscaled to keep the file small)
    _mux(movie, wav, os.path.join(HERE, "nvt_breathing.mp4"),
         reencode_video=True, scale="640:480")


# --------------------------------------------------------------------------- #
# 4. dephasing-time sweep
# --------------------------------------------------------------------------- #
def make_sweep(tmp, dur=10.0):
    n = int(dur * SR)
    frac = np.linspace(0.0, 1.0, n)
    T2 = 0.200 * (0.020 / 0.200) ** frac     # 200 -> 20 fs, log sweep (ps)
    freq = 200.0 * (0.10 / T2)               # faster dephasing -> higher pitch
    trem_rate = 2.0 + (1.0 / T2) / 3.0       # tremolo speeds up as T2 shrinks
    trem = 0.5 * (1.0 + 0.9 * np.sin(2 * np.pi * np.cumsum(trem_rate) / SR))
    tone = _osc(freq) * trem
    wav = os.path.join(tmp, "sweep.wav")
    _write_wav(wav, tone, tone)

    nfr = int(dur * FPS)
    T2_fr = 0.200 * (0.020 / 0.200) ** np.linspace(0, 1, nfr)
    tt = np.linspace(0, 0.6, 400)            # ps
    fig, ax = plt.subplots(figsize=(8, 5.2))
    fig.suptitle(r"Sweeping the dephasing time $T_2^*$: 200 $\to$ 20 fs",
                 fontsize=12, weight="bold")
    for T2v in [0.200, 0.120, 0.060, 0.030, 0.020]:
        ax.plot(tt, np.exp(-tt / T2v), color="0.8", lw=1)
    (curve,) = ax.plot(tt, np.exp(-tt / T2_fr[0]), color="#d62728", lw=2.5)
    txt = ax.text(0.6, 0.85, "", fontsize=11)
    ax.set_xlabel("time (ps)"); ax.set_ylabel(r"coherence $\propto e^{-t/T_2^*}$")
    ax.set_xlim(0, 0.6); ax.set_ylim(0, 1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    def upd(i):
        curve.set_ydata(np.exp(-tt / T2_fr[i]))
        txt.set_text(rf"$T_2^* = {T2_fr[i]*1000:5.0f}$ fs")
        return curve, txt

    anim = animation.FuncAnimation(fig, upd, frames=nfr, blit=False)
    silent = os.path.join(tmp, "sweep_silent.mp4")
    _save_anim(anim, silent)
    _mux(silent, wav, os.path.join(HERE, "parameter_sweep.mp4"))


BUILDERS = {
    "1": ("decoherence", make_decoherence),
    "2": ("davydov_chord", make_davydov),
    "3": ("nvt_breathing", make_nvt_breathing),
    "4": ("parameter_sweep", make_sweep),
}


def main(argv):
    which = argv or list(BUILDERS)
    with tempfile.TemporaryDirectory() as tmp:
        for key in which:
            name, fn = BUILDERS[key]
            print(f"[sonify] building {name}.mp4 ...")
            fn(tmp)
            print(f"[sonify] wrote {name}.mp4")


if __name__ == "__main__":
    main(sys.argv[1:])
