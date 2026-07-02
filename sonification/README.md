# Sonifications

Short **sonified videos** of the paper's results for the
*openquantumsonification* channel. In every clip the audio is *physically*
driven by the same quantities shown on screen — the sound is a rendering of the
data, not a soundtrack laid over it.

Regenerate everything with:

```bash
python sonify.py            # build all four clips into this folder
python sonify.py 1 3        # build only clips 1 and 3
```

Needs `numpy`, `scipy`, `matplotlib` and `ffmpeg` (all in `../environment.yml`);
no GPU. Audio is synthesised in pure NumPy and muxed with ffmpeg.

| Clip | You see | You hear | Driven by |
|------|---------|----------|-----------|
| `decoherence.mp4` | One stochastic (SSE) quantum trajectory: coherence, site and adiabatic populations | The excitation pans between the two chromophore "speakers"; the tone shimmers while coherent and goes dull as it localises (sub-100 fs) | site populations → L/R gain; \|ρ₁₂(t)\| → shimmer/brightness (`open_quantum_dynamics.py`) |
| `davydov_chord.mp4` | The two exciton eigenstates E±J(t) and their absorption peaks converging | Two tones whose detuning **is** the coupling; as the solvent relaxes (optical→static dielectric) the interval collapses from audible beating into a unison | J(t) = Debye-screened coupling → the interval width |
| `nvt_breathing.mp4` | The real restrained-NVT dimer trajectory (200 frames) | The carrier pitch follows the excitonic coupling J(t); broadband "thermal" noise swells with the potential energy | per-frame J (STEOM ensemble) → pitch; per-frame NVT potential energy → noise amplitude |
| `parameter_sweep.mp4` | The coherence-decay curve steepening as T₂\* scans 200→20 fs | The pitch glides up and the tremolo speeds up: faster dephasing = a shriller, more agitated tone | T₂\* sweep → pitch and tremolo rate |

All quantities come from the released pipeline (`../open_quantum_dynamics.py`,
`../coupling_paper_steom_thermal/coupling_samples.csv`, `../nvt_restrained.log`,
`../videos/nvt_restrained.mp4`).
