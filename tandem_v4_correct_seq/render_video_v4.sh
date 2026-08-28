#!/usr/bin/env bash
# PyMOL/ffmpeg render of the corrected-sequence OpenMM run.
# Every quality setting is run_nvt.py's DOCUMENTED default, not a guess:
#   ray tracing ON, 1920x1440, 12 fps, x264 CRF 16, preset slow,
#   align_states on, zoom_target=context buffer -8.0 A, qm_residue CR2.
# 1920x1440 @ 12 fps matches the committed production movie recorded in
# coupling_nvt_production_cr2_1000_20260721/README.md.
# niced so the gas-phase STEOM job (submission blocker) keeps CPU priority.
set -u
export PATH=/home/robson/anaconda3/envs/TeraChem/bin:$PATH
cd /home/robson/PetaChem || exit 1
S=tandem_v4_correct_seq/render_status.txt
echo "[$(date '+%F %H:%M:%S')] render start" > $S
nice -n 15 /home/robson/anaconda3/envs/TeraChem/bin/python - <<'PY' >> tandem_v4_correct_seq/render.log 2>&1
import sys
from pathlib import Path
sys.path.insert(0, '/home/robson/PetaChem')
sys.argv = ['run_nvt.py']
import run_nvt
code = run_nvt.render_video(
    traj_path=Path('/home/robson/PetaChem/tandem_v4_correct_seq/tandem_v4_traj_300.pdb'),
    output_path=Path('/home/robson/PetaChem/videos/tandem_v4_correctseq_1ns.mp4'),
    frames_dir=Path('/home/robson/PetaChem/tandem_v4_correct_seq/render_frames'),
    fps=12, width=1920, height=1440,
    align_states=True, zoom_target='context', zoom_buffer_a=-8.0, ray=True,
    qm_residue='CR2', qm_protein_cutoff_a=2.65, qm_nearest_waters=5,
    rotate_90='none', crf=16, preset='slow',
    keep_frames=False, allow_gif_fallback=False, pymol_launch='headless',
)
print('render_video returned', code)
PY
echo "[$(date '+%F %H:%M:%S')] render exit rc=$?" >> $S
