import json
import os
import glob
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

def build_holstein_hamiltonian(e1, e2, J, omega, S, n_max):
    n_vib = n_max + 1
    dim = 2 * n_vib * n_vib
    H = np.zeros((dim, dim))
    for el in range(2):
        for v1 in range(n_vib):
            for v2 in range(n_vib):
                i = el * n_vib * n_vib + v1 * n_vib + v2
                H[i, i] = (e1 if el == 0 else e2) + omega * (v1 + v2)
                v_exc = v1 if el == 0 else v2
                if v_exc + 1 < n_vib:
                    j = el * n_vib * n_vib + (v1+1 if el==0 else v1) * n_vib + (v2+1 if el==1 else v2)
                    H[i, j] = omega * np.sqrt(S) * np.sqrt(v_exc + 1)
                    H[j, i] = H[i, j]
                other_el = 1 - el
                j_transfer = other_el * n_vib * n_vib + v1 * n_vib + v2
                if el == 0:
                    H[i, j_transfer] = J
                    H[j_transfer, i] = J
    return H, n_vib

def lorentzian(grid, center, hwhm):
    return hwhm / np.pi / ((grid - center)**2 + hwhm**2)

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(base_dir, 'results')
    
    frame_dirs = sorted(glob.glob(os.path.join(results_dir, 'hi_camb3lyp_frame_*')))
    
    grid = np.linspace(15000, 30000, 2000)
    avg_abs_spec = np.zeros_like(grid)
    avg_cd_spec = np.zeros_like(grid)
    
    hwhm = 400.0
    J = 32.8 
    omega = 1700.0
    S = 0.5
    n_max = 4
    
    # We will use an effective geometry triple product to mimic the CD orientation.
    # From investigation, typical rotational strength scaling
    triple = 20.0 
    
    valid_frames = 0
    for frame_dir in frame_dirs:
        try:
            with open(os.path.join(frame_dir, 'site_A', 'energy_summary.json')) as f:
                data_a = json.load(f)
            with open(os.path.join(frame_dir, 'site_B', 'energy_summary.json')) as f:
                data_b = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            continue
            
        root_a = data_a['roots'][0]
        root_b = data_b['roots'][0]
        
        e1 = root_a['energy_cm-1']
        e2 = root_b['energy_cm-1']
        mu1 = np.array(root_a['transition_dipole_au'])
        mu2 = np.array(root_b['transition_dipole_au'])
        
        # We need a proper geometric triple product. 
        # r2 - r1 cross (mu1, mu2). 
        # For simplicity since we don't parse the PDB here, we'll use a fixed approximate triple.
        # But wait! We can compute the actual triple if we read the geometry.xyz centroids!
        geom_a = np.loadtxt(os.path.join(frame_dir, 'site_A', 'geometry.xyz'), skiprows=2, usecols=(1,2,3))
        geom_b = np.loadtxt(os.path.join(frame_dir, 'site_B', 'geometry.xyz'), skiprows=2, usecols=(1,2,3))
        cent_a = np.mean(geom_a, axis=0)
        cent_b = np.mean(geom_b, axis=0)
        
        triple_actual = float(np.dot(cent_b - cent_a, np.cross(mu1, mu2)))
        
        H, n_vib = build_holstein_hamiltonian(e1, e2, J, omega, S, n_max)
        evals, evecs = scipy.linalg.eigh(H)
        
        idx_1_00 = 0
        idx_2_00 = 1 * n_vib * n_vib
        
        frame_abs = np.zeros_like(grid)
        frame_cd = np.zeros_like(grid)
        
        for k in range(len(evals)):
            energy = evals[k]
            c1 = evecs[idx_1_00, k]
            c2 = evecs[idx_2_00, k]
            
            exciton_mu = c1 * mu1 + c2 * mu2
            dip_strength = np.dot(exciton_mu, exciton_mu)
            rot_strength = -np.pi * energy * c1 * c2 * triple_actual
            
            shape = lorentzian(grid, energy, hwhm)
            frame_abs += dip_strength * shape
            frame_cd += rot_strength * shape
            
        avg_abs_spec += frame_abs
        avg_cd_spec += frame_cd
        valid_frames += 1

    if valid_frames > 0:
        avg_abs_spec /= valid_frames
        avg_cd_spec /= valid_frames
        
        # Plotting against nm
        wavelengths = 1e7 / grid
        
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(wavelengths, avg_abs_spec, lw=2, color='C0')
        plt.title(f'Ensemble Holstein Absorption (n={valid_frames})')
        plt.xlabel('Wavelength (nm)')
        plt.xlim(350, 700)
        
        plt.subplot(1, 2, 2)
        plt.plot(wavelengths, avg_cd_spec, lw=2, color='C1')
        plt.axhline(0, color='k', linestyle='-', lw=1, alpha=0.5)
        plt.title(f'Ensemble Holstein CD (n={valid_frames})')
        plt.xlabel('Wavelength (nm)')
        plt.xlim(350, 700)
        
        plt.tight_layout()
        plot_path = os.path.join(base_dir, 'ensemble_holstein.png')
        plt.savefig(plot_path, dpi=300)
        print(f"Successfully processed {valid_frames} frames.")
        print(f"Saved ensemble results to {plot_path}")
    else:
        print("No valid frames found.")
