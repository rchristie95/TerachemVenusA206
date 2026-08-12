import json
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import os

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
    
    with open(os.path.join(base_dir, 'results/hi_camb3lyp_frame_0000/site_A/energy_summary.json')) as f:
        data_a = json.load(f)
    with open(os.path.join(base_dir, 'results/hi_camb3lyp_frame_0000/site_B/energy_summary.json')) as f:
        data_b = json.load(f)
        
    root_a = data_a['roots'][0]
    root_b = data_b['roots'][0]
    
    e1 = root_a['energy_cm-1']
    e2 = root_b['energy_cm-1']
    mu1 = np.array(root_a['transition_dipole_au'])
    mu2 = np.array(root_b['transition_dipole_au'])
    
    triple = 100.0 # Arbitrary scaling for plotting
    
    J = 32.8 
    omega = 1700.0
    S = 0.5
    n_max = 5
    
    H, n_vib = build_holstein_hamiltonian(e1, e2, J, omega, S, n_max)
    evals, evecs = scipy.linalg.eigh(H)
    
    grid = np.linspace(min(e1, e2) - 1000, max(e1, e2) + 6000, 1000)
    abs_spec = np.zeros_like(grid)
    cd_spec = np.zeros_like(grid)
    hwhm = 400.0
    
    idx_1_00 = 0 * n_vib * n_vib + 0 * n_vib + 0
    idx_2_00 = 1 * n_vib * n_vib + 0 * n_vib + 0
    
    for k in range(len(evals)):
        energy = evals[k]
        c1 = evecs[idx_1_00, k]
        c2 = evecs[idx_2_00, k]
        
        exciton_mu = c1 * mu1 + c2 * mu2
        dip_strength = np.dot(exciton_mu, exciton_mu)
        rot_strength = -np.pi * energy * c1 * c2 * triple
        
        shape = lorentzian(grid, energy, hwhm)
        abs_spec += dip_strength * shape
        cd_spec += rot_strength * shape
        
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(1e7/grid, abs_spec)
    plt.title('Holstein Absorption')
    plt.xlabel('Wavelength (nm)')
    
    plt.subplot(1, 2, 2)
    plt.plot(1e7/grid, cd_spec)
    plt.axhline(0, color='k', linestyle='--')
    plt.title('Holstein CD (Relative)')
    plt.xlabel('Wavelength (nm)')
    plt.tight_layout()
    plot_path = os.path.join(base_dir, 'holstein_frame0000.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Saved plot to {plot_path}")
