import os
import json
import numpy as np
import scipy.linalg
from scipy.fft import fft, fftfreq, fftshift

# =====================================================================
# Point 1: Vibronic Structure (Holstein Exciton Model)
# =====================================================================
def build_holstein_hamiltonian(e1, e2, J, omega, S, n_max):
    """
    Build the 2-site Holstein Hamiltonian in the localized multi-quanta basis.
    Basis: |el_site, v_1, v_2>
    
    e1, e2: Site energies (cm^-1)
    J: Electronic coupling (cm^-1)
    omega: Mode frequency (cm^-1)
    S: Huang-Rhys factor
    n_max: Max vibrational quanta per site
    
    Returns: H matrix and basis list
    """
    n_vib = n_max + 1
    dim = 2 * n_vib * n_vib
    H = np.zeros((dim, dim))
    
    # Pre-calculate Franck-Condon overlap integrals between displaced and undisplaced oscillators
    # For a displacement d = sqrt(2S), the overlap <m|n_disp> can be calculated recursively or analytically.
    # We will use a simplified linear displacement in the second-quantization form instead:
    # H = sum_i E_i |i><i| + J (|1><2| + |2><1|) + omega * sum_j b_j^+ b_j 
    #     + omega * sqrt(S) * sum_i |i><i| (b_i^+ + b_i)
    
    # Basis: el in {0, 1}, v1, v2
    for el in range(2):
        for v1 in range(n_vib):
            for v2 in range(n_vib):
                i = el * n_vib * n_vib + v1 * n_vib + v2
                
                # Diagonal terms
                H[i, i] = (e1 if el == 0 else e2) + omega * (v1 + v2)
                
                # Vibrational coupling on the excited site
                v_exc = v1 if el == 0 else v2
                
                if v_exc + 1 < n_vib:
                    # b^\dagger term
                    j = el * n_vib * n_vib + (v1+1 if el==0 else v1) * n_vib + (v2+1 if el==1 else v2)
                    H[i, j] = omega * np.sqrt(S) * np.sqrt(v_exc + 1)
                    H[j, i] = H[i, j] # Hermitian conjugate (b term)
                    
                # Exciton transfer (J)
                # J only transfers the electronic excitation, leaving vibrational quanta unchanged
                # because the basis |v1, v2> is the *undisplaced* basis for both states.
                other_el = 1 - el
                j_transfer = other_el * n_vib * n_vib + v1 * n_vib + v2
                if el == 0:
                    H[i, j_transfer] = J
                    H[j_transfer, i] = J
                    
    return H

def get_holstein_spectrum(e1, e2, J, omega, S, n_max=4, kT=200):
    H = build_holstein_hamiltonian(e1, e2, J, omega, S, n_max)
    evals, evecs = scipy.linalg.eigh(H)
    return evals, evecs


# =====================================================================
# Point 2: Electrostatic Descriptor Surrogate 
# =====================================================================
def calc_electrostatic_descriptor(qm_coords, qm_charges_gs, qm_charges_es, mm_coords, mm_charges):
    """
    Calculate the first-order detuning surrogate: sum_i (delta q_i) * V_i
    where V_i is the electrostatic potential from MM point charges at QM atom i.
    """
    # delta_q = q_es - q_gs
    delta_q = qm_charges_es - qm_charges_gs
    
    V_qm = np.zeros(len(qm_coords))
    for i, qm_pos in enumerate(qm_coords):
        dist = np.linalg.norm(mm_coords - qm_pos, axis=1)
        # Convert to Bohr/au if necessary, assuming coordinates are in Angstrom and charges in e
        # 1 e / Angstrom = 14.3996 eV
        V_qm[i] = np.sum(14.3996 * mm_charges / dist)
        
    detuning = np.sum(delta_q * V_qm)
    return detuning

def fit_descriptor_to_qm(qm_detunings, descriptor_values):
    """
    Fit descriptor to QM detunings to obtain a scaling factor.
    detuning_qm = a * descriptor + b
    """
    A = np.vstack([descriptor_values, np.ones(len(descriptor_values))]).T
    a, b = np.linalg.lstsq(A, qm_detunings, rcond=None)[0]
    return a, b


# =====================================================================
# Point 3: Second-Order Cumulant Lineshape
# =====================================================================
def cumulant_lineshape(C_t, dt, T, max_t):
    """
    Calculate the lineshape function g(t) from the energy gap correlation function C(t).
    g(t) = \\int_0^t dt' \\int_0^{t'} dt'' C(t'')
    """
    times = np.arange(0, max_t, dt)
    # Double cumulative sum (trapezoidal approximation)
    C_t_sampled = np.interp(times, np.arange(len(C_t))*dt, C_t)
    
    int1 = np.cumsum(C_t_sampled) * dt
    g_t = np.cumsum(int1) * dt
    
    # The dipole auto-correlation is exp(-g(t))
    dipole_corr = np.exp(-g_t)
    
    # Spectrum is real part of Fourier transform
    freqs = fftfreq(len(times), d=dt)
    spectrum = np.real(fft(dipole_corr))
    
    return fftshift(freqs), fftshift(spectrum)


# =====================================================================
# Point 4: Polarizable Screening
# =====================================================================
def apply_optical_screening(detunings, eps_opt=1.78):
    """
    Apply empirical optical dielectric screening to the detunings.
    """
    return detunings / eps_opt


# =====================================================================
# Point 5: Detuning Scale Calibration (STEOM)
# =====================================================================
def apply_steom_calibration(tddft_detunings, steom_detunings):
    """
    Calibrate TDDFT detunings against STEOM detunings.
    """
    a, b = np.linalg.lstsq(np.vstack([tddft_detunings, np.ones(len(tddft_detunings))]).T, steom_detunings, rcond=None)[0]
    return a, b


# =====================================================================
# Point 6: Absolute Molar Ellipticity (ORCA Input Generation)
# =====================================================================
def generate_orca_steom_input(coords, symbols, mm_charges=None):
    """
    Generate ORCA input for STEOM-CCSD to get transition magnetic dipoles.
    """
    lines = [
        "! DLPNO-STEOM-CCSD def2-SVP def2/J def2-SVP/C",
        "%mdci",
        "  nroots 4",
        "  dograd true # to get densities/dipoles",
        "end",
        "* xyz 0 1"
    ]
    for sym, coord in zip(symbols, coords):
        lines.append(f"{sym} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
    lines.append("*")
    
    if mm_charges is not None:
        lines.append("%pointcharges \"charges.pc\"")
        
    return "\n".join(lines)

if __name__ == "__main__":
    print("Quantitative CD Toolkit Module loaded.")
    
    # Test Point 1
    print("Testing Point 1: Holstein Hamiltonian")
    H = build_holstein_hamiltonian(e1=20000, e2=20100, J=50, omega=1700, S=0.5, n_max=3)
    evals, _ = scipy.linalg.eigh(H)
    print(f"Top eigenvalues: {evals[:5]}")
    
    # Test Point 2 stub
    print("Testing Point 2 stub: Electrostatic descriptor")
    qm = np.array([[0,0,0], [1,0,0]])
    mm = np.array([[5,0,0]])
    dq = np.array([-0.1, 0.1])
    mm_q = np.array([1.0])
    det = calc_electrostatic_descriptor(qm, dq, dq, mm, mm_q)
    print(f"Calculated surrogate detuning: {det:.2f} eV")
