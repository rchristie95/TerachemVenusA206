import numpy as np

def parse_pdb_atoms(pdb_file):
    atoms = []
    lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line)
            lines.append(line)
    return atoms, lines

def parse_biomt(pdb_file):
    transforms = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            # Look for SMTRY (Crystallographic) or BIOMT (Biological)
            # Prioritize SMTRY to explore crystal packing
            if "REMARK 290   SMTRY" in line:
                parts = line.split()
                # Format: REMARK 290 SMTRYn idx ...
                row = int(parts[2][-1]) - 1
                idx = int(parts[3])
                
                if idx not in transforms:
                    transforms[idx] = np.eye(4)
                
                vals = [float(x) for x in parts[4:8]]
                transforms[idx][row, :] = vals
    return transforms

def parse_scale_matrix(pdb_file):
    S = np.eye(4)
    found = False
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("SCALE1"):
                vals = [float(x) for x in line.split()[1:5]]
                S[0, :] = vals
                found = True
            elif line.startswith("SCALE2"):
                vals = [float(x) for x in line.split()[1:5]]
                S[1, :] = vals
            elif line.startswith("SCALE3"):
                vals = [float(x) for x in line.split()[1:5]]
                S[2, :] = vals
    return S if found else None

def get_centroid(atoms):
    coords = []
    for line in atoms:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])
    return np.mean(coords, axis=0)

def apply_transform(atom_line, matrix, chain_id):
    # Extract XYZ
    x = float(atom_line[30:38])
    y = float(atom_line[38:46])
    z = float(atom_line[46:54])
    
    vec = np.array([x, y, z, 1.0])
    new_vec = np.dot(matrix, vec)
    
    # Format new line
    # Preserve formatting strictly
    new_line = list(atom_line)
    
    # Update XYZ
    new_x = f"{new_vec[0]:8.3f}"
    new_y = f"{new_vec[1]:8.3f}"
    new_z = f"{new_vec[2]:8.3f}"
    
    new_line[30:38] = list(new_x)
    new_line[38:46] = list(new_y)
    new_line[46:54] = list(new_z)
    
    # Update Chain ID
    new_line[21] = chain_id
    
    return "".join(new_line)

def main():
    pdb_file = "1MYW.pdb"
    out_file = "venus_dimer.pdb"
    
    atoms, raw_lines = parse_pdb_atoms(pdb_file)
    transforms = parse_biomt(pdb_file)
    scale_mat = parse_scale_matrix(pdb_file)
    
    if not transforms:
        print("[!] No transforms found.")
        return

    ref_center = get_centroid(atoms)
    ref_center_hom = np.array([ref_center[0], ref_center[1], ref_center[2], 1.0])
    
    # Invert Scale Matrix to go Frac -> Cart
    # S maps Cart -> Frac. S_inv maps Frac -> Cart.
    inv_scale = np.eye(4)
    if scale_mat is not None:
        inv_scale = np.linalg.inv(scale_mat)
        print("Using SCALE matrix for lattice search.")
        
    closest_dist = 9999.0
    best_M = None
    best_shift = None # (dx, dy, dz) in Cartesian
    
    for idx, M in transforms.items():
        # 1. Apply Symmetry in Cartesian
        p_sym = np.dot(M, ref_center_hom)
        
        if scale_mat is not None:
            # 2. Convert to Fractional
            p_frac = np.dot(scale_mat, p_sym)
            
            # 3. Search lattice shifts (i,j,k)
            for i in range(-3, 4):
                for j in range(-3, 4):
                    for k in range(-3, 4):
                        shift_frac = np.array([float(i), float(j), float(k), 0.0])
                        candidate_frac = p_frac + shift_frac
                        
                        # Convert back to Cartesian
                        candidate_cart = np.dot(inv_scale, candidate_frac)
                        
                        dist = np.linalg.norm(candidate_cart[:3] - ref_center)
                        
                        if dist > 0.1 and dist < closest_dist:
                            closest_dist = dist
                            best_M = M
                            shift_cart_vec = np.dot(inv_scale, shift_frac)
                            best_shift = shift_cart_vec
        else:
            # No lattice info, just check raw transform
            dist = np.linalg.norm(p_sym[:3] - ref_center)
            if dist > 0.1 and dist < closest_dist:
                closest_dist = dist
                best_M = M
                best_shift = np.zeros(4)
                
    print(f"Closest Neighbor Distance: {closest_dist:.2f} A")
    
    if best_M is not None:
        with open(out_file, 'w') as f:
            f.write("HEADER    GENERATED CLOSEST DIMER WITH LATTICE SHIFT\\n")
            
            # Chain A
            for line in atoms:
                f.write(apply_transform(line, np.eye(4), 'A'))
                
            # Chain B
            effective_M = best_M.copy()
            effective_M[:, 3] += best_shift
            
            for line in atoms:
                f.write(apply_transform(line, effective_M, 'B'))
                
            f.write("END\\n")
        print(f"Wrote {out_file}")
    else:
        print("Could not find a valid neighbor.")

if __name__ == "__main__":
    main()