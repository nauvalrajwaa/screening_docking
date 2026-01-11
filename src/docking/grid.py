from Bio.PDB import PDBParser
import numpy as np
import warnings

def calculate_center_from_residues(pdb_path, active_residues_str):
    """
    Calculates the geometric center of specific residues in a PDB/PDBQT file.
    
    Args:
        pdb_path (str): Path to PDB/PDBQT file.
        active_residues_str (str): Comma-separated list of residues, e.g. "A:41,A:145"
                                   Format: Chain:ResNum
                                   
    Returns:
        tuple: (x, y, z) coordinates of the center.
    """
    try:
        # PDBQT is PDB-compatible enough for PDBParser generally (ignores partial charges col)
        # Suppress PDB warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', pdb_path)
            
        target_residues = []
        # Parse residue string: "A:41,B:20" -> [('A', 41), ('B', 20)]
        for item in active_residues_str.split(','):
            parts = item.strip().split(':')
            if len(parts) == 2:
                chain_id = parts[0]
                try:
                    res_id = int(parts[1])
                    target_residues.append((chain_id, res_id))
                except ValueError:
                    print(f"Warning: Invalid residue number in '{item}'. Skipping.")
            else:
                 print(f"Warning: Invalid residue format '{item}'. Expected Chain:ResNum (e.g., A:41).")

        atom_coords = []
        found_residues = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check if this residue is in our target list
                    # residue.id is usually (' ', resseq, ' ')
                    res_seq = residue.id[1]
                    chain_id = chain.id
                    
                    if (chain_id, res_seq) in target_residues:
                        found_residues += 1
                        # Use all atoms or CA? Let's use all atoms for better pocket center
                        for atom in residue:
                            atom_coords.append(atom.get_coord())
                            
        if not atom_coords:
            print("Error: No matching residues found in structure.")
            return None
            
        # Calculate centroid
        coords_array = np.array(atom_coords)
        center = np.mean(coords_array, axis=0)
        
        print(f"Calculated center based on {found_residues} residues ({len(atom_coords)} atoms).")
        return tuple(center)

    except Exception as e:
        print(f"Error calculating grid center: {e}")
        return None
