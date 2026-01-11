import subprocess
import os

def convert_smiles_to_pdbqt(smiles, output_path, name="ligand"):
    """
    Converts a SMILES string to a PDBQT file using OpenBabel.
    1. Generate 3D coordinates (--gen3d).
    2. Add hydrogens (-h).
    3. Protonate at pH 7.4 (-p 7.4).
    """
    try:
        # Construct obabel command
        # obabel -:"<smiles>" -O <output_path> --gen3d -h -p 7.4
        cmd = [
            "obabel",
            f"-:{smiles}",
            "-O", output_path,
            "--gen3d",
            "-h",
            "-p", "7.4"
        ]
        
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        print(f"Error converting {name}: {e}")
        return False

def prepare_receptor(pdb_file, output_pdbqt):
    """
    Converts a PDB receptor file to PDBQT format using OpenBabel.
    Adds polar hydrogens and handles partial charges (standard Obabel behavior for receptors).
    
    Args:
        pdb_file (str): Path to input PDB file.
        output_pdbqt (str): Path to output PDBQT file.
    
    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        # obabel -ipdb <input> -opdbqt -O <output> -xr -xn -xp
        # -xr: output rigid fragments (not really needed for receptor but safe)
        # -xn: ignore waters (usually good for docking unless water is critical)
        # -xp: add polar hydrogens
        # We generally want to add hydrogens (-h) if they are missing.
        
        cmd = [
            "obabel",
            "-ipdb", pdb_file,
            "-opdbqt",
            "-O", output_pdbqt,
            "-xn", # No waters
            "-h",  # Add hydrogens
            "-p", "7.4" # Protonate residues
        ]
        
        print(f"Converting receptor {pdb_file} to PDBQT...")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # Post-processing: Remove ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF lines
        if os.path.exists(output_pdbqt):
            with open(output_pdbqt, 'r') as f:
                lines = f.readlines()
            
            with open(output_pdbqt, 'w') as f:
                for line in lines:
                    if not line.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                        f.write(line)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Obabel error converting receptor: {e}")
        return False
    except Exception as e:
        print(f"Error preparing receptor: {e}")
        return False
