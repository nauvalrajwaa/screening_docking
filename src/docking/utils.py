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
    Adds polar hydrogens and Gasteiger partial charges (required for AutoGrid4).
    
    Args:
        pdb_file (str): Path to input PDB file.
        output_pdbqt (str): Path to output PDBQT file.
    
    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        cmd = [
            "obabel",
            "-ipdb", pdb_file,
            "-opdbqt",
            "-O", output_pdbqt,
            "-xn",  # No waters
            "-h",   # Add hydrogens
            "-p", "7.4",  # Protonate residues
            "--partialcharge", "gasteiger"  # Add Gasteiger charges for AutoGrid4
        ]
        
        print(f"Converting receptor {pdb_file} to PDBQT with Gasteiger charges...")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # Post-processing: Remove ROOT/ENDROOT/BRANCH/ENDBRANCH/TORSDOF lines
        if os.path.exists(output_pdbqt):
            with open(output_pdbqt, 'r') as f:
                lines = f.readlines()
            
            with open(output_pdbqt, 'w') as f:
                for line in lines:
                    if not any(keyword in line for keyword in ['ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF']):
                        f.write(line)
            
            print(f"Receptor prepared: {output_pdbqt}")
            return True
        else:
            print(f"Failed to generate {output_pdbqt}")
            return False
            
    except Exception as e:
        print(f"Error preparing receptor {pdb_file}: {e}")
        return False
