from Bio.PDB import PDBList, PDBParser, MMCIFParser
import os
import shutil
import warnings
import sys
from rdkit import Chem

def download_pdb(pdb_id, output_dir="data/pdb"):
    """
    Downloads a PDB file.
    Args:
        pdb_id (str): PDB ID.
        output_dir (str): Directory.
    Returns:
        str: Path to downloaded file or None.
    """
    pdb_id = pdb_id.lower()
    os.makedirs(output_dir, exist_ok=True)
    pdbl = PDBList(verbose=False)
    try:
        ent_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=output_dir, overwrite=True)
        # Rename .ent to .pdb for consistency if needed, or keep as is.
        # BioPython saves as pdb<id>.ent
        final_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        if os.path.exists(ent_file):
             shutil.move(ent_file, final_path)
             return final_path
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        return None

def fetch_ligand_from_pdb(pdb_id, output_dir="data/controls_pdb"):
    """
    Downloads a PDB structure, extracts HETATM ligands, and converts to SMILES.
    
    Args:
        pdb_id (str): PDB ID (e.g., '3HTB')
        output_dir (str): Directory to save downloaded files.
        
    Returns:
        list: List of dictionaries [{'name': '3HTB_LIG', 'smiles': '...', 'ligand_id': '...'}, ...]
    """
    pdb_id = pdb_id.lower()
    
    # 1. Download PDB
    ent_file = download_pdb(pdb_id, output_dir)
    if not ent_file:
         return []

    # 2. Parse and Find Ligands
    extracted_ligands = []
    
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, ent_file)

        # Iterate over all residues to find HETATMs (excluding water)
        ligands_found = {} # Key: (ResName, Chain, ID) -> Residue Object
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    # HETATM check
                    # residue.id[0] starts with 'H_' for hetero residues
                    if residue.id[0].startswith('H_'):
                        res_name = residue.resname
                        if res_name not in ['HOH', 'WAT']: # Skip water
                            # Store unique ligands
                            # We might want to save them as separate PDB files to convert to SMILES
                            key = (res_name, chain.id, residue.id[1])
                            if key not in ligands_found:
                                ligands_found[key] = residue

        # 3. Save Ligand to PDB and Convert to SMILES
        # We need to save the ligand atoms to a temporary PDB file, then use RDKit/Obabel to read it.
        
        from Bio.PDB import PDBIO, Select
        
        class ResSelect(Select):
            def __init__(self, target_residue):
                self.target_residue = target_residue
            def accept_residue(self, residue):
                return residue == self.target_residue

        io = PDBIO()
        io.set_structure(structure)
        
        for (res_name, chain, res_id), residue in ligands_found.items():
            lig_filename = f"{pdb_id}_{res_name}_{chain}_{res_id}.pdb"
            lig_path = os.path.join(output_dir, lig_filename)
            
            io.save(lig_path, ResSelect(residue))
            
            # Convert to SMILES using RDKit
            try:
                mol = Chem.MolFromPDBFile(lig_path)
                if mol:
                    smi = Chem.MolToSmiles(mol)
                    if smi:
                        extracted_ligands.append({
                            'name': f"{pdb_id.upper()}_{res_name}",
                            'smiles': smi,
                            'origin': f"PDB:{pdb_id}"
                        })
            except Exception as e:
                print(f"Failed to convert {lig_filename} to SMILES: {e}")
                
    except Exception as e:
        print(f"Error extracting ligands from {pdb_id}: {e}")
        
    return extracted_ligands
