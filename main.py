import argparse
import pandas as pd
import os
import sys
from tqdm import tqdm
from datetime import datetime

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.utils.data_loader import load_csv, clean_data
from src.utils.report import generate_html_report
from src.descriptors.chemical import get_mol, calculate_ecfp, calculate_bro5, calculate_tanimoto, calculate_dice
from src.descriptors.llm import LLMDescriptor
from src.screening.deeppurpose_module import DeepPurposeScreener
from src.docking import VinaDocker, AutoDockDocker, convert_smiles_to_pdbqt, calculate_center_from_residues, prepare_receptor
from src.utils.pdb import fetch_ligand_from_pdb, download_pdb
import re

class Logger(object):
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        self.terminal.flush()
        self.log.flush()

def main():
    parser = argparse.ArgumentParser(description="LLM Prediction & Screening Tool")
    # Input Files
    parser.add_argument('--compounds', type=str, required=True, help="Path to CSV file containing compounds")
    parser.add_argument('--controls', type=str, help="Path to CSV file containing control compounds (optional if --pdb_controls is set)")
    parser.add_argument('--output', type=str, default="experiment_results", help="Base name for output files")
    parser.add_argument('--model', type=str, default="seyonec/ChemBERta-zinc-base-v1", 
                        help="HuggingFace model name or alias (chemberta-base, chemberta-77m, chemberta-mtr, chemberta-mlm)")
    parser.add_argument('--device', type=str, default="cpu", choices=['cpu', 'cuda', 'mps'], help="Device to use (cpu, cuda, mps)")
    parser.add_argument('--use_deeppurpose', action='store_true', help="Enable DeepPurpose screening")
    parser.add_argument('--dp_target', type=str, help="Target sequence (Amino Acid) for DeepPurpose DTI screening")
    parser.add_argument('--dp_model', type=str, default="MPNN_CNN_BindingDB", help="DeepPurpose pre-trained model name")
    parser.add_argument('--pdb_controls', type=str, help="Comma-separated PDB IDs to fetch ligands from (e.g., '3HTB,1ABC')")
    
    # Docking Arguments
    parser.add_argument('--docking_mode', type=str, default='none', choices=['none', 'vina', 'autodock'], help="Docking mode")
    parser.add_argument('--top_n_dock', type=int, default=10, help="Number of top candidates to dock")
    parser.add_argument('--receptor', type=str, help="Path to receptor PDBQT/Maps")
    parser.add_argument('--center_x', type=float, help="Grid center X")
    parser.add_argument('--center_y', type=float, help="Grid center Y")
    parser.add_argument('--center_z', type=float, help="Grid center Z")
    parser.add_argument('--size_x', type=float, default=20, help="Grid size X")
    parser.add_argument('--size_y', type=float, default=20, help="Grid size Y")
    parser.add_argument('--size_z', type=float, default=20, help="Grid size Z")
    parser.add_argument('--vina_bin', type=str, default='vina', help="Path to Vina binary")
    parser.add_argument('--autodock_bin', type=str, default='autodock-gpu', help="Path to AutoDock-GPU binary")
    parser.add_argument('--active_residues', type=str, help="Comma-separated active residues (e.g., 'A:41,A:145') to center grid")
    
    args = parser.parse_args()
    
    # --- Smart Receptor Logic ---
    # Check if receptor is a PDB ID and handle it EARLY so it can be used for controls
    if args.receptor and re.match(r'^[a-zA-Z0-9]{4}$', args.receptor) and not os.path.exists(args.receptor):
        pdb_id = args.receptor
        print(f"Detected PDB ID '{pdb_id}' as receptor. Downloading...")
        downloaded_path = download_pdb(pdb_id, output_dir="data/receptors")
        
        if downloaded_path:
            print(f"Receptor downloaded to: {downloaded_path}")
            args.receptor = downloaded_path # Update arg to file path
            
            # Auto-sync with pdb_controls if empty
            if not args.pdb_controls:
                print(f"Auto-adding '{pdb_id}' to --pdb_controls for screening comparison.")
                args.pdb_controls = pdb_id
            elif pdb_id not in args.pdb_controls:
                # Optional: Add it anyway? User said "automate sync".
                # If user specified SOME controls, maybe they want those AND the receptor?
                # For safety, let's append it.
                print(f"Adding '{pdb_id}' to existing --pdb_controls.")
                args.pdb_controls += f",{pdb_id}"
        else:
            print(f"Warning: Failed to download PDB {pdb_id}. Will assume it is a file path that is missing.")

    # 0. Setup Run Directory
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join("runs", run_id)
    os.makedirs(run_dir, exist_ok=True)
    
    # Setup Logger
    log_path = os.path.join(run_dir, "run.log")
    sys.stdout = Logger(log_path)
    
    print(f"--- Run ID: {run_id} | Output Directory: {run_dir} ---")
    print(f"Command: {' '.join(sys.argv)}")

    # 1. Load Data
    print("--- Loading Data ---")
    try:
        df_compounds = load_csv(args.compounds)
        
        df_compounds = clean_data(df_compounds)
        
        print(f"Compounds: {len(df_compounds)}")
    except Exception as e:
        print(f"Error loading compounds data: {e}")
        return

    # 2. Initialize Models
    print(f"--- Initializing Models (Device: {args.device}) ---")
    llm_desc = LLMDescriptor(model_name=args.model, device=args.device)
    
    dp_screener = None
    if args.use_deeppurpose:
        dp_screener = DeepPurposeScreener(target_seq=args.dp_target, model_name=args.dp_model)

    # 3. Process Controls
    print("--- Pre-processing Controls ---")
    control_data = []
    
    # 3a. CSV Controls (Optional)
    if args.controls:
        df_control = load_csv(args.controls)
        df_control = clean_data(df_control)
        
        for i, row in df_control.iterrows():
            smi = row['smiles']
            name = row.get('nama_kontrol', row.get('nama', f'Ctrl_{i+1}'))
            mol = get_mol(smi)
            if mol:
                fp = calculate_ecfp(mol)
                emb = llm_desc.get_embedding(smi)
                control_data.append({
                    'name': name,
                    'smiles': smi,
                    'fp': fp,
                    'emb': emb
                })

    # 3b. Fetch PDB Controls
    if args.pdb_controls:
        print(f"--- Fetching PDB Controls: {args.pdb_controls} ---")
        pdb_ids = [pid.strip() for pid in args.pdb_controls.split(',')]
        for pid in pdb_ids:
            fetched_ligands = fetch_ligand_from_pdb(pid)
            for lig in fetched_ligands:
                print(f"  + Added control: {lig['name']}")
                mol = get_mol(lig['smiles'])
                if mol:
                    fp = calculate_ecfp(mol)
                    emb = llm_desc.get_embedding(lig['smiles'])
                    control_data.append({
                        'name': lig['name'],
                        'smiles': lig['smiles'],
                        'fp': fp,
                        'emb': emb
                    })

    if len(control_data) == 0:
        print("Error: No controls provided. Please specify --controls or --pdb_controls.")
        return

    print(f"Valid Controls: {len(control_data)}")

    # 4. Process Compounds (Batch Processing)
    print("--- Processing Compounds ---")
    results = []
    
    # Prepare batches
    batch_size = 32
    compounds_list = df_compounds.to_dict('records')
    
    for i in tqdm(range(0, len(compounds_list), batch_size), desc="Processing Batches"):
        batch = compounds_list[i:i+batch_size]
        batch_smiles = [item['smiles'] for item in batch]
        
        # 1. Batch LLM Embeddings (GPU Accelerated)
        batch_embeddings = llm_desc.get_batch_embeddings(batch_smiles, batch_size=len(batch))
        
        # 2. Process each compound in the batch
        for j, item in enumerate(batch):
            smi_comp = item['smiles']
            mol_comp = get_mol(smi_comp)
            
            if not mol_comp:
                continue
                
            # Chemical Descriptors
            props = calculate_bro5(mol_comp) # Now uses the upgraded calculate_properties
            fp_comp = calculate_ecfp(mol_comp)
            emb_comp = batch_embeddings[j]
            
            entry = {
                'SMILES': smi_comp,
                **props
            }
            
            # Similarity to Controls
            for ctrl in control_data:
                # 1. Chemical Similarities
                tanimoto = calculate_tanimoto(fp_comp, ctrl['fp'])
                dice = calculate_dice(fp_comp, ctrl['fp'])
                
                # 2. LLM Similarities
                cosine = llm_desc.calculate_similarity(emb_comp, ctrl['emb'], metric='cosine')
                euclidean = llm_desc.calculate_similarity(emb_comp, ctrl['emb'], metric='euclidean')
                
                # 3. Hybrid Score (Weighted Average)
                # Default: 50% Chemical (Tanimoto) + 50% LLM (Cosine)
                hybrid_score = (0.5 * tanimoto) + (0.5 * cosine)

                entry[f"Tanimoto_{ctrl['name']}"] = round(tanimoto, 3)
                entry[f"Dice_{ctrl['name']}"] = round(dice, 3)
                entry[f"LLM_Cos_{ctrl['name']}"] = round(cosine, 3)
                entry[f"LLM_Euc_{ctrl['name']}"] = round(euclidean, 3)
                entry[f"Hybrid_{ctrl['name']}"] = round(hybrid_score, 3)
            
            # DeepPurpose Screening (Optional)
            if dp_screener and dp_screener.available:
                # Predict property/affinity
                pred_score = dp_screener.predict(smi_comp)
                entry['DeepPurpose_Score'] = pred_score

            results.append(entry)

    # 5. Save Results
    df_results = pd.DataFrame(results)
    
    # Calculate Max Hybrid Score for sorting
    hybrid_cols = [c for c in df_results.columns if c.startswith('Hybrid_')]
    if hybrid_cols:
        df_results['Max_Hybrid_Score'] = df_results[hybrid_cols].max(axis=1)
        df_results = df_results.sort_values(by='Max_Hybrid_Score', ascending=False)
        
        print("\n" + "="*60)
        print("TOP 5 CANDIDATES (Based on Hybrid Score)")
        print("="*60)
        cols_to_show = ['SMILES', 'bRo5_Status', 'Max_Hybrid_Score'] + hybrid_cols[:2] # Show first 2 controls
        print(df_results[cols_to_show].head(5).to_string(index=False))
        print("="*60 + "\n")

    # 6. Docking (Optional)
    if args.docking_mode != 'none':
        if not args.receptor:
            print("Error: Docking mode requires --receptor argument.")
            return

        receptor_file = args.receptor
        center_calc_file = args.receptor # Default to whatever the user provided
        
        # Auto-Receptor Preparation
        if receptor_file.endswith('.pdb'):
            center_calc_file = receptor_file # Explicitly use the PDB for center calc
            print(f"Detected .pdb receptor. Converting to .pdbqt...")
            pdbqt_file = receptor_file.replace('.pdb', '.pdbqt')
            if prepare_receptor(receptor_file, pdbqt_file):
                receptor_file = pdbqt_file
            else:
                print("Failed to convert receptor. Proceeding with original file check...")
        
        if not os.path.exists(receptor_file):
            print(f"Error: Receptor file {receptor_file} not found.")
            return
            
        print(f"--- Starting Docking ({args.docking_mode}) on Top {args.top_n_dock} Candidates ---")
        
        # Setup Docker
        docker = None
        if args.docking_mode == 'vina':
            docker = VinaDocker(args.vina_bin)
        elif args.docking_mode == 'autodock':
            docker = AutoDockDocker(args.autodock_bin)
            
        # Filter top N
        docking_candidates = df_results.head(args.top_n_dock).copy()
            
        # Determine Grid Center
        center = (args.center_x, args.center_y, args.center_z)
        if args.active_residues:
            print(f"Calculating grid center from residues: {args.active_residues}")
            center = calculate_center_from_residues(center_calc_file, args.active_residues) # Use PDB if available
            if center:
                print(f"Calculated Center: {center}")
            else:
                print("Error: Could not calculate center from residues. Falling back to manual coordinates.")
                center = (args.center_x, args.center_y, args.center_z)

        docking_scores = []
        
        docking_dir = os.path.join(run_dir, "docking")
        os.makedirs(docking_dir, exist_ok=True)
        
        for idx, row in tqdm(docking_candidates.iterrows(), total=len(docking_candidates), desc="Docking"):
            smi = row['SMILES']
            name = f"cand_{idx}"
            ligand_pdbqt = os.path.join(docking_dir, f"{name}.pdbqt")
            output_pdbqt = os.path.join(docking_dir, f"{name}_out.pdbqt")
            
            # Convert to PDBQT
            if convert_smiles_to_pdbqt(smi, ligand_pdbqt):
                # Dock
                size = (args.size_x, args.size_y, args.size_z)
                    
                if center[0] is None:
                    print("Error: Grid center coordinates required.")
                    break

                score = docker.dock(
                    ligand_pdbqt, 
                    receptor_file, 
                    center, 
                    size, 
                    output_pdbqt
                )
                docking_scores.append(score)
            else:
                docking_scores.append(None)
                    
            # Add scores back to results
            # Note: This aligns because we iterated over the slice copy. 
            # To merge back robustly to the main dataframe:
            docking_candidates['Docking_Score'] = docking_scores
            
            # Update original dataframe
            # We explicitly update the rows that were docked
            df_results.loc[docking_candidates.index, 'Docking_Score'] = docking_scores
            
            # Re-sort if we have docking scores? Maybe just leave as is but output new top list
            print("\n" + "="*60)
            print("TOP DOCKING RESULTS")
            print("="*60)
            print(df_results.dropna(subset=['Docking_Score']).sort_values(by='Docking_Score')[['SMILES', 'Max_Hybrid_Score', 'Docking_Score']].to_string(index=False))

    # CSV
    csv_path = os.path.join(run_dir, f"{args.output}.csv")

    df_results.to_csv(csv_path, index=False)
    print(f"Results saved to {csv_path}")
    
    # HTML Report
    html_path = os.path.join(run_dir, f"{args.output}.html")
    generate_html_report(df_results, output_path=html_path)

if __name__ == "__main__":
    main()
