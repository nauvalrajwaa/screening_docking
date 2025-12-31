import argparse
import pandas as pd
import os
import sys
from tqdm import tqdm

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.utils.data_loader import load_csv, clean_data
from src.utils.report import generate_html_report
from src.descriptors.chemical import get_mol, calculate_ecfp, calculate_bro5, calculate_tanimoto
from src.descriptors.llm import LLMDescriptor
from src.screening.deeppurpose_module import DeepPurposeScreener

def main():
    parser = argparse.ArgumentParser(description="LLM Prediction & Screening Tool")
    parser.add_argument('--compounds', type=str, required=True, help="Path to compounds CSV")
    parser.add_argument('--controls', type=str, required=True, help="Path to controls CSV")
    parser.add_argument('--model', type=str, default="seyonec/ChemBERTa-zinc-base-v1", 
                        help="HuggingFace model name or alias (chemberta-base, chemberta-77m, chemberta-mtr, chemberta-mlm)")
    parser.add_argument('--device', type=str, default="cpu", choices=['cpu', 'cuda', 'mps'], help="Device to use (cpu, cuda, mps)")
    parser.add_argument('--output', type=str, default="results", help="Output filename prefix")
    parser.add_argument('--use_deeppurpose', action='store_true', help="Enable DeepPurpose screening")
    
    args = parser.parse_args()

    # 1. Load Data
    print("--- Loading Data ---")
    try:
        df_compounds = load_csv(args.compounds)
        df_control = load_csv(args.controls)
        
        df_compounds = clean_data(df_compounds)
        df_control = clean_data(df_control)
        
        print(f"Compounds: {len(df_compounds)}")
        print(f"Controls: {len(df_control)}")
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # 2. Initialize Models
    print(f"--- Initializing Models (Device: {args.device}) ---")
    llm_desc = LLMDescriptor(model_name=args.model, device=args.device)
    
    dp_screener = None
    if args.use_deeppurpose:
        dp_screener = DeepPurposeScreener()

    # 3. Pre-process Controls
    print("--- Pre-processing Controls ---")
    control_data = []
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
                # Tanimoto (Chemical)
                tanimoto = calculate_tanimoto(fp_comp, ctrl['fp'])
                # Cosine (LLM)
                cosine = llm_desc.calculate_similarity(emb_comp, ctrl['emb'])
                
                entry[f"Tanimoto_{ctrl['name']}"] = round(tanimoto, 3)
                entry[f"LLM_{ctrl['name']}"] = round(cosine, 3)
            
            # DeepPurpose Screening (Optional)
            if dp_screener and dp_screener.available:
                # Placeholder for actual screening logic
                preds = dp_screener.predict_properties([smi_comp])
                entry['DeepPurpose_Pred'] = preds[0]

            results.append(entry)

    # 5. Save Results
    df_results = pd.DataFrame(results)
    
    # CSV
    csv_path = f"{args.output}.csv"
    df_results.to_csv(csv_path, index=False)
    print(f"Results saved to {csv_path}")
    
    # HTML Report
    html_path = f"{args.output}.html"
    generate_html_report(df_results, output_path=html_path)

if __name__ == "__main__":
    main()
