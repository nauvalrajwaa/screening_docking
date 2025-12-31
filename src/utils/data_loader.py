import pandas as pd
import os

def load_csv(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    try:
        df = pd.read_csv(filepath)
        # Standardize columns
        df.columns = df.columns.str.lower().str.strip()
        return df
    except Exception as e:
        raise Exception(f"Error reading CSV: {e}")

def clean_data(df, smiles_col='smiles'):
    if smiles_col not in df.columns:
        raise ValueError(f"Column '{smiles_col}' not found in dataframe")
    
    df = df.dropna(subset=[smiles_col])
    df[smiles_col] = df[smiles_col].astype(str)
    return df
