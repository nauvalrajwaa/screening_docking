import torch
from torch.utils.data import Dataset
import pandas as pd

class MoleculeDataset(Dataset):
    def __init__(self, encodings, labels=None):
        self.encodings = encodings
        self.labels = labels

    def __getitem__(self, idx):
        item = {key: torch.tensor(val[idx]) for key, val in self.encodings.items()}
        if self.labels is not None:
            item['labels'] = torch.tensor(self.labels[idx])
        return item

    def __len__(self):
        return len(self.encodings['input_ids'])

def prepare_dataset(csv_path, tokenizer, smiles_col='smiles', label_col='label', max_length=512):
    """
    Loads data from CSV and tokenizes it for training.
    Assumes binary classification (0/1) for now.
    """
    df = pd.read_csv(csv_path)
    
    # Basic cleaning
    df.columns = df.columns.str.lower().str.strip()
    smiles_col = smiles_col.lower()
    label_col = label_col.lower()
    
    if smiles_col not in df.columns or label_col not in df.columns:
        raise ValueError(f"CSV must contain '{smiles_col}' and '{label_col}' columns")
        
    df = df.dropna(subset=[smiles_col, label_col])
    
    smiles = df[smiles_col].astype(str).tolist()
    labels = df[label_col].astype(int).tolist()
    
    # Tokenize
    encodings = tokenizer(smiles, truncation=True, padding=True, max_length=max_length)
    
    return MoleculeDataset(encodings, labels)
