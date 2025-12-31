import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
import numpy as np

class LLMDescriptor:
    # Pre-defined models for easy access
    AVAILABLE_MODELS = {
        "chemberta-base": "seyonec/ChemBERTa-zinc-base-v1",
        "chemberta-77m": "seyonec/ChemBERTa-zinc-deepchem-77m",
        "chemberta-mtr": "DeepChem/ChemBERTa-77M-MTR",
        "chemberta-mlm": "DeepChem/ChemBERTa-77M-MLM",
        "molbert": "samsmasling/MolBERT" # Example, verify if exists or use another known one
    }

    def __init__(self, model_name="seyonec/ChemBERTa-zinc-base-v1", device="cpu"):
        self.device = torch.device(device)
        
        # Check if model_name is a short alias
        if model_name.lower() in self.AVAILABLE_MODELS:
            self.model_name = self.AVAILABLE_MODELS[model_name.lower()]
            print(f"--- Model alias '{model_name}' resolved to '{self.model_name}' ---")
        else:
            self.model_name = model_name
            
        print(f"--- Loading AI Model ({self.model_name}) on {self.device}... ---")
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
            self.model = AutoModel.from_pretrained(self.model_name).to(self.device)
        except Exception as e:
            print(f"Error loading model {model_name}: {e}")
            raise e

    def get_embedding(self, smiles):
        """Generates embedding for a single SMILES string"""
        if not smiles or not isinstance(smiles, str):
            return None
        return self.get_batch_embeddings([smiles])[0]

    def get_batch_embeddings(self, smiles_list, batch_size=32):
        """Generates embeddings for a list of SMILES strings in batches"""
        embeddings = []
        
        for i in range(0, len(smiles_list), batch_size):
            batch = smiles_list[i:i+batch_size]
            
            # Tokenize batch
            inputs = self.tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=512)
            inputs = {k: v.to(self.device) for k, v in inputs.items()}
            
            with torch.no_grad():
                out = self.model(**inputs)
            
            # Mean pooling and move to CPU
            batch_embeddings = out.last_hidden_state.mean(dim=1).cpu().numpy()
            embeddings.extend(batch_embeddings)
            
        return embeddings

    def calculate_similarity(self, emb1, emb2, metric='cosine'):
        """Calculates similarity/distance between two embeddings"""
        if emb1 is None or emb2 is None:
            return 0.0
            
        if metric == 'cosine':
            return cosine_similarity(emb1, emb2)[0][0]
        elif metric == 'euclidean':
            # Return negative distance so higher is better (closer) or just raw distance
            # Usually for distance, lower is better. Here we return raw distance.
            return euclidean_distances(emb1, emb2)[0][0]
        else:
            raise ValueError(f"Unknown metric: {metric}")
