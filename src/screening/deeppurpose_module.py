import warnings
import pandas as pd
import numpy as np

try:
    from DeepPurpose import utils, CompoundPred
    DEEPPURPOSE_AVAILABLE = True
except ImportError:
    DEEPPURPOSE_AVAILABLE = False

class DeepPurposeScreener:
    def __init__(self, target_seq=None, model_name='MPNN_CNN_BindingDB'):
        self.available = DEEPPURPOSE_AVAILABLE
        self.target_seq = target_seq
        self.model_name = model_name
        self.model = None
        
        if not self.available:
            print("Warning: DeepPurpose library not found. Screening will be skipped.")
            return
            
        print(f"--- Initializing DeepPurpose ({model_name}) ---")
        try:
            from DeepPurpose import models
            # Load pre-trained model
            # Note: This will download the model if not present (~100MB+)
            if 'BindingDB' in model_name or 'DAVIS' in model_name or 'KIBA' in model_name:
                self.model = models.model_pretrained(model=model_name)
                self.mode = 'DTI'
                if not self.target_seq:
                    print("Warning: DTI model loaded but no target sequence provided via --dp_target.")
            else:
                # Assume ADMET or other property prediction
                # DeepPurpose has specific functions for ADMET, this is a general placeholder
                self.mode = 'ADMET'
                print("Note: ADMET screening requires specific model loading logic.")
                
        except Exception as e:
            print(f"Error loading DeepPurpose model: {e}")
            self.available = False

    def predict(self, smiles):
        """
        Runs prediction for a single SMILES string.
        """
        if not self.available or not self.model:
            return "N/A"
            
        try:
            from DeepPurpose import utils
            
            if self.mode == 'DTI':
                if not self.target_seq:
                    return "No Target"
                
                # Data Process for DTI
                X_pred = utils.data_process(X_drug=[smiles], X_target=[self.target_seq], y=[0], 
                                          drug_encoding=self.model.drug_encoding, 
                                          target_encoding=self.model.target_encoding, 
                                          split_method='no_split',
                                          frac=[0,0,0])
                
                # Predict
                y_pred = self.model.predict(X_pred)
                return round(y_pred[0], 3)
                
            else:
                return "N/A (ADMET Not Impl)"
                
        except Exception as e:
            print(f"DeepPurpose Error: {e}")
            return "Error"
