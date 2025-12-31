import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer, TrainingArguments, Trainer
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
from .dataset import prepare_dataset
import os

def compute_metrics(pred):
    labels = pred.label_ids
    preds = pred.predictions.argmax(-1)
    precision, recall, f1, _ = precision_recall_fscore_support(labels, preds, average='binary')
    acc = accuracy_score(labels, preds)
    return {
        'accuracy': acc,
        'f1': f1,
        'precision': precision,
        'recall': recall
    }

class ModelFinetuner:
    def __init__(self, model_name="seyonec/ChemBERTa-zinc-base-v1", device="cpu"):
        self.model_name = model_name
        self.device = device
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        
    def train(self, train_csv, output_dir, epochs=3, batch_size=16, learning_rate=2e-5, val_csv=None):
        print(f"--- Starting Fine-tuning for {self.model_name} ---")
        
        # Load Model (For Classification - 2 labels: Active/Inactive)
        model = AutoModelForSequenceClassification.from_pretrained(self.model_name, num_labels=2)
        model.to(self.device)
        
        # Prepare Data
        print("Loading Training Data...")
        train_dataset = prepare_dataset(train_csv, self.tokenizer)
        
        eval_dataset = None
        if val_csv:
            print("Loading Validation Data...")
            eval_dataset = prepare_dataset(val_csv, self.tokenizer)
            
        # Training Arguments
        training_args = TrainingArguments(
            output_dir=output_dir,
            num_train_epochs=epochs,
            per_device_train_batch_size=batch_size,
            per_device_eval_batch_size=batch_size,
            warmup_steps=500,
            weight_decay=0.01,
            logging_dir=os.path.join(output_dir, 'logs'),
            logging_steps=10,
            evaluation_strategy="epoch" if eval_dataset else "no",
            save_strategy="epoch",
            learning_rate=learning_rate,
            load_best_model_at_end=True if eval_dataset else False,
            use_mps_device=(self.device == 'mps')
        )
        
        trainer = Trainer(
            model=model,
            args=training_args,
            train_dataset=train_dataset,
            eval_dataset=eval_dataset,
            compute_metrics=compute_metrics
        )
        
        # Train
        print("Training...")
        trainer.train()
        
        # Save
        print(f"Saving model to {output_dir}...")
        model.save_pretrained(output_dir)
        self.tokenizer.save_pretrained(output_dir)
        print("Fine-tuning Complete!")
