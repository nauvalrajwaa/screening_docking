import argparse
import os
import sys

# Add src to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.finetuning.trainer import ModelFinetuner

def main():
    parser = argparse.ArgumentParser(description="Fine-tune LLM for Specific Target (e.g., MDR)")
    parser.add_argument('--train_data', type=str, required=True, help="Path to training CSV (must have 'smiles' and 'label' columns)")
    parser.add_argument('--val_data', type=str, help="Path to validation CSV (optional)")
    parser.add_argument('--model', type=str, default="seyonec/ChemBERTa-zinc-base-v1", help="Base model to fine-tune")
    parser.add_argument('--output_dir', type=str, default="finetuned_model", help="Directory to save the fine-tuned model")
    parser.add_argument('--epochs', type=int, default=3, help="Number of training epochs")
    parser.add_argument('--batch_size', type=int, default=16, help="Batch size")
    parser.add_argument('--lr', type=float, default=2e-5, help="Learning rate")
    parser.add_argument('--device', type=str, default="cpu", choices=['cpu', 'cuda', 'mps'], help="Device to use")
    
    args = parser.parse_args()
    
    # Initialize Finetuner
    finetuner = ModelFinetuner(model_name=args.model, device=args.device)
    
    # Run Training
    finetuner.train(
        train_csv=args.train_data,
        val_csv=args.val_data,
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr
    )

if __name__ == "__main__":
    main()
