# LLM-Based Compound Prediction & Screening Tool

This project is a comprehensive tool for analyzing chemical compounds using both traditional chemical descriptors (ECFP, MACCS, bRo5, QED) and modern Large Language Models (ChemBERTa). It includes features for similarity searching, property calculation, DeepPurpose screening, and model fine-tuning.

## Features

*   **Hybrid Descriptors**: Combines traditional chemical fingerprints (Morgan/ECFP, MACCS Keys) with deep learning embeddings (ChemBERTa).
*   **Chemical Property Analysis**: Calculates Molecular Weight, LogP, HBD, HBA, PSA, Rotatable Bonds, QED, Fsp3, and checks bRo5 compliance.
*   **High-Performance**: Supports GPU acceleration and batch processing for rapid analysis of large libraries.
*   **Model Flexibility**: Built-in support for multiple pre-trained models (ChemBERTa variants) and custom fine-tuned models.
*   **Fine-Tuning Module**: Train models on your specific biological targets (e.g., MDR, specific receptors) using your own data.
*   **Interactive Reporting**: Generates detailed HTML reports with interactive Plotly charts (Chemical Space, Compliance Status).
*   **DeepPurpose Integration**: Optional integration with DeepPurpose for DTI/ADMET screening.

## Installation

1.  **Clone the repository**:
    ```bash
    git clone <repository-url>
    cd Ahmed
    ```

2.  **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### 1. Basic Prediction & Analysis
Run the main analysis script to compare your compounds against a set of controls.

```bash
python main.py --compounds "path/to/compounds.csv" --controls "path/to/controls.csv"
```

**Arguments:**
*   `--compounds`: CSV file containing compounds to analyze (must have a `smiles` column).
*   `--controls`: CSV file containing control compounds (must have a `smiles` column).
*   `--model`: (Optional) HuggingFace model name or alias. Default: `chemberta-base`.
*   `--device`: (Optional) `cpu`, `cuda` (NVIDIA GPU), or `mps` (Mac Metal). Default: `cpu`.
*   `--output`: (Optional) Prefix for output files. Default: `results`.
*   `--use_deeppurpose`: (Optional) Flag to enable DeepPurpose screening.

### 2. Using GPU and Specific Models
To speed up processing using a GPU and select a specific model variant:

```bash
python main.py --compounds "data.csv" --controls "ctrl.csv" --device cuda --model chemberta-77m
```

**Supported Model Aliases:**
*   `chemberta-base`: `seyonec/ChemBERTa-zinc-base-v1`
*   `chemberta-77m`: `seyonec/ChemBERTa-zinc-deepchem-77m`
*   `chemberta-mtr`: `DeepChem/ChemBERTa-77M-MTR`
*   `chemberta-mlm`: `DeepChem/ChemBERTa-77M-MLM`

### 3. Fine-Tuning for Specific Targets
Train the LLM on your own dataset (e.g., for MDR activity) to create a specialized model.

**Data Format**: Your training CSV must have `smiles` and `label` (0 or 1) columns.

```bash
python train.py --train_data "mdr_data.csv" --output_dir "mdr_model" --epochs 5 --device cuda
```

**Using the Fine-Tuned Model**:
After training, use the saved model directory in the main script:

```bash
python main.py --compounds "compounds.csv" --controls "controls.csv" --model "./mdr_model"
```

## Project Structure

```
Ahmed/
├── main.py                     # Main CLI for prediction and analysis
├── train.py                    # CLI for fine-tuning models
├── requirements.txt            # Python dependencies
├── README.md                   # Project documentation
└── src/
    ├── descriptors/
    │   ├── chemical.py         # RDKit descriptors (ECFP, MACCS, Properties)
    │   └── llm.py              # LLM embedding generation (Batch & GPU support)
    ├── screening/
    │   └── deeppurpose_module.py # DeepPurpose integration
    ├── finetuning/
    │   ├── trainer.py          # Training loop logic
    │   └── dataset.py          # Data preparation
    └── utils/
        ├── data_loader.py      # CSV handling
        └── report.py           # HTML/Plotly report generation
```

## Output

The tool generates two files:
1.  **`results.csv`**: A raw data file containing all calculated properties, similarity scores, and embeddings.
2.  **`results.html`**: An interactive HTML report visualizing the chemical space and analysis summary.
