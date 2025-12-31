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
*   `--dp_target`: (Required if using DeepPurpose) Amino acid sequence of the target protein.
*   `--dp_model`: (Optional) DeepPurpose pre-trained model name. Default: `MPNN_CNN_BindingDB`.

### 2. DeepPurpose Screening (DTI)
To screen your compounds against a specific biological target (e.g., MDR1/P-gp) using Deep Learning:

```bash
python main.py \
  --compounds "compounds.csv" \
  --controls "controls.csv" \
  --use_deeppurpose \
  --dp_target "MVS..." \
  --dp_model "MPNN_CNN_BindingDB"
```

**Available DeepPurpose Models:**

| Model Name | Drug Encoder | Target Encoder | Dataset | Best For |
| :--- | :--- | :--- | :--- | :--- |
| `MPNN_CNN_BindingDB` | Graph (MPNN) | CNN | BindingDB | **General Purpose (Recommended)** |
| `CNN_CNN_BindingDB` | CNN | CNN | BindingDB | High Speed |
| `Morgan_AAC_BindingDB` | Fingerprint | Amino Acid Comp. | BindingDB | Baseline / Classical |
| `MPNN_CNN_DAVIS` | Graph (MPNN) | CNN | DAVIS | **Kinase Targets** |
| `MPNN_CNN_KIBA` | Graph (MPNN) | CNN | KIBA | Kinase Targets |

### 3. Using GPU and Specific LLM Models
To speed up processing using a GPU and select a specific LLM variant:

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
Screening_Docking/
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

## Input Data Format

The tool accepts CSV files for compounds and controls. The CSV files must contain specific columns for the tool to work correctly. Column names are case-insensitive.

### 1. Compounds CSV (`--compounds`)
*   **Required Column**: `smiles` (The SMILES string of the compound).
*   **Optional Columns**: Any other columns will be preserved but not used for calculation.

**Example:**
| smiles | ID |
| :--- | :--- |
| CCO | Comp_1 |
| c1ccccc1 | Comp_2 |

### 2. Controls CSV (`--controls`)
*   **Required Column**: `smiles` (The SMILES string of the control compound).
*   **Optional Column**: `nama` or `nama_kontrol` (Name of the control). If missing, controls will be named `Ctrl_1`, `Ctrl_2`, etc.

**Example:**
| smiles | nama_kontrol |
| :--- | :--- |
| CC(=O)OC1=CC=CC=C1C(=O)O | Aspirin |
| CN1C=NC2=C1C(=O)N(C(=O)N2C)C | Caffeine |

### 3. Training Data CSV (for Fine-Tuning)
*   **Required Columns**:
    *   `smiles`: The chemical structure.
    *   `label`: Binary label (`0` for Inactive, `1` for Active).

## Output

The tool generates two files:
1.  **`results.csv`**: A raw data file containing all calculated properties, similarity scores, and embeddings.
2.  **`results.html`**: An interactive HTML report visualizing the chemical space and analysis summary.
