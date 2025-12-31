import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, QED, MACCSkeys

def get_mol(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except:
        return None

def calculate_ecfp(mol, radius=2, nBits=2048):
    """Calculates Morgan Fingerprint (ECFP)"""
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)

def calculate_maccs(mol):
    """Calculates MACCS Keys Fingerprint"""
    if mol is None:
        return None
    return MACCSkeys.GenMACCSKeys(mol)

def calculate_properties(mol):
    """Calculates comprehensive physicochemical properties including bRo5 and QED"""
    if mol is None:
        return None
    
    # Basic Properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    psa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    
    # Advanced Properties
    qed = QED.qed(mol)
    csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    
    # bRo5 conditions
    # MW <= 1000, -2 <= LogP <= 10, HBD <= 6, HBA <= 15, PSA <= 250, RotB <= 20
    cond = [mw<=1000, -2<=logp<=10, hbd<=6, hba<=15, psa<=250, rotb<=20]
    violations = 6 - sum(cond)
    status = "PASS" if violations == 0 else f"FAIL({violations})"
    
    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': hbd,
        'HBA': hba,
        'PSA': round(psa, 2),
        'RotB': rotb,
        'QED': round(qed, 3),
        'Fsp3': round(csp3, 3),
        'NumRings': rings,
        'NumAromatic': aromatic_rings,
        'bRo5_Status': status
    }

# Alias for backward compatibility
calculate_bro5 = calculate_properties

def calculate_tanimoto(fp1, fp2):
    """Calculates Tanimoto similarity between two fingerprints"""
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def calculate_bulk_tanimoto(query_fp, list_fps):
    """
    Calculates Tanimoto similarity between a query fingerprint and a list of fingerprints.
    Much faster than looping in Python.
    """
    if query_fp is None or not list_fps:
        return []
    return DataStructs.BulkTanimotoSimilarity(query_fp, list_fps)

def calculate_dice(fp1, fp2):
    """Calculates Dice Similarity (gives more weight to matches)"""
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.DiceSimilarity(fp1, fp2)

def calculate_tversky(fp1, fp2, alpha=0.8, beta=0.2):
    """
    Calculates Tversky Similarity (Asymmetric).
    Good for substructure matching.
    alpha=1, beta=0 -> Asymmetric containment (is fp1 inside fp2?)
    """
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TverskySimilarity(fp1, fp2, alpha, beta)
