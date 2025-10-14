# 5_planner/01_mcs_template_selector.py
from pathlib import Path
import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, DataStructs

# Input/output paths (relative to 5_planner root)
BANK   = Path("./data/template_bank.csv")
TARGET = Path("../cases/orthosporin/target.smi")
OUTJS  = Path("./data/selector_out.json")

def load_target_smiles(p: Path) -> str:
    """Load target SMILES (first line)"""
    text = p.read_text(encoding="utf-8").strip()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        raise ValueError(f"Empty target file: {p}")
    return lines[0]

def morgan_fp(mol, radius=2, nbits=2048):
    """Generate Morgan fingerprint"""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)

def score_pair(target_mol, core_mol):
    """
    Tanimoto similarity + MCS atoms 계산
    
    Returns:
        (tanimoto, mcs_atoms)
    """
    # Tanimoto similarity
    fp_t = morgan_fp(target_mol)
    fp_c = morgan_fp(core_mol)
    tanim = DataStructs.TanimotoSimilarity(fp_t, fp_c)
    
    # MCS (Maximum Common Substructure)
    try:
        mcs = rdFMCS.FindMCS([target_mol, core_mol], timeout=5)
        mcs_atoms = mcs.numAtoms if mcs else 0
    except Exception:
        mcs_atoms = 0
    
    return tanim, mcs_atoms

def main():
    print("[START] MCS template selector...")
    
    # Load inputs
    bank = pd.read_csv(BANK)
    print(f"[INFO] Loaded template bank: {len(bank)} BGCs")
    
    target_smi = load_target_smiles(TARGET)
    print(f"[INFO] Target SMILES: {target_smi}")
    
    tmol = Chem.MolFromSmiles(target_smi)
    if tmol is None:
        raise ValueError(f"Invalid target SMILES: {target_smi}")
    
    # Compare each BGC core with target
    results = []
    for _, row in bank.iterrows():
        cmol = Chem.MolFromSmiles(row["core_smiles"])
        if cmol is None:
            print(f"[WARN] Invalid core SMILES for {row['bgc_id']}: {row['core_smiles']}")
            continue
        
        tanim, mcs_atoms = score_pair(tmol, cmol)
        results.append({
            "bgc_id": str(row["bgc_id"]),
            "tanimoto": round(float(tanim), 6),
            "mcs_atoms": int(mcs_atoms),
            "core_smiles": row["core_smiles"],
            "rules_seq": row["rules_seq"],
            "n_ext": int(row["n_ext"]),
            "has_kr": bool(row["has_kr"]),
            "has_dh": bool(row["has_dh"]),
            "closure_type": str(row.get("closure_type", "")),
            "at_substrate_set": str(row.get("at_substrate_set", "")),
        })
    
    if not results:
        raise RuntimeError("No valid BGC cores found in template bank")
    
    # Deterministic sorting: n_ext>=5 highest priority, then (-tanimoto, -mcs_atoms, bgc_id)
    # Prioritize n_ext>=5 candidates for hexaketide targets like Orthosporin
    required_n_ext = 5
    results.sort(key=lambda x: (
        0 if x["n_ext"] >= required_n_ext else 1,  # n_ext>=5 first
        -x["tanimoto"],
        -x["mcs_atoms"],
        x["bgc_id"]
    ))
    
    # Output
    OUTJS.parent.mkdir(parents=True, exist_ok=True)
    output = {
        "target": target_smi,
        "ranked": results,
        "top1": results[0] if results else None
    }
    
    with open(OUTJS, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"[DONE] Selector out written: {OUTJS}")
    if results:
        top = results[0]
        print(f"       Top-1: {top['bgc_id']} | Tanimoto={top['tanimoto']:.4f} | MCS_atoms={top['mcs_atoms']}")

if __name__ == "__main__":
    main()
