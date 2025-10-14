#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import warnings

# Suppress FutureWarnings from SynTemp
warnings.filterwarnings("ignore", category=FutureWarning, module="syntemp")

# SynTemp standard chemistry utils
from syntemp.SynChemistry.deionize import Deionize
from syntemp.SynChemistry.tautomerize import Tautomerize

INCSV  = Path("../../1_preprocessing/data/syntemp_input.csv")  # reaction_id, reactant_smiles, product_smiles
OUTCSV = Path("../data/syntemp_input.standardized.csv")
OUTCSV.parent.mkdir(parents=True, exist_ok=True)

def stdz(smiles: str) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return smiles
    try:
        s = Deionize(smiles)
        s = Tautomerize(s)
        return s
    except Exception:
        return smiles

def main():
    if not INCSV.exists():
        raise SystemExit(f"[ERROR] Missing {INCSV}")

    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)

    needed = {"reaction_id", "reactant_smiles", "product_smiles"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] Required columns missing: {missing}")

    df["reactant_smiles_std"] = df["reactant_smiles"].map(stdz)
    df["product_smiles_std"]  = df["product_smiles"].map(stdz)

    # Remove duplicates (same reaction_id, same standardized reactant/product)
    df = df.drop_duplicates(subset=["reaction_id", "reactant_smiles_std", "product_smiles_std"])

    df.to_csv(OUTCSV, index=False)
    print(f"[OK] wrote {OUTCSV} (rows={len(df)})")

if __name__ == "__main__":
    main()
