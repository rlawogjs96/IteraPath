#!/usr/bin/env python3
"""
syntemp_input.csv --> minimal, unbalanced reaction set exported by pre-processing. 
Applying SynTemp's Deionize and Tautomerize ensures consistent SMILES before stoichiometry augmentation and balance checks. 
SynTemp emphasizes that reaction balancing is a pre-requisite because it ensures that the AAM and the ITS are chemically consistent. 

Unbalanced reactions—those missing atoms (especially hydrogens, charges, or counterions)—lead to incorrect or incomplete AAMs.
These, in turn, distort the reaction center (the subgraph representing changing bonds) and break the morphisms needed for Double Pushout (DPO) rewriting.

The authors explicitly note that "balanced reactions are required so that the mapping between reactant and product atoms can be one-to-one," which guarantees mass conservation and allows SynTemp to construct correct ITS graphs.

Because DPO rewriting treats reactions as formal graph transformations, an unbalanced equation would correspond to a non-invertible or undefined morphism — effectively, a broken category-theoretic map.

SynTemp's hydrogen-completion routines (section 2.4 of the paper) further explain that missing hydrogens make it impossible to determine whether a reaction conserves charge and mass, so they augment incomplete AAMs to include explicit hydrogen counts before rule extraction
"""
from pathlib import Path
import pandas as pd
import warnings

# Suppress FutureWarnings from SynTemp
warnings.filterwarnings("ignore", category=FutureWarning, module="syntemp")

# SynTemp standard chemistry utils
from syntemp.SynChemistry.deionize import Deionize
from syntemp.SynChemistry.tautomerize import Tautomerize

INCSV  = Path("../../1_preprocessing/data/syntemp_input.csv")  # reaction_id, reactant_smiles, product_smiles
#INCSV  = Path("../../1_preprocessing/data/syntemp_input_with_mal.csv")
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

    # Keep only reaction_id and standardized SMILES columns
    df = df[["reaction_id", "reactant_smiles_std", "product_smiles_std"]]

    # Remove duplicates (same reaction_id, same standardized reactant/product)
    df = df.drop_duplicates(subset=["reaction_id", "reactant_smiles_std", "product_smiles_std"])

    df.to_csv(OUTCSV, index=False)
    print(f"[OK] wrote {OUTCSV} (rows={len(df)})")

if __name__ == "__main__":
    main()
