#!/usr/bin/env python3
"""
Export minimal SynTemp-ready input:
- Keep rows with BOTH reactant_smiles and product_smiles present
- Include reaction_id, reactant_smiles, product_smiles
- Also write a sidecar metadata table for later (domains, AT, closure, etc.)
Input : 1_preprocessing/iPKS_rxn.reactions.labeled.csv
Output: 1_preprocessing/syntemp_input.csv
        1_preprocessing/syntemp_input_meta.csv
"""

from pathlib import Path
import pandas as pd

INCSV   = Path("../data/iPKS_rxn.reactions.labeled.csv")
OUTCSV  = Path("../data/syntemp_input.csv")
OUTMETA = Path("../data/syntemp_input_meta.csv")

def main():
    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)

    # basic filter: keep only rows that have both sides
    mask = df["reactant_smiles"].astype(str).str.len().gt(0) & df["product_smiles"].astype(str).str.len().gt(0)
    keep = df[mask].copy()

    # recommended: keep only known step_types (EXT/KR/DH/CLOSURE). Exclude OTHER
    """
        - Keep only known step_types (EXT/KR/DH/CLOSURE). 
        - Exclude OTHER
    """
    allowed = {"EXT","KR","DH","ER","CLOSURE"}
    keep = keep[keep["step_type"].isin(allowed)].copy()

    # sanitize metadata values like literal 'nan'
    if "at_substrate" in keep.columns:
        keep["at_substrate"] = keep["at_substrate"].replace({"nan":"", "NaN":"", "NONE":"", "None":""})

    # minimal SynTemp input
    syn = keep[["reaction_id","reactant_smiles","product_smiles"]].drop_duplicates().reset_index(drop=True)
    syn.to_csv(OUTCSV, index=False)
    print(f"[OK] wrote {OUTCSV} (rows={len(syn)})")

    # sidecar metadata (for mapping/rule tags later)
    meta_cols = [c for c in [
        "reaction_id","bgc_id","module_from","module_to",
        "domains","domains_norm","domains_norm_no_acp",
        "at_substrate","kr_annotation","dh_annotation","er_annotation","te_annotation",
        "step_type","closure"
    ] if c in keep.columns]
    meta = keep[meta_cols].drop_duplicates().reset_index(drop=True)
    meta.to_csv(OUTMETA, index=False)
    print(f"[OK] wrote {OUTMETA} (rows={len(meta)})")

if __name__ == "__main__":
    main()
