#!/usr/bin/env python3
"""
- Purpose: Link consecutive modules (n -> n+1) into reaction pairs for each BGC. 
- Actions:
    * Input: iPKS_rxn.cleaned.csv  (from 00_validate_and_normalize)
    * Output: iPKS_rxn.reactions.csv
    * Each row = one reaction step: reactant_smiles, product_smiles, module_from, module_to, domains
    * Canonicalizes column names using aliases (bgc_id, step_idx, domains, etc.) 
    * Deduplicates rows and writes: iPKS_rxn.reactions.csv
"""

import pandas as pd
from pathlib import Path

IN_CSV  = Path("../data/iPKS_rxn.cleaned.csv")
OUT_CSV = Path("../data/iPKS_rxn.reactions.csv")

def main():
    df = pd.read_csv(IN_CSV, dtype=str)
    # column aliases
    col_id  = "bgc_id" if "bgc_id" in df.columns else "BGC"
    col_mod = "step_idx" if "step_idx" in df.columns else "module_num"
    col_mid = "intermediate_smiles" if "intermediate_smiles" in df.columns else "smiles"
    
    # Sort
    df[col_mod] = df[col_mod].astype(int)
    df = df.sort_values([col_id, col_mod])
    
    records = []
    for bgc, sub in df.groupby(col_id):
        sub = sub.sort_values(col_mod).reset_index(drop=True)
        for i in range(len(sub)-1):
            r = sub.loc[i, col_mid]
            p = sub.loc[i+1, col_mid]
            if pd.isna(r) or pd.isna(p): 
                continue
            rec = {
                "bgc_id": bgc,
                "module_from": int(sub.loc[i, col_mod]),
                "module_to": int(sub.loc[i+1, col_mod]),
                "reactant_smiles": r,
                "product_smiles": p,
                "domains": sub.loc[i+1, "domains"] if "domains" in sub.columns else "",
                "at_substrate": sub.loc[i+1, "at_annotation"] if "at_annotation" in sub.columns else "",
                "kr_annotation": sub.loc[i+1, "kr_annotation"] if "kr_annotation" in sub.columns else "",
                "dh_annotation": sub.loc[i+1, "dh_annotation"] if "dh_annotation" in sub.columns else "",
                "er_annotation": sub.loc[i+1, "er_annotation"] if "er_annotation" in sub.columns else "",
                "te_annotation": sub.loc[i+1, "te_annotation"] if "te_annotation" in sub.columns else "",
            }
            records.append(rec)

    out = pd.DataFrame(records)
    out.to_csv(OUT_CSV, index=False)
    print(f"[OK] Wrote {OUT_CSV} ({len(out)} reactions)")

if __name__ == "__main__":
    main()
