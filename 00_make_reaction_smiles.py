#!/usr/bin/env python3
# 3_templates/scripts/00_make_reaction_smiles.py
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent
INCSV  = BASE.parent.parent / "2_balancing" / "data" / "syntemp_ready.csv"
OUTDIR = BASE.parent / "data"
OUTCSV = OUTDIR / "reactions.csv"

def main():
    if not INCSV.exists():
        raise SystemExit(f"Missing {INCSV}")
    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)
    need = {"reaction_id","reactant_smiles","product_smiles"}
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"Missing columns: {miss}")
    df["reactions"] = df["reactant_smiles"] + ">>" + df["product_smiles"]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df[["reaction_id","reactions"]].to_csv(OUTCSV, index=False)
    print(f"[OK] wrote {OUTCSV} rows={len(df)}")

if __name__ == "__main__":
    main()
