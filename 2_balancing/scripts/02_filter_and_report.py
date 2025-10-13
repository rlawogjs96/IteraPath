#!/usr/bin/env python3
# 2_balancing/scripts/02_filter_and_report.py
# - inputs:  ../data/balance_results.csv, ../data/syntemp_input.standardized.csv
#            ../../1_preprocessing/data/syntemp_input_meta.csv
# - outputs: ../data/balanced_only.csv, ../data/syntemp_ready.csv, ../data/filter_summary.md

from pathlib import Path
import pandas as pd

# ---- resolve paths relative to this script ----
BASE_DIR = Path(__file__).resolve().parent            # .../2_balancing/scripts
DATA_DIR = BASE_DIR.parent / "data"                   # .../2_balancing/data
PREP_DATA_DIR = BASE_DIR.parent.parent / "1_preprocessing" / "data"  # .../1_preprocessing/data

BALCSV   = DATA_DIR / "balance_results.csv"
STDCSV   = DATA_DIR / "syntemp_input.standardized.csv"
META     = PREP_DATA_DIR / "syntemp_input_meta.csv"

OUTOK    = DATA_DIR / "balanced_only.csv"
OUTREADY = DATA_DIR / "syntemp_ready.csv"
OUTMD    = DATA_DIR / "filter_summary.md"

def main():
    if not BALCSV.exists() or not STDCSV.exists():
        raise SystemExit("[ERROR] Run 00_standardize_with_syntemp.py and 01_balance_with_syntemp.py first.")

    bal = pd.read_csv(BALCSV, dtype=str, keep_default_na=False)
    std = pd.read_csv(STDCSV, dtype=str, keep_default_na=False)

    ok_ids = set(bal.loc[bal["balanced"].astype(bool), "reaction_id"])

    ok_std = std[std["reaction_id"].isin(ok_ids)].copy()

    if META.exists():
        meta = pd.read_csv(META, dtype=str, keep_default_na=False)
        ok_std = ok_std.merge(meta, on="reaction_id", how="left")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ok_std.to_csv(OUTOK, index=False)

    ready = (
        ok_std[["reaction_id", "reactant_smiles_std", "product_smiles_std"]]
        .rename(columns={
            "reactant_smiles_std": "reactant_smiles",
            "product_smiles_std":  "product_smiles"
        })
        .drop_duplicates()
    )
    ready.to_csv(OUTREADY, index=False)

    by_step = ok_std["step_type"].value_counts(dropna=False).to_dict() if "step_type" in ok_std.columns else {}
    OUTMD.write_text(
        "\n".join([
            "# Filter Summary",
            f"- Balanced rows exported: {len(ok_std)}",
            f"- By step_type: {by_step}",
            "",
            "Files:",
            f"- {OUTOK}",
            f"- {OUTREADY}",
        ]),
        encoding="utf-8"
    )

    print(f"[OK] balanced_only: {len(ok_std)} rows -> {OUTOK}")
    print(f"[OK] syntemp_ready : {len(ready)} rows -> {OUTREADY}")
    print(f"[OK] wrote {OUTMD}")

if __name__ == "__main__":
    main()
