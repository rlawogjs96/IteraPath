#!/usr/bin/env python3
# 5_planner/scripts/00_build_template_bank.py
from pathlib import Path
import pandas as pd

# Resolve paths relative to this script location (repo-root safe)
BASE     = Path(__file__).resolve().parent
BALANCED = BASE.parent.parent / "2_balancing" / "data" / "balanced_only.csv"
RULEMETA = BASE.parent.parent / "4_rulemeta" / "data" / "rulemeta.csv"
OUTBANK  = BASE.parent / "data" / "template_bank.csv"

def load_tables():
    """Load input files"""
    bal = pd.read_csv(BALANCED)
    meta = pd.read_csv(RULEMETA)
    # Convert reaction_id to string
    for df in (bal, meta):
        if "reaction_id" in df.columns:
            df["reaction_id"] = df["reaction_id"].astype(str)
    return bal, meta

def build_bank(bal: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Build template bank per BGC
    
    For each BGC:
    - Final core SMILES (product of CLOSURE step)
    - Rule sequence (reaction_id list)
    - Meta summary (n_ext, has_kr, has_dh, closure_type, at_substrate_set)
    """
    # Start with balanced_only as base
    ready = bal.copy()
    
    # Merge step_type, closure info from rulemeta
    if not meta.empty and "reaction_id" in meta.columns:
        meta_cols = ["reaction_id", "step_type", "closure"]
        available_cols = [c for c in meta_cols if c in meta.columns]
        ready = ready.merge(
            meta[available_cols].drop_duplicates("reaction_id"),
            on="reaction_id",
            how="left",
            suffixes=("_bal", "_meta")
        )
        
        # Unify step_type: _meta takes priority (_bal may not exist in balanced_only)
        if "step_type_meta" in ready.columns:
            ready["step_type"] = ready["step_type_meta"]
        elif "step_type_bal" in ready.columns:
            ready["step_type"] = ready["step_type_bal"]
    
    # Infer if step_type missing
    if "step_type" not in ready.columns:
        ready["step_type"] = "UNKNOWN"
    
    # Extract final core function
    def final_core(dfbgc: pd.DataFrame) -> str:
        # CLOSURE step's product_smiles takes priority
        if "step_type" in dfbgc.columns:
            clo = dfbgc[dfbgc["step_type"] == "CLOSURE"]
            if not clo.empty:
                sorted_clo = clo.sort_values(["module_from", "module_to"], na_position='last')
                # Prefer standardized column if present
                if "product_smiles_std" in sorted_clo.columns:
                    return sorted_clo.iloc[-1]["product_smiles_std"]
                if "product_smiles" in sorted_clo.columns:
                    return sorted_clo.iloc[-1]["product_smiles"]
        # fallback: last product_smiles in BGC
        if "module_from" in dfbgc.columns and "module_to" in dfbgc.columns:
            sorted_bgc = dfbgc.sort_values(["module_from", "module_to"], na_position='last')
        else:
            sorted_bgc = dfbgc
        if "product_smiles_std" in sorted_bgc.columns:
            return sorted_bgc.iloc[-1]["product_smiles_std"]
        return sorted_bgc.iloc[-1].get("product_smiles", "")
    
    rows = []
    for bgc, dfb in ready.groupby("bgc_id"):
        dfb = dfb.sort_values(["module_from", "module_to"], na_position='last').reset_index(drop=True)
        
        # Final core SMILES
        core = final_core(dfb)
        
        # Rule sequence
        rules_seq = dfb["reaction_id"].astype(str).tolist()
        
        # AT substrate set
        at_set = []
        if "at_substrate" in dfb.columns:
            at_set = sorted(set(dfb["at_substrate"].dropna().astype(str).tolist()))
            at_set = [x for x in at_set if x and x.lower() not in ['nan', 'none', '']]
        
        # EXT count (enhanced debugging)
        n_ext = 0
        if "step_type" in dfb.columns:
            ext_mask = (dfb["step_type"] == "EXT")
            n_ext = int(ext_mask.sum())
            # Debugging: fallback if all step_type are NaN/UNKNOWN
            if n_ext == 0:
                # Estimate: rules_seq length minus 1 CLOSURE
                total_steps = len(dfb)
                closure_count = int((dfb["step_type"] == "CLOSURE").sum())
                if closure_count > 0:
                    n_ext = total_steps - closure_count
                else:
                    # Last resort fallback: total steps - 1
                    n_ext = max(0, total_steps - 1)
        
        # Check KR/DH presence
        has_kr = False
        has_dh = False
        if "domains" in dfb.columns:
            has_kr = bool(dfb["domains"].fillna("").astype(str).str.contains("KR", case=False).any())
            has_dh = bool(dfb["domains"].fillna("").astype(str).str.contains("DH", case=False).any())
        
        # Closure type
        closure_type = "PT_OR_TE"
        if "closure" in dfb.columns:
            closure_vals = dfb["closure"].dropna().astype(str).unique().tolist()
            closure_vals = [x for x in closure_vals if x and x.lower() not in ['nan', 'none', '']]
            if closure_vals:
                closure_type = ";".join(sorted(closure_vals))
        
        rows.append({
            "bgc_id": bgc,
            "core_smiles": core,
            "rules_seq": "|".join(rules_seq),
            "n_ext": n_ext,
            "has_kr": has_kr,
            "has_dh": has_dh,
            "closure_type": closure_type,
            "at_substrate_set": "|".join(at_set) if at_set else "",
        })
    
    return pd.DataFrame(rows)

def main():
    print("[START] Building template bank...")
    OUTBANK.parent.mkdir(parents=True, exist_ok=True)
    
    bal, meta = load_tables()
    print(f"[INFO] Loaded: balanced_only={len(bal)} rows, rulemeta={len(meta)} rows")
    
    bank = build_bank(bal, meta)
    bank.to_csv(OUTBANK, index=False)
    
    print(f"[DONE] Template bank written: {OUTBANK}")
    print(f"       Rows: {len(bank)}")
    print(f"       BGCs: {', '.join(bank['bgc_id'].astype(str).tolist())}")

if __name__ == "__main__":
    main()
