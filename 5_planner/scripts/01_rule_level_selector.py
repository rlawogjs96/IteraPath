#!/usr/bin/env python3
# 5_planner/01_rule_level_selector.py
from pathlib import Path
import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, DataStructs

# Resolve paths relative to this script
BASE      = Path(__file__).resolve().parent
BALANCED  = BASE.parent.parent / "2_balancing" / "data" / "balanced_only.csv"
RULEMETA  = BASE.parent.parent / "4_rulemeta" / "data" / "rulemeta.csv"
BANK      = BASE.parent / "data" / "template_bank.csv"
TARGET_SMI= BASE.parent.parent / "cases" / "orthosporin" / "target.smi"
OUTJS     = BASE.parent / "data" / "selector_out.json"

required_n_ext = 5

def load_target() -> str:
    if TARGET_SMI.exists():
        txt = TARGET_SMI.read_text(encoding="utf-8").strip().splitlines()
        for line in txt:
            s = line.strip()
            if s:
                return s
    # fallback: Orthosporin from constraints in deterministic_planner
    return "CC(CC1=CC2=CC(=CC(=C2C(=O)O1)O)O)O"

def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)

def score_pair(target_mol, core_mol):
    fp_t = morgan_fp(target_mol)
    fp_c = morgan_fp(core_mol)
    tani = DataStructs.TanimotoSimilarity(fp_t, fp_c)
    try:
        mcs = rdFMCS.FindMCS([target_mol, core_mol], timeout=5)
        mcs_atoms = mcs.numAtoms if mcs else 0
    except Exception:
        mcs_atoms = 0
    return float(tani), int(mcs_atoms)

def main():
    print("[START] Rule-level selector (closure-driven)...")
    bal  = pd.read_csv(BALANCED)
    meta = pd.read_csv(RULEMETA)
    bank = pd.read_csv(BANK) if BANK.exists() else pd.DataFrame()
    for df in (bal, meta):
        if "reaction_id" in df.columns:
            df["reaction_id"] = df["reaction_id"].astype(str)

    target_smi = load_target()
    tmol = Chem.MolFromSmiles(target_smi)
    if tmol is None:
        raise ValueError(f"Invalid target SMILES: {target_smi}")

    # Focus on closure steps; prefer standardized product SMILES
    cols = set(bal.columns)
    prod_col = "product_smiles_std" if "product_smiles_std" in cols else "product_smiles"
    if prod_col not in cols:
        raise RuntimeError("balanced_only.csv missing product_smiles/_std column")

    # Merge meta to access step_type/closure per reaction
    merged = bal.merge(meta.drop_duplicates("reaction_id"), on="reaction_id", how="left", suffixes=("", "_m"))
    if "step_type" not in merged.columns:
        # try meta suffix
        if "step_type_m" in merged.columns:
            merged["step_type"] = merged["step_type_m"]
        else:
            merged["step_type"] = ""

    closures = merged[merged["step_type"].astype(str).str.upper().eq("CLOSURE")].copy()
    if closures.empty:
        raise RuntimeError("No CLOSURE rows found; cannot perform closure-driven selection")

    results = []
    for _, row in closures.iterrows():
        core_smi = str(row.get(prod_col, ""))
        cmol = Chem.MolFromSmiles(core_smi) if core_smi else None
        if cmol is None:
            continue
        tani, mcs_atoms = score_pair(tmol, cmol)
        bgc = str(row.get("bgc_id", ""))
        # EXT count for this BGC
        dfb = merged[merged["bgc_id"] == bgc]
        ext_n = int((dfb["step_type"].astype(str).str.upper()=="EXT").sum())
        # rules_seq for this BGC (ordered)
        sort_cols = [c for c in ["module_from","module_to"] if c in dfb.columns]
        if sort_cols:
            dfb = dfb.sort_values(sort_cols)
        rules_seq = "|".join(dfb["reaction_id"].astype(str).tolist())
        # KR/DH detection from domains
        has_kr = bool(dfb.get("domains", pd.Series(dtype=str)).fillna("").astype(str).str.contains("KR", case=False).any())
        has_dh = bool(dfb.get("domains", pd.Series(dtype=str)).fillna("").astype(str).str.contains("DH", case=False).any())
        results.append({
            "bgc_id": bgc,
            "tanimoto": round(tani,6),
            "mcs_atoms": int(mcs_atoms),
            "core_smiles": core_smi,
            "rules_seq": rules_seq,
            "n_ext": ext_n,
            "has_kr": has_kr,
            "has_dh": has_dh,
            "closure_type": str(row.get("closure","")),
            "at_substrate_set": "|".join(sorted(set(dfb.get("at_substrate", pd.Series(dtype=str)).dropna().astype(str)))) if "at_substrate" in dfb.columns else "",
        })

    if not results:
        raise RuntimeError("No closure candidates scored")

    # Score: require hexaketide first, then rank by Tanimoto, then MCS
    def sort_key(r):
        return (int(r.get("n_ext",0))==required_n_ext, r.get("tanimoto",0.0), r.get("mcs_atoms",0))
    results.sort(key=sort_key, reverse=True)

    OUTJS.parent.mkdir(parents=True, exist_ok=True)
    out = {"target": target_smi, "ranked": results, "top1": results[0]}
    with open(OUTJS, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[DONE] selector_out.json written (rule-level): {OUTJS}")
    top = results[0]
    print(f"       Top-1: {top['bgc_id']} | n_ext={top['n_ext']} | Tanimoto={top['tanimoto']} | MCS_atoms={top['mcs_atoms']}")

if __name__ == "__main__":
    main()


