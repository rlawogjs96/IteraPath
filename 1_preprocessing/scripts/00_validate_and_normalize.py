#!/usr/bin/env python3
"""
Validate & normalize iPKS_rxn.csv
- Locks a canonical schema via alias mapping
- Trims/normalizes values
- Tokenizes domain list -> KS/AT/KR/DH/ER/TE/PT
- Optionally validates/canonicalizes SMILES (if RDKit available)
Outputs:
  - iPKS_rxn.cleaned.csv
  - qc_report.md
"""

import sys, re, csv, os
from pathlib import Path
import pandas as pd

# ---- Config ----
INFILE  = Path("../data/iPKS_rxn.csv")
OUT_CSV = Path("../data/iPKS_rxn.cleaned.csv")
OUT_QC  = Path("../data/qc_report.md")

# Column aliases -> canonical names
# (필요시 여기에 별칭을 계속 추가해도 됨)
ALIAS = {
    "bgc": "bgc_id", "bgc_id": "bgc_id", "bgc accession": "bgc_id",

    "reaction_id": "reaction_id", "rxn_id": "reaction_id",
    "step_idx": "step_idx",                 # 이미 step_idx로 들어온 경우 유지
    "module_num": "step_idx",               # ★ 숫자 인덱스는 step_idx로
    "module": "module_name",                # ★ 모듈 '이름'은 module_name으로

    "reactant": "reactant_smiles", "reactant_smiles": "reactant_smiles", "substrate_smiles": "reactant_smiles",
    "product": "product_smiles", "product_smiles": "product_smiles",

    "domains": "domains", "domain_list": "domains", "domain": "domains",

    "at": "at_substrate", "at_substrate": "at_substrate", "extender": "at_substrate",

    "closure": "closure", "pt": "closure", "ring_closure": "closure",

    "notes": "notes", "comment": "notes"
}


# Canonical column order (없어도 진행되지만, 가능하면 맞추자)
CANON = [
    "bgc_id", "reaction_id", "step_idx", "module_name",
    "reactant_smiles", "product_smiles",
    "domains", "at_substrate", "closure",
    "notes"
]

# Allowed domain tokens
DOMAIN_VOCAB = {"KS","AT","KR","DH","ER","TE","PT"}

# Try RDKit (optional)
try:
    from rdkit import Chem
    from rdkit.Chem import MolToSmiles
    RDKit_OK = True
except Exception:
    RDKit_OK = False


def _strip_norm(x: str):
    if x is None: return ""
    x = str(x).strip()
    if x.upper() in {"NA","N/A","NONE","NULL"}: return ""
    return x

def normalize_domains(cell: str):
    """Split on comma/semicolon/space, normalize tokens -> KS/AT/KR/DH/ER/TE/PT"""
    if not cell: return ""
    # split on common separators
    toks = re.split(r"[;,/ ]+", cell.strip())
    toks = [t for t in toks if t]
    norm = []
    for t in toks:
        u = t.upper().strip()
        # quick aliasing
        if u in {"ACP"}: continue  # carrier는 제외
        if u in {"THIOESTERASE"}: u = "TE"
        if u in {"PRODUCTTEMPLATE","PTD","PT_DOMAIN"}: u = "PT"
        if u in {"KETOREDUCTASE"}: u = "KR"
        if u in {"DEHYDRATASE"}: u = "DH"
        if u in {"ENRED","ENOYLREDUCTASE"}: u = "ER"
        if u in {"BETA-KETOACYLSYNTHASE","KSQ"}: u = "KS"
        if u in {"ACYLTRANSFERASE"}: u = "AT"
        if u in DOMAIN_VOCAB:
            norm.append(u)
    # dedup, keep stable order
    seen, out = set(), []
    for u in norm:
        if u not in seen:
            out.append(u); seen.add(u)
    return ",".join(out)

def canon_smiles(smi: str):
    if not smi or not RDKit_OK:
        return smi, False
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: return smi, True
        can = MolToSmiles(mol, canonical=True)
        return can, False
    except Exception:
        return smi, True

def main():
    if not INFILE.exists():
        print(f"[ERROR] Missing {INFILE}. Put iPKS_rxn.csv there.")
        sys.exit(1)

    df = pd.read_csv(INFILE, dtype=str, keep_default_na=False)
    # 1) Column alias -> canonical
    new_cols = {}
    for c in df.columns:
        key = c.strip()
        low = key.lower()
        new_cols[c] = ALIAS.get(low, c.strip().lower())  # fallback: lower
    df = df.rename(columns=new_cols)

    # 2) Trim and normalize values
    for c in df.columns:
        df[c] = df[c].map(_strip_norm)

    # 3) Try to ensure canonical columns exist (create empty if missing)
    for c in CANON:
        if c not in df.columns:
            df[c] = ""

    # 4) Domain normalization
    if "domains" in df.columns:
        df["domains"] = df["domains"].map(normalize_domains)

    # 5) RDKit-based SMILES canon (optional)
    r_errs, p_errs = 0, 0
    if RDKit_OK:
        if "reactant_smiles" in df.columns:
            cans, errs = [], 0
            for s in df["reactant_smiles"].tolist():
                can, bad = canon_smiles(s)
                cans.append(can)
                errs += int(bad)
            df["reactant_smiles"] = cans
            r_errs = errs
        if "product_smiles" in df.columns:
            cans, errs = [], 0
            for s in df["product_smiles"].tolist():
                can, bad = canon_smiles(s)
                cans.append(can)
                errs += int(bad)
            df["product_smiles"] = cans
            p_errs = errs

    # 6) Drop exact duplicate rows
    before = len(df)
    df = df.drop_duplicates()
    dropped = before - len(df)

    # 7) Reorder columns (canonical first)
    ordered = [c for c in CANON if c in df.columns] + [c for c in df.columns if c not in CANON]
    df = df[ordered]

    # 8) Basic QC stats
    lines = []
    lines.append("# iPKS_rxn QC Report\n")
    lines.append(f"- Input rows: {before}")
    lines.append(f"- Output rows (deduped): {len(df)}  (dropped: {dropped})")
    # null rates on key cols
    key_cols = [c for c in ["bgc_id","reaction_id","step_idx","reactant_smiles","product_smiles","domains","at_substrate","closure"] if c in df.columns]
    if key_cols:
        lines.append("\n## Null rates (key columns)\n")
        nulls = (df[key_cols] == "").mean().sort_values(ascending=False)
        for c,v in nulls.items():
            lines.append(f"- {c}: {v:.2%}")

    # domain vocab check
    if "domains" in df.columns:
        lines.append("\n## Domain token distribution\n")
        from collections import Counter
        tok = []
        for cell in df["domains"].tolist():
            if not cell: continue
            tok += cell.split(",")
        cnt = Counter(tok)
        for k,v in cnt.most_common():
            lines.append(f"- {k}: {v}")

    # SMILES errors
    if RDKit_OK:
        lines.append("\n## SMILES canonicalization")
        lines.append(f"- reactant parse errors: {r_errs}")
        lines.append(f"- product  parse errors: {p_errs}")
    else:
        lines.append("\n## SMILES canonicalization")
        lines.append("- RDKit not available; skipped.")

    # 9) Save
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False, quoting=csv.QUOTE_MINIMAL)
    OUT_QC.write_text("\n".join(lines), encoding="utf-8")

    print(f"[OK] Wrote {OUT_CSV}")
    print(f"[OK] Wrote {OUT_QC}")

if __name__ == "__main__":
    main()
