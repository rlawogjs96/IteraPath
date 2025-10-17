#!/usr/bin/env python3
# 5_planner/02_rule_level_planner.py
from __future__ import annotations

from pathlib import Path
import json
from typing import List, Dict, Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


# ──────────────────────────────────────────────────────────────────────────────
# Paths (resolved relative to this script)
# ──────────────────────────────────────────────────────────────────────────────
BASE        = Path(__file__).resolve().parent
CASES_DIR   = BASE.parent.parent / "cases" / "orthosporin"
CONSTRAINTS = CASES_DIR / "constraints.json"

BALANCED    = BASE.parent.parent / "2_balancing" / "data" / "balanced_only.csv"
RULEMETA    = BASE.parent.parent / "4_rulemeta" / "data" / "rulemeta.csv"
RULES_IDX   = BASE.parent.parent / "3_templates" / "data" / "rules_index.csv"
RULES_DIR   = BASE.parent.parent / "3_templates" / "data" / "rules"

OUTDIR      = CASES_DIR
OUT_ROUTE   = OUTDIR / "route_cross.json"
OUT_AUDIT   = OUTDIR / "audit_cross.md"


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def _load_constraints() -> Dict[str, Any]:
    if CONSTRAINTS.exists():
        with open(CONSTRAINTS, encoding="utf-8") as f:
            c = json.load(f)
        return c
    # Fallback defaults (Orthosporin)
    return {
        "target": {
            "name": "Orthosporin",
            "smiles": "CC(CC1=CC2=CC(=CC(=C2C(=O)O1)O)O)O",
        },
        "core_pathway": {
            "pks_type": "I_NR",
            "starter_unit": "acetyl-CoA",
            "at_substrate": "mal|malonyl",
            "required_n_ext": 5,
            "er_allowed": False,
            "kr_mode": "late_or_once",
            "dh": "optional",
        },
        "closure": {
            "require_pt": True,
            "pt_mode_allowed": ["c2c7"],
            "release": "TE_lactonization",
        },
    }


def _load_target_smiles(constraints: Dict[str, Any]) -> str:
    tgt = constraints.get("target", {}).get("smiles")
    if isinstance(tgt, str) and tgt.strip():
        return tgt.strip()
    return "CC(CC1=CC2=CC(=CC(=C2C(=O)O1)O)O)O"


def _morgan_fp(mol, radius: int = 2, nbits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def _tanimoto(a: Chem.Mol, b: Chem.Mol) -> float:
    return float(DataStructs.TanimotoSimilarity(_morgan_fp(a), _morgan_fp(b)))


def _join_rule_index(meta: pd.DataFrame, rules_idx: pd.DataFrame) -> pd.DataFrame:
    # Expect rules_idx to provide file_stem; construct full path under RULES_DIR
    idx = rules_idx.copy()
    if "file_stem" not in idx.columns:
        # best-effort: infer from reaction_id
        if "reaction_id" in idx.columns:
            idx["file_stem"] = idx["reaction_id"].astype(str)
        else:
            idx["file_stem"] = ""
    idx["gml_file"] = idx["file_stem"].astype(str) + ".gml"
    idx["gml_path"] = idx["gml_file"].apply(lambda s: str(RULES_DIR / s))
    if "reaction_id" in idx.columns:
        idx["reaction_id"] = idx["reaction_id"].astype(str)
    if "reaction_id" in meta.columns:
        meta["reaction_id"] = meta["reaction_id"].astype(str)
    # Prefer a left join on both reaction_id and file_stem when present
    on_cols = [c for c in ["reaction_id", "file_stem"] if c in meta.columns and c in idx.columns]
    if not on_cols:
        on_cols = ["reaction_id"] if "reaction_id" in meta.columns and "reaction_id" in idx.columns else None
    if on_cols is None:
        # As a fallback, just attach gml_path by file_stem uniqueness
        merged = meta.merge(idx[["file_stem", "gml_path"]].drop_duplicates("file_stem"), on="file_stem", how="left")
    else:
        merged = meta.merge(idx[[*on_cols, "gml_path"]].drop_duplicates(on_cols), on=on_cols, how="left")
    return merged


def _sanitize_smiles(s: str) -> str:
    return (s or "").strip()


def _pick_closure(meta_joined: pd.DataFrame, target_smi: str) -> Dict[str, Any]:
    tmol = Chem.MolFromSmiles(target_smi)
    if tmol is None:
        raise ValueError(f"Invalid target SMILES: {target_smi}")
    df = meta_joined.copy()
    df["step_type"] = df.get("step_type", "").astype(str)
    closures = df[df["step_type"].str.upper().eq("CLOSURE")].copy()
    if closures.empty:
        raise RuntimeError("No CLOSURE rules available in rulemeta.csv")
    # Prefer standardized product column
    prod_col = "product_smiles_std" if "product_smiles_std" in closures.columns else (
        "product_smiles" if "product_smiles" in closures.columns else None
    )
    if prod_col is None:
        # score all as zero and pick first with gml
        closures["tanimoto"] = 0.0
    else:
        scores: List[float] = []
        for s in closures[prod_col].astype(str).tolist():
            mol = Chem.MolFromSmiles(_sanitize_smiles(s)) if s else None
            scores.append(_tanimoto(tmol, mol) if mol is not None else 0.0)
        closures["tanimoto"] = scores
    # rank by: hexaketide first (n_ext==5 if available), then tanimoto
    n_ext = closures.get("n_ext") if "n_ext" in closures.columns else None
    if n_ext is not None:
        closures["_is_hex"] = (closures["n_ext"].astype(str) == "5").astype(int)
        closures = closures.sort_values(["_is_hex", "tanimoto"], ascending=[False, False])
    else:
        closures = closures.sort_values(["tanimoto"], ascending=[False])
    row = closures.iloc[0].to_dict()
    return {
        "reaction_id": str(row.get("reaction_id", "")),
        "bgc_id": str(row.get("bgc_id", "")),
        "tanimoto": float(row.get("tanimoto", 0.0)),
        "gml_path": str(row.get("gml_path", "")),
        "product_smiles": str(row.get("product_smiles_std") or row.get("product_smiles") or ""),
    }


def _select_ext_rules(meta_joined: pd.DataFrame, constraints: Dict[str, Any], k: int) -> List[Dict[str, Any]]:
    df = meta_joined.copy()
    df["step_type"] = df.get("step_type", "").astype(str)
    exts = df[df["step_type"].str.upper().eq("EXT")].copy()

    # AT substrate filter
    allowed_at = set(str(constraints.get("core_pathway", {}).get("at_substrate", "mal")).split("|"))
    if "at_substrate" in exts.columns:
        exts["at_norm"] = exts["at_substrate"].astype(str).str.lower().str.replace("_", "").str.replace("-", "")
        def _ok(x: str) -> bool:
            x = (x or "").lower()
            return any(tok in x for tok in allowed_at)
        exts = exts[exts["at_norm"].apply(_ok)]

    # Prefer diversity across BGCs, lowest module index first when available
    sort_cols = []
    if "module_from" in exts.columns:
        sort_cols.append("module_from")
    if "module_to" in exts.columns:
        sort_cols.append("module_to")
    if sort_cols:
        exts = exts.sort_values(sort_cols, na_position='last')

    picked: List[Dict[str, Any]] = []
    seen_bgcs: set[str] = set()
    for _, r in exts.iterrows():
        if len(picked) >= k:
            break
        bgc = str(r.get("bgc_id", ""))
        if bgc in seen_bgcs:
            continue
        picked.append({
            "reaction_id": str(r.get("reaction_id", "")),
            "bgc_id": bgc,
            "gml_path": str(r.get("gml_path", "")),
        })
        seen_bgcs.add(bgc)
    # If not enough, fill regardless of BGC diversity
    if len(picked) < k:
        remainder = []
        for _, r in exts.iterrows():
            item = {
                "reaction_id": str(r.get("reaction_id", "")),
                "bgc_id": str(r.get("bgc_id", "")),
                "gml_path": str(r.get("gml_path", "")),
            }
            if item not in picked:
                remainder.append(item)
        picked.extend(remainder[: max(0, k - len(picked))])
    return picked[:k]


def main() -> None:
    print("[START] Cross‑BGC rule‑level planner (metadata‑driven)")

    # Load inputs
    constraints = _load_constraints()
    target_smi = _load_target_smiles(constraints)
    n_ext_req = int(constraints.get("core_pathway", {}).get("required_n_ext", 5))

    # Read data tables
    if not RULEMETA.exists():
        raise FileNotFoundError(f"Missing rulemeta.csv: {RULEMETA}")
    if not RULES_IDX.exists():
        raise FileNotFoundError(f"Missing rules_index.csv: {RULES_IDX}")

    meta = pd.read_csv(RULEMETA)
    rix  = pd.read_csv(RULES_IDX)
    for df in (meta, rix):
        if "reaction_id" in df.columns:
            df["reaction_id"] = df["reaction_id"].astype(str)

    meta_j = _join_rule_index(meta, rix)

    # Pick closure rule
    closure = _pick_closure(meta_j, target_smi)
    print(f"[INFO] Chosen closure: {closure['reaction_id']} (bgc={closure['bgc_id']}, tanimoto={closure['tanimoto']:.4f})")

    # Select EXT rules from the global pool (diverse bgc set)
    ext_rules = _select_ext_rules(meta_j, constraints, n_ext_req)
    print(f"[INFO] Selected {len(ext_rules)} EXT rules across BGCs")

    # Assemble a metadata route (no application yet)
    route_steps: List[Dict[str, Any]] = []
    for i, ext in enumerate(ext_rules, 1):
        route_steps.append({
            "order": i,
            "step_type": "EXT",
            "reaction_id": ext["reaction_id"],
            "bgc_id": ext["bgc_id"],
            "gml_path": ext["gml_path"],
        })
    route_steps.append({
        "order": n_ext_req + 1,
        "step_type": "CLOSURE",
        "reaction_id": closure["reaction_id"],
        "bgc_id": closure["bgc_id"],
        "gml_path": closure["gml_path"],
        "note": "PT (c2c7) + TE lactonization expected",
    })

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = {
        "mode": "cross_bgc_rule_level",
        "target": target_smi,
        "constraints": constraints,
        "chosen_closure": closure,
        "ext_rules": ext_rules,
        "route_steps": route_steps,
    }
    with open(OUT_ROUTE, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    # Write a brief audit
    with open(OUT_AUDIT, "w", encoding="utf-8") as f:
        f.write("# Cross‑BGC Rule‑Level Planner (metadata‑driven)\n\n")
        f.write(f"Target: {target_smi}\n\n")
        f.write("Chosen closure rule:\n")
        f.write(f"- reaction_id: {closure['reaction_id']}\n")
        f.write(f"- bgc_id     : {closure['bgc_id']}\n")
        f.write(f"- tanimoto   : {closure['tanimoto']:.4f}\n")
        f.write(f"- gml        : {closure['gml_path']}\n\n")
        f.write(f"Selected EXT rules (n={len(ext_rules)}; required={n_ext_req}):\n")
        for i, r in enumerate(ext_rules, 1):
            f.write(f"  {i}. {r['reaction_id']} (bgc={r['bgc_id']}) gml={r['gml_path']}\n")
        f.write("\nNotes:\n")
        f.write("- This is a metadata‑driven assembly across BGCs without rule application.\n")
        f.write("- Next step: apply rules (reverse) with SynKit Reactor to generate intermediates.\n")

    print(f"[DONE] Wrote {OUT_ROUTE}")
    print(f"[DONE] Wrote {OUT_AUDIT}")


if __name__ == "__main__":
    main()


