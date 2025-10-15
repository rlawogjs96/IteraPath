from pathlib import Path
import pandas as pd
import ast

BASE     = Path(__file__).resolve().parent
IN_RULES = BASE.parent.parent / "3_templates" / "data" / "rules_index.csv"
IN_META  = BASE.parent.parent / "1_preprocessing" / "data" / "syntemp_input_meta.csv"
IN_BAL   = BASE.parent.parent / "2_balancing" / "data" / "balanced_only.csv"
OUT_PATH = BASE.parent / "data" / "rulemeta.csv"
IN_AAM   = BASE.parent.parent / "3_templates" / "data" / "aam.ndjson"

def _safe_list(v):
    if pd.isna(v):
        return []
    if isinstance(v, list):
        return v
    s = str(v).strip()
    try:
        x = ast.literal_eval(s)
        if isinstance(x, list):
            return [str(i) for i in x]
    except Exception:
        pass
    return [p.strip() for p in s.split(",") if p.strip()]

def load_csv(p: Path, name: str) -> pd.DataFrame:
    if not p.exists():
        raise FileNotFoundError(f"{name} not found: {p}")
    df = pd.read_csv(p)
    # reaction_id 표준화
    if "reaction_id" in df.columns:
        df["reaction_id"] = df["reaction_id"].astype(str)
    return df

def first_nonnull(*vals):
    for v in vals:
        if v is None:
            continue
        s = str(v).strip()
        if s and s.lower() != "nan":
            return s
    return ""

def main():
    print("[START] Consolidating rule metadata...")
    print(f"[INFO] Loading rules_index from: {IN_RULES.resolve()}")
    print(f"[INFO] Loading input_meta  from: {IN_META.resolve()}")
    print(f"[INFO] Loading balanced    from: {IN_BAL.resolve()}")

    rules = load_csv(IN_RULES, "rules_index")
    # Parse variant index from file_stem suffix (__rN); keep ALL rows (one per rule file)
    if "file_stem" in rules.columns:
        rules["variant"] = (
            rules["file_stem"].astype(str).str.extract(r"__r(\d+)$").astype("Int64")
        )
        rules["variant"] = rules["variant"].fillna(0).astype(int)
    meta  = load_csv(IN_META,  "input_meta")
    bal   = load_csv(IN_BAL,   "balanced_only")

    # Diagnosis: original row counts / ID sets
    rid_rules = set(rules["reaction_id"]) if "reaction_id" in rules.columns else set()
    rid_meta  = set(meta["reaction_id"]) if "reaction_id" in meta.columns else set()
    rid_bal   = set(bal["reaction_id"]) if "reaction_id" in bal.columns else set()

    print(f"[DIAG] rows: rules={len(rules)}, meta={len(meta)}, balanced={len(bal)}")
    print(f"[DIAG] unique reaction_id: rules={len(rid_rules)}, meta={len(rid_meta)}, balanced={len(rid_bal)}")

    # Baseline: reaction_id in rules_index is the reference (must not be dropped)
    # Merge
    base = rules.copy()

    meta_ren = meta.add_suffix("_im")
    base = base.merge(meta_ren, left_on="reaction_id", right_on="reaction_id_im", how="left")
    if "reaction_id_im" in base.columns:
        base = base.drop(columns=["reaction_id_im"])

    # Check if balanced column exists before merging
    # Adapt to current balancing headers
    # Prefer standardized names if present
    bal_cols = ["reaction_id"]
    if "reactant_smiles_std" in bal.columns:
        bal_cols += ["reactant_smiles_std","product_smiles_std"]
    elif set(["reactant_smiles","product_smiles"]).issubset(bal.columns):
        bal_cols += ["reactant_smiles","product_smiles"]
    if "balanced" in bal.columns:
        bal_cols.append("balanced")
    base = base.merge(bal[bal_cols], on="reaction_id", how="left")

    # Optional: merge selected mapper and AAM from AAM NDJSON if present
    if IN_AAM.exists():
        try:
            aam_df = pd.read_json(IN_AAM, lines=True)
            keep_cols = [c for c in ["reaction_id","selected_mapper","aam"] if c in aam_df.columns]
            if keep_cols:
                base = base.merge(aam_df[keep_cols], on="reaction_id", how="left")
        except Exception:
            pass

    # Normalize bgc_id: first valid value from rules / meta / balanced
    def unify_bgc(row):
        return first_nonnull(
            row.get("bgc_id"),
            row.get("bgc"),                     # In case of different name
            row.get("bgc_id_im"),
        )

    # Normalize domains
    # Note: rules_index.csv has 'domain' (singular), input_meta has 'domains' (plural)
    def unify_domains(row):
        d_rules = _safe_list(row.get("domain"))  # ← 'domain' (singular)
        d_meta  = _safe_list(row.get("domains_im"))
        dd = sorted(set([*d_rules, *d_meta]))
        return ",".join(dd) if dd else ""

    # substrate / at_substrate : at_substrate takes priority
    def unify_substrate(row):
        at = str(row.get("at_substrate_im") or "").strip()
        if at:
            return at
        sub = str(row.get("substrate") or "").strip()
        return sub

    # step_type: meta takes priority, otherwise infer from rule_type
    def unify_step_type(row):
        st = str(row.get("step_type_im") or "").strip()
        if st:
            return st
        rt = str(row.get("rule_type") or "").upper().strip()
        if rt in {"EXT","CLOSURE"}:
            return rt
        return "UNKNOWN"

    # closure: closure_im from input_meta takes priority, otherwise closure from rules
    def unify_closure(row):
        c_im = str(row.get("closure_im") or "").strip()
        if c_im and c_im.lower() != "nan":
            return c_im
        c = str(row.get("closure") or "").strip()
        return c

    # origin: default to iPKS_bacterial if missing
    def unify_origin(row):
        o = str(row.get("origin") or "").strip()
        return o if o else "iPKS_bacterial"

    base["bgc_id"]       = base.apply(unify_bgc, axis=1)
    base["domains_norm"] = base.apply(unify_domains, axis=1)
    base["at_substrate"] = base.apply(unify_substrate, axis=1)
    base["step_type"]    = base.apply(unify_step_type, axis=1)
    base["closure"]      = base.apply(unify_closure, axis=1)
    base["origin"]       = base.apply(unify_origin, axis=1)

    # If CLOSURE but closure is empty, correct it
    mask_clo = base["step_type"].str.upper() == "CLOSURE"
    base.loc[mask_clo & (base["closure"].eq("") | base["closure"].str.lower().eq("nan")), "closure"] = "PT_OR_TE"

    # ========== PT mode inference (manual mapping until GML parsing is ready) ==========
    PT_MAP = {
        # isocoumarin family
        "BGC1000000": "c2c7",  # icmM (6-hydroxy"8-hydroxy-3-methyl-3,4-dihydroisochromen-1-one")
        "BGC1000001": "c2c7",  # SACE_5532
        "BGC1000006": "c2c7",  # NcsB (naphthocyclinone)
        # tetralone family
        "BGC1000002": "c2c7",  # AziB (azinomycin)
        "BGC1000005": "c2c7",  # ChlB1 (chlorothricin)
        # other aromatic
        "BGC1000003": "c2c7",  # PorKM1
        "BGC1000004": "c2c7",  # CalO5 (calicheamicin)
        "BGC1000007": "c2c7",  # AviM (avilamycin)
    }
    
    def infer_pt_mode(row):
        """Infer PT mode for CLOSURE steps (manual mapping until GML parsing)"""
        if str(row.get("step_type", "")).upper() != "CLOSURE":
            return ""
        bgc = str(row.get("bgc_id", "")).strip()
        return PT_MAP.get(bgc, "unknown")
    
    base["pt_mode"] = base.apply(infer_pt_mode, axis=1)
    
    # Strict validation: fail if any CLOSURE has unknown pt_mode
    unknown_pt = base[(base["step_type"] == "CLOSURE") & (base["pt_mode"] == "unknown")]
    if not unknown_pt.empty:
        print("[ERROR] CLOSURE steps with unknown pt_mode detected:")
        for _, row in unknown_pt.iterrows():
            print(f"  - {row['reaction_id']} (BGC: {row.get('bgc_id', 'N/A')})")
        raise AssertionError(
            "CLOSURE step(s) with unknown pt_mode detected. "
            "Please extend PT_MAP in 00_consolidate_rule_metadata.py before proceeding."
        )
    
    print(f"[INFO] PT mode assigned to {int((base['pt_mode'] != '').sum())} CLOSURE steps")

    # Clean up output columns
    # module_from, module_to may have _im suffix
    out_cols = [
        "reaction_id","bgc_id",
        "module_from","module_from_im","module_to","module_to_im",
        "rule_type","step_type",
        "domains_norm","at_substrate","closure","pt_mode","origin",
        "reactant_smiles_std","product_smiles_std","reactant_smiles","product_smiles","balanced",
        "domain","substrate","domains_im","at_substrate_im","file_stem",
        "variant","selected_mapper","aam"
    ]
    # Select only existing columns
    existing_cols = [c for c in out_cols if c in base.columns]
    out = base[existing_cols].copy()
    
    # Convert module_from_im → module_from, module_to_im → module_to
    if "module_from_im" in out.columns and "module_from" not in out.columns:
        out["module_from"] = out["module_from_im"]
    if "module_to_im" in out.columns and "module_to" not in out.columns:
        out["module_to"] = out["module_to_im"]
    
    # Final cleanup: remove _im suffix columns
    out = out[[c for c in out.columns if not c.endswith("_im")]].copy()

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, index=False)

    # Summary (including missing values)
    st_counts = out["step_type"].value_counts(dropna=False).to_dict() if "step_type" in out.columns else {}
    print(f"[DONE] Wrote: {OUT_PATH.resolve()}")
    print(f"       Rows: {len(out)}")
    # Count bgc_id missing values too
    if "bgc_id" in out.columns:
        print(f"       Unique BGCs: {out['bgc_id'].nunique(dropna=True)} (missing={int(out['bgc_id'].isna().sum())})")
    else:
        print("       Unique BGCs: NA")
    print("       Step types (incl. UNKNOWN):", st_counts)

    # Output missing diagnostics
    rid_out = set(out["reaction_id"])
    missing_from_meta = sorted(rid_rules - rid_meta)
    missing_from_bal  = sorted(rid_rules - rid_bal)
    missing_in_out    = sorted(rid_rules - rid_out)
    if missing_from_meta:
        print(f"[WARN] reaction_id present in rules_index but missing in input_meta: {missing_from_meta}")
    if missing_from_bal:
        print(f"[WARN] reaction_id present in rules_index but missing in balanced_only: {missing_from_bal}")
    if missing_in_out:
        print(f"[WARN] reaction_id present in rules_index but missing in OUTPUT: {missing_in_out}")

if __name__ == "__main__":
    main()
