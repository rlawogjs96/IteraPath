# 5_planner/02_deterministic_planner.py
from pathlib import Path
import json
import pandas as pd

# Input/output paths (relative to 5_planner root)
SELECTOR = Path("./data/selector_out.json")
BALANCED = Path("../2_balancing/data/balanced_only.csv")
RULEMETA = Path("../4_rulemeta/data/rulemeta.csv")

OUTDIR   = Path("../cases/orthosporin")
OUT_ROUTE= OUTDIR / "route.json"
OUT_GENE = OUTDIR / "gene_spec.json"
OUT_AUD  = OUTDIR / "audit.md"

# Gate settings
ALLOWED_AT = {"mal", "malonyl"}
DISABLE_ER = True

def load_inputs():
    """Load input files"""
    with open(SELECTOR) as f:
        sel = json.load(f)
    
    bal = pd.read_csv(BALANCED)
    meta = pd.read_csv(RULEMETA)
    
    for df in (bal, meta):
        if "reaction_id" in df.columns:
            df["reaction_id"] = df["reaction_id"].astype(str)
    
    return sel, bal, meta

def enforce_gates(dfb: pd.DataFrame, meta: pd.DataFrame):
    """
    Gate validation and gene specification generation
    
    Gates:
    1. ER prohibited (bacterial aromatic iPKS)
    2. EXT: KS+AT required, AT=mal verification
    3. CLOSURE: closure metadata presence check
    4. (Optional) balanced flag check
    
    Returns:
        (ok: bool, logs: list, gene_spec: dict or None)
    """
    logs = []
    gene_domains = set()
    at_used = set()
    pt_present = False
    te_present = False
    
    # Merge metadata
    df = dfb.merge(
        meta.drop_duplicates("reaction_id"),
        on="reaction_id",
        how="left",
        suffixes=("", "_meta")
    )
    
    # Sort by module_from, module_to
    sort_cols = []
    if "module_from" in df.columns:
        sort_cols.append("module_from")
    if "module_to" in df.columns:
        sort_cols.append("module_to")
    
    if sort_cols:
        df = df.sort_values(sort_cols, na_position='last').reset_index(drop=True)
    
    # Check each reaction
    for i, row in df.iterrows():
        rid = row["reaction_id"]
        step = row.get("step_type") or row.get("step_type_meta") or "UNKNOWN"
        domains = str(row.get("domains") or row.get("domains_meta") or "")
        at_sub = str(row.get("at_substrate") or row.get("at_substrate_meta") or "").lower().strip()
        closure = str(row.get("closure") or row.get("closure_meta") or "")
        
        # ER prohibition
        if DISABLE_ER and "ER" in domains.upper():
            logs.append(f"[FAIL] {rid}: ER domain found (domains={domains})")
            return False, logs, None
        
        # EXT gate
        if step == "EXT":
            if "KS" not in domains.upper() or "AT" not in domains.upper():
                logs.append(f"[FAIL] {rid}: EXT without KS+AT (domains={domains})")
                return False, logs, None
            
            if at_sub and at_sub not in ALLOWED_AT:
                logs.append(f"[FAIL] {rid}: AT substrate '{at_sub}' not in {ALLOWED_AT}")
                return False, logs, None
            
            gene_domains.update(["KS", "AT"])
            if at_sub:
                at_used.add(at_sub)
        
        # Reduction gate (simplified)
        if "KR" in domains.upper():
            gene_domains.add("KR")
        if "DH" in domains.upper():
            gene_domains.add("DH")
        
        # Termination gate
        if step == "CLOSURE":
            if not closure or closure.lower() in ["nan", "none", ""]:
                logs.append(f"[FAIL] {rid}: CLOSURE without closure metadata")
                return False, logs, None
            
            val = closure.upper()
            if "PT" in val:
                pt_present = True
                gene_domains.add("PT")
            if "TE" in val:
                te_present = True
                gene_domains.add("TE")
        
        # (Optional) balanced flag check
        if "balanced" in df.columns:
            b = row.get("balanced")
            if pd.notna(b) and str(b).lower() == "false":
                logs.append(f"[WARN] {rid}: balanced=False in input")
        
        logs.append(f"[OK] {rid}: step={step}, domains={domains or '-'}, at={at_sub or '-'}, closure={closure or '-'}")
    
    # Final sanity check
    if not pt_present and not te_present:
        logs.append("[WARN] No PT/TE observed in pathway (check CLOSURE rules)")
    
    # Generate gene specification
    gene_spec = {
        "required_domains": sorted(gene_domains),
        "AT_substrates": sorted(at_used) if at_used else ["mal"],
        "closure": "PT (C2-C7 aldol) + TE (lactonization/release)",
        "notes": [
            "bacterial_aromatic_iPKS",
            "ER_disabled",
            "hexaketide_core (acetyl + 5×malonyl)",
            "PT_mode: C2-C7 (preferred for isocoumarin)"
        ]
    }
    
    return True, logs, gene_spec

def main():
    print("[START] Deterministic planner...")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    
    # Load inputs
    sel, bal, meta = load_inputs()
    
    if not sel.get("top1"):
        raise RuntimeError("selector_out.json has no top1. Run 01_mcs_template_selector.py first.")
    
    chosen = sel["top1"]["bgc_id"]
    print(f"[INFO] Selected BGC: {chosen}")
    
    # Filter selected BGC pathway
    dfb = bal[bal["bgc_id"] == chosen].copy()
    if dfb.empty:
        raise RuntimeError(f"No records for bgc_id={chosen} in balanced_only.csv")
    
    print(f"[INFO] Found {len(dfb)} reactions for {chosen}")
    
    # Gate validation
    ok, logs, gene_spec = enforce_gates(dfb, meta)
    
    # ▼▼▼ Added: Force EXT×5 & PT ▼▼▼
    # Sort (same criteria as route output)
    sort_cols = []
    if "module_from" in dfb.columns:
        sort_cols.append("module_from")
    if "module_to" in dfb.columns:
        sort_cols.append("module_to")
    
    if sort_cols:
        dfb_sorted = dfb.sort_values(sort_cols, na_position='last').reset_index(drop=True)
    else:
        dfb_sorted = dfb.copy()
    
    # Force EXT count (Orthosporin: hexaketide = acetyl + 5×malonyl)
    force_n_ext = 5
    ext_count = 0
    if "step_type" in dfb_sorted.columns:
        ext_count = int((dfb_sorted["step_type"] == "EXT").sum())
    
    if ext_count < force_n_ext:
        ok = False
        logs.append(f"[FAIL] EXT count {ext_count} < required {force_n_ext} (hexaketide required)")
    else:
        logs.append(f"[OK] EXT count {ext_count} >= required {force_n_ext}")
    
    # Force PT (C2-C7 aldol cyclization for isocoumarin)
    has_pt = False
    if "closure" in dfb_sorted.columns:
        has_pt = dfb_sorted["closure"].fillna("").astype(str).str.upper().str.contains("PT").any()
    # Fallback: also check domains
    if not has_pt and "domains" in dfb_sorted.columns:
        has_pt = dfb_sorted["domains"].fillna("").astype(str).str.upper().str.contains("PT").any()
    
    if not has_pt:
        ok = False
        logs.append("[FAIL] PT not found in pathway (PT required for isocoumarin C2-C7 closure)")
    else:
        logs.append("[OK] PT present (C2-C7 aldol cyclization assumed; TE co-occurs for release)")
    # ▲▲▲ End of additions ▲▲▲
    
    # Generate route.json
    sort_cols = []
    if "module_from" in dfb.columns:
        sort_cols.append("module_from")
    if "module_to" in dfb.columns:
        sort_cols.append("module_to")
    
    if sort_cols:
        dfb = dfb.sort_values(sort_cols, na_position='last')
    
    # Select columns to output
    step_cols = ["reaction_id"]
    for c in ["step_type", "reactant_smiles", "product_smiles", "domains", "at_substrate"]:
        if c in dfb.columns:
            step_cols.append(c)
    
    route = {
        "bgc_id": chosen,
        "rules_seq": dfb["reaction_id"].astype(str).tolist(),
        "steps": dfb[step_cols].fillna("").to_dict(orient="records"),
        "target": sel.get("target"),
        "selector_metrics": sel["top1"],
        "status": "PASS" if ok else "FAIL"
    }
    
    with open(OUT_ROUTE, "w") as f:
        json.dump(route, f, indent=2)
    
    # Generate gene_spec.json
    if gene_spec:
        with open(OUT_GENE, "w") as f:
            json.dump(gene_spec, f, indent=2)
    
    # Generate audit.md
    with open(OUT_AUD, "w", encoding="utf-8") as f:
        f.write("# Deterministic Planner Audit\n\n")
        f.write(f"- Selected BGC: {chosen}\n")
        f.write(f"- Target: {sel.get('target')}\n")
        f.write(f"- Selector metrics:\n")
        f.write(f"  - Tanimoto: {sel['top1'].get('tanimoto')}\n")
        f.write(f"  - MCS atoms: {sel['top1'].get('mcs_atoms')}\n\n")
        f.write("## Gate Logs\n\n")
        for line in logs:
            f.write(line + "\n")
        f.write("\n## Result\n\n")
        if ok:
            f.write("[RESULT] **PASS** - All gates satisfied\n")
        else:
            f.write("[RESULT] **FAIL** - See logs above\n")
    
    print(f"[DONE] Wrote:")
    print(f"       - {OUT_ROUTE}")
    print(f"       - {OUT_GENE}")
    print(f"       - {OUT_AUD}")
    print(f"       Status: {'PASS' if ok else 'FAIL'}")

if __name__ == "__main__":
    main()
