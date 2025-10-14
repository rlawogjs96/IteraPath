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
CONSTRAINTS = OUTDIR / "constraints.json"
CONSTRAINTS_FREEZE = OUTDIR / "constraints.freeze.json"

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

def load_constraints():
    """
    Load constraints from constraints.json
    Saves a frozen copy for reproducibility
    """
    if CONSTRAINTS.exists():
        with open(CONSTRAINTS, encoding="utf-8") as f:
            c = json.load(f)
        print(f"[INFO] Loaded constraints from {CONSTRAINTS}")
    else:
        # Default constraints (Orthosporin)
        c = {
            "target": {
                "name": "Orthosporin",
                "smiles": "CC(CC1=CC2=CC(=CC(=C2C(=O)O1)O)O)O"
            },
            "core_pathway": {
                "pks_type": "I_NR",
                "starter_unit": "acetyl-CoA",
                "at_substrate": "mal|malonyl",
                "required_n_ext": 5,
                "er_allowed": False,
                "kr_mode": "late_or_once",
                "dh": "optional"
            },
            "closure": {
                "require_pt": True,
                "pt_mode_allowed": ["c2c7"],
                "release": "TE_lactonization"
            },
            "tailoring": {
                "steps": []
            }
        }
        print("[INFO] Using default constraints (no constraints.json found)")
    
    # Save frozen copy for reproducibility
    with open(CONSTRAINTS_FREEZE, "w", encoding="utf-8") as fz:
        json.dump(c, fz, indent=2, ensure_ascii=False)
    print(f"[INFO] Saved frozen constraints to {CONSTRAINTS_FREEZE}")
    
    return c

def enforce_gates(dfb: pd.DataFrame, meta: pd.DataFrame, constraints: dict):
    """
    Gate validation and gene specification generation
    
    Uses constraints to determine:
    - allowed AT substrates
    - ER allowed/prohibited
    - required n_ext
    - PT/TE requirements
    
    Returns:
        (ok: bool, logs: list, gene_spec: dict or None)
    """
    # Extract constraints
    core = constraints.get("core_pathway", {})
    closure_spec = constraints.get("closure", {})
    tailoring = constraints.get("tailoring", {}).get("steps", [])
    
    # Parse constraints
    at_substrate_raw = core.get("at_substrate", "mal")
    allowed_at = set(at_substrate_raw.split("|"))
    er_allowed = bool(core.get("er_allowed", False))
    kr_mode = core.get("kr_mode", "late_or_once")
    dh_state = core.get("dh", "optional")
    pt_mode_list = closure_spec.get("pt_mode_allowed", ["c2c7"])
    release_type = closure_spec.get("release", "TE_lactonization")
    
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
        
        # ER prohibition check
        if not er_allowed and "ER" in domains.upper():
            logs.append(f"[FAIL] {rid}: ER domain found but er_allowed=False (domains={domains})")
            return False, logs, None
        
        # EXT gate
        if step == "EXT":
            if "KS" not in domains.upper() or "AT" not in domains.upper():
                logs.append(f"[FAIL] {rid}: EXT without KS+AT (domains={domains})")
                return False, logs, None
            
            if at_sub and at_sub not in allowed_at:
                logs.append(f"[FAIL] {rid}: AT substrate '{at_sub}' not in {allowed_at}")
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
    pt_mode_str = "/".join(pt_mode_list) if pt_mode_list else "unknown"
    gene_spec = {
        "required_domains": sorted(gene_domains),
        "AT_substrates": sorted(at_used) if at_used else sorted(allowed_at),
        "ER_allowed": er_allowed,
        "KR_mode": kr_mode,
        "DH": dh_state,
        "closure": f"PT ({pt_mode_str}) + {release_type}",
        "tailoring": tailoring,
        "notes": [
            f"pks_type: {core.get('pks_type', 'I_NR')}",
            f"ER: {'enabled' if er_allowed else 'disabled'}",
            f"PT_mode: {pt_mode_str}"
        ]
    }
    
    return True, logs, gene_spec

def main():
    print("[START] Deterministic planner...")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    
    # Load inputs
    sel, bal, meta = load_inputs()
    constraints = load_constraints()
    
    if not sel.get("top1"):
        raise RuntimeError("selector_out.json has no top1. Run 01_mcs_template_selector.py first.")
    
    chosen = sel["top1"]["bgc_id"]
    print(f"[INFO] Selected BGC: {chosen}")
    
    # Filter selected BGC pathway
    dfb = bal[bal["bgc_id"] == chosen].copy()
    if dfb.empty:
        raise RuntimeError(f"No records for bgc_id={chosen} in balanced_only.csv")
    
    print(f"[INFO] Found {len(dfb)} reactions for {chosen}")
    
    # Gate validation with constraints
    ok, logs, gene_spec = enforce_gates(dfb, meta, constraints)
    
    # ▼▼▼ Added: Force EXT count & PT validation ▼▼▼
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
    
    # Force EXT count from constraints
    force_n_ext = constraints.get("core_pathway", {}).get("required_n_ext", 5)
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
    
    # Build steps with mechanism annotation for CLOSURE
    steps = []
    for _, row in dfb.iterrows():
        step = {
            "reaction_id": str(row["reaction_id"]),
        }
        # Add available columns
        for col in ["step_type", "reactant_smiles", "product_smiles", "domains", "at_substrate"]:
            if col in dfb.columns:
                step[col] = str(row.get(col, ""))
        
        # Add mechanism annotation for CLOSURE steps
        if row.get("step_type") == "CLOSURE":
            pt_mode = "/".join(constraints.get("closure", {}).get("pt_mode_allowed", ["c2c7"]))
            release = constraints.get("closure", {}).get("release", "TE_lactonization")
            step["mechanism"] = f"PT ({pt_mode} aldol) + {release}"
        
        steps.append(step)
    
    route = {
        "bgc_id": chosen,
        "rules_seq": dfb["reaction_id"].astype(str).tolist(),
        "steps": steps,
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
