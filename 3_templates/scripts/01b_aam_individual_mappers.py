#!/usr/bin/env python3
# 3_templates/scripts/01b_aam_individual_mappers.py
#
# Run multiple AAM tools independently and save each mapper's results separately.
# This allows manual inspection and selection of the best mapping for each reaction,
# especially useful for specialized chemistry where consensus may fail.
#
from pathlib import Path
import json
import pandas as pd
import sys

# SynAAM
def _ensure_local_syntemp():
    try:
        import syntemp  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    here = Path(__file__).resolve()
    repo_root = here.parents[2]
    local_pkg = repo_root / "SynTemp"
    if local_pkg.exists():
        sys.path.insert(0, str(local_pkg))

try:
    from syntemp.SynAAM.atom_map_consensus import AAMConsensus
except ModuleNotFoundError:
    _ensure_local_syntemp()
    from syntemp.SynAAM.atom_map_consensus import AAMConsensus

BASE = Path(__file__).resolve().parent
DATA = BASE.parent / "data"
INCSV = DATA / "reactions.csv"
OUT_DIR = DATA / "individual_mappers"

# Configuration
# Using only RDT (rule-based MCS) because ML mappers strip cofactors
# RDT sees full multi-component reactions and performed best on biochemical data
MAPPERS_REQUESTED = ["rdt"]

def _ensure_reactions_csv() -> pd.DataFrame:
    """Ensure data/reactions.csv exists with columns [reaction_id, reactions]."""
    if INCSV.exists():
        df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)
        cols = set(df.columns)
        if {"reaction_id","reactions"}.issubset(cols):
            return df[["reaction_id","reactions"]].copy()
        if {"reaction_id","reactant_smiles","product_smiles"}.issubset(cols):
            tmp = df[["reaction_id","reactant_smiles","product_smiles"]].copy()
            tmp["reactions"] = tmp["reactant_smiles"].astype(str) + ">>" + tmp["product_smiles"].astype(str)
            return tmp[["reaction_id","reactions"]]
        raise SystemExit(f"{INCSV} exists but missing required columns.")
    
    # Try to build from 2_balancing output
    alt = (BASE.parent.parent / "2_balancing" / "data" / "balance_results.csv").resolve()
    if alt.exists():
        sdf = pd.read_csv(alt, dtype=str, keep_default_na=False)
        need = {"reaction_id","reactions"}
        miss = need - set(sdf.columns)
        if miss:
            raise SystemExit(f"Cannot construct reactions.csv: missing {miss} in {alt}")
        out = sdf[["reaction_id","reactions"]].drop_duplicates().copy()
        DATA.mkdir(parents=True, exist_ok=True)
        out.to_csv(INCSV, index=False)
        print(f"[INFO] Built {INCSV} from {alt} ({len(out)} reactions)")
        return out
    raise SystemExit("Missing reactions input. Provide 3_templates/data/reactions.csv or run 2_balancing.")

def main():
    # Prepare input reactions
    df = _ensure_reactions_csv()
    df = df.dropna(subset=["reaction_id","reactions"]).copy()
    df = df[df["reactions"].astype(str).str.contains(">>")].drop_duplicates("reaction_id")
    
    print(f"[INFO] Processing {len(df)} reactions")
    
    # Detect available mappers
    detected = []
    for name in MAPPERS_REQUESTED:
        try:
            _tmp = AAMConsensus([], mappers=[name])
            _tmp.import_mapper(name)
            detected.append(name)
            print(f"[INFO] Detected mapper: {name}")
        except Exception as e:
            print(f"[WARN] Mapper {name} not available: {e}")
            continue
    
    if not detected:
        raise SystemExit(f"[ERROR] No available mappers found from: {MAPPERS_REQUESTED}")
    
    MAPPERS = detected
    print(f"[INFO] Using mappers: {MAPPERS}")
    
    # Prepare batch data
    batch = [{"reaction_id": rid, "reactions": rsmi, "R-id": rid} for rid, rsmi in df.values]
    
    # Handle RDT separately (needs special parameters)
    if "rdt" in MAPPERS:
        from syntemp.SynAAM.rdt_wrapper import map_with_rdt
        import importlib.resources
        
        rdt_jar = str(importlib.resources.files("syntemp.SynAAM").joinpath("RDT_2.4.1.jar"))
        working_dir = str(DATA / "rdt_temp")
        
        print(f"[INFO] Running RDT mapper on {len(batch)} reactions...")
        for rec in batch:
            rsmi = rec["reactions"]
            try:
                mapped = map_with_rdt(rsmi, rdt_jar, working_dir)
                rec["rdt"] = mapped
            except Exception as e:
                print(f"[WARN] RDT failed for {rec['reaction_id']}: {e}")
                rec["rdt"] = rsmi
        
        # Remove rdt from MAPPERS list since we handled it manually
        MAPPERS = [m for m in MAPPERS if m != "rdt"]
    
    # Run other mappers through consensus if any remain
    if MAPPERS:
        cons = AAMConsensus([], mappers=MAPPERS)
        batch = cons.batch_consensus(batch)
    
    # Add rdt back for output processing
    if "rdt" in detected:
        MAPPERS = detected
    
    # Create output directory
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save results for each mapper separately
    for mapper_name in MAPPERS:
        output_file = OUT_DIR / f"{mapper_name}_mappings.csv"
        output_json = OUT_DIR / f"{mapper_name}_mappings.json"
        
        results = []
        json_results = []
        
        for rec in batch:
            rid = rec.get("reaction_id") or rec.get("R-id", "N/A")
            original_rsmi = rec.get("reactions", "")
            mapped_rsmi = rec.get(mapper_name, "")
            
            # Check if mapping succeeded
            success = isinstance(mapped_rsmi, str) and ">>" in mapped_rsmi
            
            results.append({
                "reaction_id": rid,
                "original_reaction": original_rsmi,
                "mapped_reaction": mapped_rsmi if success else "",
                "success": success
            })
            
            json_results.append({
                "reaction_id": rid,
                "original_reaction": original_rsmi,
                "mapped_reaction": mapped_rsmi if success else "",
                "success": success,
                "mapper": mapper_name
            })
        
        # Save CSV
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, index=False)
        
        # Save JSON (more detailed)
        with output_json.open("w", encoding="utf-8") as f:
            for result in json_results:
                f.write(json.dumps(result, ensure_ascii=False) + "\n")
        
        success_count = sum(1 for r in results if r["success"])
        print(f"[OK] {mapper_name}: {success_count}/{len(results)} successful -> {output_file}")
    
    # Create a comparison file showing where mappers agree/disagree
    comparison_file = OUT_DIR / "mapper_comparison.csv"
    comparison_results = []
    
    for rec in batch:
        rid = rec.get("reaction_id") or rec.get("R-id", "N/A")
        original_rsmi = rec.get("reactions", "")
        
        comp_row = {
            "reaction_id": rid,
            "original_reaction": original_rsmi
        }
        
        # Add each mapper's result
        mapper_results = {}
        for mapper_name in MAPPERS:
            mapped = rec.get(mapper_name, "")
            success = isinstance(mapped, str) and ">>" in mapped
            comp_row[f"{mapper_name}_mapped"] = mapped if success else "FAILED"
            comp_row[f"{mapper_name}_success"] = success
            if success:
                mapper_results[mapper_name] = mapped
        
        # Simple agreement check: all successful mappers gave same result
        if mapper_results:
            unique_mappings = set(mapper_results.values())
            comp_row["mappers_agree"] = len(unique_mappings) == 1
            comp_row["num_successful_mappers"] = len(mapper_results)
        else:
            comp_row["mappers_agree"] = False
            comp_row["num_successful_mappers"] = 0
        
        comparison_results.append(comp_row)
    
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df.to_csv(comparison_file, index=False)
    
    agree_count = sum(1 for r in comparison_results if r.get("mappers_agree", False))
    print(f"\n[OK] Comparison: {agree_count}/{len(comparison_results)} reactions with mapper agreement -> {comparison_file}")
    print(f"[INFO] All individual mapper results saved to: {OUT_DIR}")

if __name__ == "__main__":
    main()

