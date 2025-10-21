#!/usr/bin/env python3
# 3_templates/scripts/01_aam_ensemble_consensus.py
from pathlib import Path
import json
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import sys
"""
Currently not working because the ML toolkits RXNMapper and other ones are chem based not bio synthetic based. 
Trained on the wrong dataset. 
"""
# SynAAM
# Ensure local SynTemp package is importable when running from repo
def _ensure_local_syntemp():
    try:
        import syntemp  # noqa: F401
        return
    except ModuleNotFoundError:
        pass
    here = Path(__file__).resolve()
    # repo_root ≈ .../IteraPath
    repo_root = here.parents[2]
    local_pkg = repo_root / "SynTemp"
    if local_pkg.exists():
        sys.path.insert(0, str(local_pkg))

try:
    from syntemp.SynAAM.atom_map_consensus import AAMConsensus
except ModuleNotFoundError:
    _ensure_local_syntemp()
    from syntemp.SynAAM.atom_map_consensus import AAMConsensus
# Support both class name variants across SynTemp versions
try:
    from syntemp.SynAAM.aam_postprocess import AAMPostprocessor as _Post
except ImportError:
    from syntemp.SynAAM.aam_postprocess import AMMPostprocessor as _Post

# SynITS
from syntemp.SynITS.its_extraction import ITSExtraction

BASE   = Path(__file__).resolve().parent
DATA   = BASE.parent / "data"
INCSV  = DATA / "reactions.csv"
OUT_AAM= DATA / "aam.ndjson"
OUT_IR = DATA / "its_raw.ndjson"  # ITS before H-correction, with isomorphism/representative selection applied

# Configuration
MAPPERS_REQUESTED = ["rxn_mapper", "graphormer", "local_mapper"]  # Add "rdt" if needed
REQUIRE_EQ = 2   # Minimum number of agreeing tools (e.g., 2/3)

def _to_json_safe(obj):
    """Convert objects (including networkx.Graph, tuples, lists, and dicts) into
    JSON-serializable structures. Graphs -> node-link JSON.
    """
    # Basic JSON-serializable types
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    # NetworkX Graph types (Graph, DiGraph, MultiGraph, MultiDiGraph)
    if isinstance(obj, (nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph)):
        try:
            return json_graph.node_link_data(obj)
        except Exception:
            return str(obj)
    # Tuples/lists: convert recursively
    if isinstance(obj, (list, tuple)):
        return [_to_json_safe(x) for x in obj]
    # Sets
    if isinstance(obj, set):
        return [_to_json_safe(x) for x in obj]
    # Dicts: ensure values are serializable; convert keys to str just in case
    if isinstance(obj, dict):
        return {str(k): _to_json_safe(v) for k, v in obj.items()}
    # Fallback
    try:
        json.dumps(obj)
        return obj
    except Exception:
        return str(obj)

def _ensure_reactions_csv() -> pd.DataFrame:
    """Ensure data/reactions.csv exists with columns [reaction_id, reactions].
    If missing, attempt to build it from 2_balancing/data/balance_results.csv.
    """
    if INCSV.exists():
        df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)
        cols = set(df.columns)
        if {"reaction_id","reactions"}.issubset(cols):
            return df[["reaction_id","reactions"]].copy()
        if {"reaction_id","reactant_smiles","product_smiles"}.issubset(cols):
            tmp = df[["reaction_id","reactant_smiles","product_smiles"]].copy()
            tmp["reactions"] = tmp["reactant_smiles"].astype(str) + ">>" + tmp["product_smiles"].astype(str)
            return tmp[["reaction_id","reactions"]]
        raise SystemExit(f"{INCSV} exists but missing required columns [reaction_id, reactions] or split SMILES.")

    # Try to build from 2_balancing output (balance_results.csv has the augmented reactions)
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
    raise SystemExit("Missing reactions input. Provide 3_templates/data/reactions.csv or run 2_balancing to produce balance_results.csv.")

def main():
    # Prepare input reactions
    df = _ensure_reactions_csv()
    # sanitize
    df = df.dropna(subset=["reaction_id","reactions"]).copy()
    df = df[df["reactions"].astype(str).str.contains(">>")].drop_duplicates("reaction_id")

    # 1) Run AAM ensemble - detect available mapper wrappers
    detected = []
    for name in MAPPERS_REQUESTED:
        try:
            _tmp = AAMConsensus([], mappers=[name])
            _tmp.import_mapper(name)
            detected.append(name)
        except Exception:
            continue
    if not detected:
        raise SystemExit(f"[ERROR] No available mappers found from: {MAPPERS_REQUESTED}")
    print(f"[INFO] Using available mappers: {detected}")
    MAPPERS = detected
    cons = AAMConsensus([], mappers=MAPPERS)
    # batch_consensus: adds { mapper_name: mapped_rsmi, ... } fields to each reaction
    batch = [{"reaction_id": rid, "reactions": rsmi, "R-id": rid} for rid, rsmi in df.values]
    batch = cons.batch_consensus(batch)  # Internally calls each mapper wrapper
    
    # Create a lookup for original reaction data (ID -> SMILES)
    original_data = {rid: rsmi for rid, rsmi in df.values}

    # 2) AAM post-processing (first filtering by number set/atom count consistency check)
    # AMMPostprocessor is static; validate and annotate but keep all records
    batch = _Post.parallel_postprocess(batch, MAPPERS, threshold=len(MAPPERS), n_jobs=4, verbose=0)

    # 3) Mapper results → ITS/RC generation & isomorphism calculation
    #    ITSExtraction.parallel_process_smiles creates graphs per mapper and returns 'equivariant' via RC isomorphism
    try:
        proc_correct, proc_wrong = ITSExtraction.parallel_process_smiles(
            batch, MAPPERS, n_jobs=4, check_method="RC", export_full=False
        )
        proc = list(proc_correct) + list(proc_wrong)
    except Exception as e:
        print(f"[WARN] ITSExtraction.parallel_process_smiles failed in batch: {e}. Falling back to per-record.")
        proc = []
        for rec in batch:
            try:
                sub_correct, sub_wrong = ITSExtraction.parallel_process_smiles(
                    [rec], MAPPERS, n_jobs=1, check_method="RC", export_full=False
                )
                proc.extend(sub_correct)
                proc.extend(sub_wrong)
            except Exception:
                continue

    # 4) Ensemble consensus filtering: Only accept reactions where mappers agree
    #    According to SynTemp paper (supp lines 201-203): "Instances in which the tools disagreed 
    #    were considered unsuccessful, and no prediction of the AAM was returned."
    consensus_records = []
    rejected_count = 0
    for rec in proc:
        equivariant = rec.get("equivariant", 0)
        # equivariant indicates how many mappers produced isomorphic reaction centers
        if equivariant >= REQUIRE_EQ:
            consensus_records.append(rec)
        else:
            rejected_count += 1
    
    print(f"[INFO] Consensus filter: {len(consensus_records)} accepted, {rejected_count} rejected (equivariant < {REQUIRE_EQ})")
    
    # Only process reactions that passed consensus filtering and have ITS
    chosen_records = [rec for rec in consensus_records if rec.get("ITSGraph")]
    
    # Debug: Check reaction IDs in processed records
    print(f"[DEBUG] Processed {len(chosen_records)} records with ITS and consensus")
    if chosen_records:
        sample = chosen_records[0]
        print(f"[DEBUG] Sample record keys: {list(sample.keys())}")
        print(f"[DEBUG] Sample reaction_id: {sample.get('reaction_id')} / R-id: {sample.get('R-id')}")

    # 5) Write: AAM (selected from agreeing mappers) + ITS_raw
    DATA.mkdir(parents=True, exist_ok=True)
    
    # Preference order for selecting which agreeing mapper to use
    mapper_preference = [m for m in ["graphormer", "rxn_mapper", "local_mapper"] if m in MAPPERS]
    
    with OUT_AAM.open("w", encoding="utf-8") as fa, OUT_IR.open("w", encoding="utf-8") as fi:
        for rec in chosen_records:
            # Try multiple ways to get the reaction ID
            rid = rec.get("reaction_id") or rec.get("R-id", "N/A")
            # If still "N/A", try to match from original data using reaction SMILES
            if rid == "N/A":
                rec_rsmi = rec.get("reactions") or rec.get("rxn_mapper") or rec.get("graphormer", "")
                for orig_id, orig_smi in original_data.items():
                    if orig_smi in rec_rsmi or rec_rsmi in orig_smi:
                        rid = orig_id
                        break
            
            # Get original reaction SMILES from lookup, fallback to mapped versions
            rsmi = original_data.get(rid) or rec.get("reactions") or rec.get("rxn_mapper") or rec.get("graphormer")
            
            # AAM: Select from one of the agreeing mappers (those that passed consensus)
            # Use preference order to pick which agreeing mapper's AAM to use
            selected_map, chosen_aam = None, None
            for mapper_name in mapper_preference:
                aam_val = rec.get(mapper_name)
                if isinstance(aam_val, str) and ">>" in aam_val:
                    selected_map, chosen_aam = mapper_name, aam_val
                    break
            
            # Ensure chosen_aam is JSON serializable
            chosen_aam_serializable = _to_json_safe(chosen_aam)
                
            aam_record = _to_json_safe({
                "reaction_id": rid,
                "reactions": rsmi,
                "selected_mapper": selected_map,
                "aam": chosen_aam_serializable
            })
            fa.write(json.dumps(aam_record, ensure_ascii=False) + "\n")

            # ITS: serialize representative ITSGraph provided by its_extraction/its_arbitrary (save as string/dictionary format as-is)
            its_payload = rec.get("ITSGraph") or rec.get("ITS") or rec.get("ITS_list")
            if its_payload is None:
                continue
            
            # Convert Graph/tuple/list objects to JSON-safe structures
            its_serializable = _to_json_safe(its_payload)
                
            its_record = _to_json_safe({
                "reaction_id": rid,
                "reactions": rsmi,
                "ITS": its_serializable
            })
            fi.write(json.dumps(its_record, ensure_ascii=False) + "\n")

    print(f"[OK] wrote {OUT_AAM}")
    print(f"[OK] wrote {OUT_IR}")

if __name__ == "__main__":
    main()
