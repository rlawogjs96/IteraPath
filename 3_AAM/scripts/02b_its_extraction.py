#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys
import json
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
try:
    # local helper utilities for GML emission
    from rule_emit_utils import build_rule_gml, sanitize_rule_stem, write_rule_triplet
except Exception:
    # allow running from other CWDs
    from .rule_emit_utils import build_rule_gml, sanitize_rule_stem, write_rule_triplet  # type: ignore


def _ensure_local_syntemp():
    """Ensure local SynTemp package is importable when running from repo checkout."""
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
    from syntemp.SynITS.its_extraction import ITSExtraction
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting
except ModuleNotFoundError:
    _ensure_local_syntemp()
    from syntemp.SynITS.its_extraction import ITSExtraction
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting


BASE = Path(__file__).resolve().parent
DATA = BASE.parent / "data"
DEFAULT_IN = DATA / "manual_mapping.csv"
DEFAULT_OUT_AAM = DATA / "aam.manual.ndjson"
DEFAULT_OUT_ITS = DATA / "its_raw.manual.ndjson"
RULE_DIR = DATA / "rules"
ORDER_TO_LABEL = {1: "-", 1.5: ":", 2: "=", 3: "#"}


def _to_json_safe(obj):
    """Convert objects (incl. networkx graphs) to JSON-serializable structures."""
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, (nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph)):
        try:
            return json_graph.node_link_data(obj)
        except Exception:
            return str(obj)
    if isinstance(obj, (list, tuple)):
        return [_to_json_safe(x) for x in obj]
    if isinstance(obj, set):
        return [_to_json_safe(x) for x in obj]
    if isinstance(obj, dict):
        return {str(k): _to_json_safe(v) for k, v in obj.items()}
    try:
        json.dumps(obj)
        return obj
    except Exception:
        return str(obj)


def load_manual_aam(in_csv: Path, id_col: str, rsmi_col: str, aam_col: str) -> pd.DataFrame:
    if not in_csv.exists():
        raise SystemExit(f"Missing input CSV: {in_csv}")
    df = pd.read_csv(in_csv, dtype=str, keep_default_na=False)
    cols = set(df.columns)

    # Heuristics for common column names
    if id_col not in cols:
        for cand in ["reaction_id", "R-id", "id"]:
            if cand in cols:
                id_col = cand
                break
    if rsmi_col not in cols:
        for cand in ["reactions_original", "reactions", "rsmi", "reaction_smiles"]:
            if cand in cols:
                rsmi_col = cand
                break
    if aam_col not in cols:
        for cand in ["AAM", "mapped", "manual", "manual_aam"]:
            if cand in cols:
                aam_col = cand
                break

    missing = [c for c in [id_col, rsmi_col, aam_col] if c not in df.columns]
    if missing:
        raise SystemExit(
            f"Input {in_csv} missing required columns. Looked for id='{id_col}', rsmi='{rsmi_col}', aam='{aam_col}'."
        )

    # Normalize column names to standard keys
    out = df[[id_col, rsmi_col, aam_col]].rename(
        columns={id_col: "R-id", rsmi_col: "reactions", aam_col: "AAM"}
    )
    out = out.dropna(subset=["R-id", "reactions", "AAM"]).drop_duplicates("R-id")
    return out


def build_its_from_manual(
    df: pd.DataFrame,
    n_jobs: int = 4,
    check_method: str = "RC",
    ignore_aromaticity: bool = False,
    sanitize: bool = True,
):
    # Prepare batch for ITSExtraction: treat manual AAM as a single mapper named "manual"
    def _sanitize_reaction_aam(rsmi: str) -> str:
        # Remove stray whitespace and trailing dots on each side; remove empty fragments
        s = (rsmi or "").strip()
        if ">>" not in s:
            return s
        left, right = s.split(">>", 1)
        def fix_side(side: str) -> str:
            side = side.strip().strip('.')
            parts = [p.strip() for p in side.split('.') if p.strip()]
            return '.'.join(parts)
        left_f = fix_side(left)
        right_f = fix_side(right)
        if not left_f or not right_f:
            return s
        return f"{left_f}>>{right_f}"

    rows = []
    for _, r in df.iterrows():
        rows.append({"R-id": r["R-id"], "manual": _sanitize_reaction_aam(r["AAM"])})

    correct, wrong = ITSExtraction.parallel_process_smiles(
        mapped_smiles_list=rows,
        mapper_names=["manual"],
        n_jobs=n_jobs,
        check_method=check_method,
        export_full=False,
        ignore_aromaticity=ignore_aromaticity,
        confident_mapper="manual",
        sanitize=sanitize,
    )
    return correct, wrong


def hydrogen_adjust(
    results: list,
    n_jobs: int = 4,
    ignore_aromaticity: bool = False,
    balance_its: bool = True,
    safe: bool = True,
    job_timeout: int = 5,
):
    if not results:
        return []
    hadj = ITSHAdjuster()
    adjusted = hadj.process_graph_data_parallel(
        graph_data_list=results,
        column="ITSGraph",
        n_jobs=n_jobs,
        verbose=0,
        ignore_aromaticity=ignore_aromaticity,
        balance_its=balance_its,
        job_timeout=job_timeout,
        safe=safe,
        get_priority_graph=False,
    )
    return adjusted


def write_outputs(
    df: pd.DataFrame,
    correct: list,
    out_aam: Path,
    out_its: Path,
    write_aam: bool = True,
    write_its: bool = True,
):
    id_to_rsmi = {r["R-id"]: r["reactions"] for _, r in df.iterrows()}
    id_to_aam = {r["R-id"]: r["AAM"] for _, r in df.iterrows()}

    if write_aam:
        with out_aam.open("w", encoding="utf-8") as fa:
            for rid, rsmi in id_to_rsmi.items():
                rec = _to_json_safe(
                    {
                        "reaction_id": rid,
                        "reactions": rsmi,
                        "selected_mapper": "manual",
                        "aam": id_to_aam.get(rid),
                    }
                )
                fa.write(json.dumps(rec, ensure_ascii=False) + "\n")

    if write_its:
        with out_its.open("w", encoding="utf-8") as fi:
            for d in correct:
                rid = d.get("R-id")
                its_tuple = d.get("ITSGraph")
                if not its_tuple or not isinstance(its_tuple, tuple) or len(its_tuple) != 3:
                    continue
                rsmi = id_to_rsmi.get(rid)
                graph_rules = d.get("GraphRules")
                record = {
                    "reaction_id": rid,
                    "reactions": rsmi,
                    "ITSGraph": _to_json_safe(its_tuple),
                }
                if graph_rules and isinstance(graph_rules, tuple) and len(graph_rules) == 3:
                    record["GraphRules"] = _to_json_safe(graph_rules)
                fi.write(json.dumps(record, ensure_ascii=False) + "\n")


def _graph_tuple_from_payload(payload):
    if isinstance(payload, (list, tuple)) and len(payload) == 3:
        graphs = []
        for item in payload:
            if isinstance(item, dict):
                graphs.append(json_graph.node_link_graph(item))
            elif isinstance(item, nx.Graph):
                graphs.append(item)
            else:
                return None
        return tuple(graphs)
    return None


def emit_rules_from_its_json(its_path: Path, rules_dir: Path, index_path: Path, max_radius: int = 3):
    rules_dir.mkdir(parents=True, exist_ok=True)
    rule_index = []
    cnt = 0
    with its_path.open("r", encoding="utf-8") as f:
        for line in f:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            rid = rec.get("reaction_id") or rec.get("R-id") or "unknown"
            rsmi = rec.get("reactions")
            its_payload = rec.get("ITSGraph") or rec.get("ITS")
            its_tuple = _graph_tuple_from_payload(its_payload) if its_payload else None
            if its_tuple is None:
                continue
            safe_rid = sanitize_rule_stem(rid)
            # Emit r0..max_radius (increasing neighborhood around RC)
            max_r = max(0, int(max_radius))
            for radius in range(0, max_r + 1):
                try:
                    L, R, K = RuleExtraction.extract_reaction_rules(*its_tuple, extend=True, n_knn=radius)
                except Exception:
                    continue
                if not isinstance(L, nx.Graph) or not isinstance(R, nx.Graph) or not isinstance(K, nx.Graph):
                    continue
                stem = f"{safe_rid}__r{radius}"
                try:
                    changed_node_ids = RuleWriting.find_changed_nodes(L, R, attributes=["charge", "hcount", "element"])
                    # Write GML/JSON/meta using shared helper
                    write_rule_triplet(
                        L, R, K, stem, rules_dir, changed_node_ids, rid=rid, rsmi=rsmi
                    )
                    cnt += 1
                    rule_index.append({"reaction_id": rid, "file_stem": stem})
                except Exception:
                    continue

    # write index
    if rule_index:
        pd.DataFrame(rule_index).to_csv(index_path, index=False)
    return cnt, len(rule_index)


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Build ITS/RC graphs from manually mapped AAM CSV using SynTemp's ITS pipeline;"
            " optionally perform hydrogen completion."
        )
    )
    ap.add_argument("--input", type=Path, default=DEFAULT_IN)
    ap.add_argument("--id-column", type=str, default="reaction_id")
    ap.add_argument("--rsmi-column", type=str, default="reactions_original")
    ap.add_argument("--aam-column", type=str, default="AAM")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--check-method", choices=["RC", "ITS"], default="RC")
    ap.add_argument("--ignore-aromaticity", action="store_true")
    ap.add_argument("--no-sanitize", action="store_true")
    ap.add_argument("--no-hadjust", action="store_true")
    ap.add_argument("--hadjust-timeout", type=int, default=5)
    ap.add_argument("--out-aam", type=Path, default=DEFAULT_OUT_AAM)
    ap.add_argument("--out-its", type=Path, default=DEFAULT_OUT_ITS)
    ap.add_argument("--no-write-aam", action="store_true")
    ap.add_argument("--no-write-its", action="store_true")
    ap.add_argument("--emit-rules", action="store_true", help="Emit GML rules from ITS output")
    ap.add_argument("--its-input", type=Path, default=None, help="Override ITS input path for rule emission")
    ap.add_argument("--rules-dir", type=Path, default=RULE_DIR)
    ap.add_argument("--max-radius", type=int, default=3)
    args = ap.parse_args()

    df = load_manual_aam(args.input, args.id_column, args.rsmi_column, args.aam_column)
    print(f"[INFO] Loaded {len(df)} manual AAM records from {args.input}")

    correct, wrong = build_its_from_manual(
        df,
        n_jobs=args.n_jobs,
        check_method=args.check_method,
        ignore_aromaticity=args.ignore_aromaticity,
        sanitize=not args.no_sanitize,
    )
    print(f"[INFO] ITS/RC extraction: {len(correct)} accepted; {len(wrong)} rejected (pre-H adjust)")

    if not args.no_hadjust:
        correct = hydrogen_adjust(
            correct,
            n_jobs=args.n_jobs,
            ignore_aromaticity=args.ignore_aromaticity,
            balance_its=True,
            safe=True,
            job_timeout=args.hadjust_timeout,
        )
        # Filter to those that still have valid ITSGraph after adjustment
        kept = [d for d in correct if isinstance(d.get("ITSGraph", (None, None, None))[2], nx.Graph)]
        print(f"[INFO] Hydrogen adjustment complete: {len(kept)} retained with valid ITS")
        correct = kept

    DATA.mkdir(parents=True, exist_ok=True)
    write_outputs(
        df=df,
        correct=correct,
        out_aam=args.out_aam,
        out_its=args.out_its,
        write_aam=not args.no_write_aam,
        write_its=not args.no_write_its,
    )
    # Optional rule emission
    if args.emit_rules:
        its_in = args.its_input if args.its_input else args.out_its
        index_csv = DATA / "rules_index.csv"
        emitted, indexed = emit_rules_from_its_json(its_in, args.rules_dir, index_csv, max_radius=args.max_radius)
        print(f"[OK] Rules emitted: {emitted}; index: {index_csv}")
    print("\n" + "=" * 70)
    print("ITS FROM MANUAL AAM - SUMMARY")
    print("=" * 70)
    print(f"Input manual AAM:     {args.input}")
    print(f"Records loaded:       {len(df)}")
    print(f"ITS accepted:         {len(correct)}")
    if not args.no_write_aam:
        print(f"AAM output:           {args.out_aam}")
    if not args.no_write_its:
        print(f"ITS output:           {args.out_its}")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()


