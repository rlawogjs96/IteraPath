#!/usr/bin/env python3
# 3_AAM/scripts/02_h_adjust_and_rules.py

from pathlib import Path
import argparse
import json
import sys
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph


def _ensure_local_utils():
    try:
        from rule_emit_utils import build_rule_gml, sanitize_rule_stem  # noqa: F401
        return
    except Exception:
        pass
    here = Path(__file__).resolve().parent
    if str(here) not in sys.path:
        sys.path.insert(0, str(here))


_ensure_local_utils()
from rule_emit_utils import build_rule_gml, sanitize_rule_stem  # type: ignore


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
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting
except ModuleNotFoundError:
    _ensure_local_syntemp()
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting


BASE = Path(__file__).resolve().parent
DATA = BASE.parent / "data"
IN_IR = DATA / "its_raw.manual.ndjson"
OUT_IH = DATA / "its_hfixed.ndjson"
RULE_DIR = DATA / "rules"
RULE_DIR.mkdir(parents=True, exist_ok=True)


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


def main():
    ap = argparse.ArgumentParser(description="Hydrogen adjust ITS and extract rules (GML).")
    ap.add_argument("--its-input", type=Path, default=IN_IR, help="Path to ITS ndjson (can be manual)")
    ap.add_argument("--out-hfixed", type=Path, default=OUT_IH)
    ap.add_argument("--rules-dir", type=Path, default=RULE_DIR)
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--max-radius", type=int, default=3)
    args = ap.parse_args()

    if not args.its_input.exists():
        raise SystemExit(f"Missing ITS input: {args.its_input}")

    # Load NDJSON records
    records = []
    with args.its_input.open("r", encoding="utf-8") as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except Exception:
                continue

    # Optional hydrogen adjustment if ITS is stored as (G,H,ITS) tuples
    sample_its = None
    for r in records:
        v = r.get("ITS") or r.get("ITSGraph")
        if v is not None:
            sample_its = v
            break
    its_is_tuple = isinstance(sample_its, tuple) and len(sample_its) == 3

    if its_is_tuple:
        hadj = ITSHAdjuster()
        try:
            records = hadj.process_graph_data_parallel(
                graph_data_list=records,
                column="ITS",
                n_jobs=args.n_jobs,
                verbose=0,
                max_hydrogen=4,
                get_priority_graph=True,
            )
        except Exception as e:
            print(f"[WARN] H-adjustment failed: {e}; continuing without.")

    # Emit rules r0..max-radius
    rule_index = []
    max_r = max(0, int(args.max_radius))
    for rec in records:
        rid = rec.get("reaction_id") or rec.get("R-id")
        if not rid:
            continue
        rsmi = rec.get("reactions")
        its_payload = rec.get("ITSGraph") or rec.get("ITS")
        its_tuple = _graph_tuple_from_payload(its_payload) if its_payload else None
        if its_tuple is None:
            continue
        # r0
        try:
            L0, R0, K0 = RuleExtraction.extract_reaction_rules(*its_tuple, extend=False, n_knn=1)
        except Exception:
            continue
        safe = sanitize_rule_stem(rid)
        for r in range(0, max_r + 1):
            try:
                if r == 0:
                    L, R, K = L0, R0, K0
                else:
                    L, R, K = RuleExtraction.extract_reaction_rules(*its_tuple, extend=True, n_knn=r)
                changed = RuleWriting.find_changed_nodes(L, R, attributes=["charge", "hcount", "element"])
                stem = f"{safe}__r{r}"
                # Write GML/JSON/meta
                gml = build_rule_gml(L, R, K, stem, changed)
                (args.rules_dir / f"{stem}.gml").write_text(gml, encoding="utf-8")
                rule_json = {
                    "rule_id": stem,
                    "left": json_graph.node_link_data(L),
                    "right": json_graph.node_link_data(R),
                    "context": json_graph.node_link_data(K),
                }
                (args.rules_dir / f"{stem}.json").write_text(json.dumps(rule_json, ensure_ascii=False), encoding="utf-8")
                meta = {
                    "reaction_id": rid,
                    "reactions": rsmi,
                    "num_nodes": K.number_of_nodes(),
                    "num_left_edges": L.number_of_edges(),
                    "num_right_edges": R.number_of_edges(),
                    "changed_nodes": changed,
                }
                (args.rules_dir / f"{stem}.meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
                rule_index.append({"reaction_id": rid, "file_stem": stem})
            except Exception:
                continue

    # Index
    if rule_index:
        import pandas as pd  # local import

        pd.DataFrame(rule_index).to_csv(DATA / "rules_index.csv", index=False)
    print(f"[OK] rules saved under {args.rules_dir}")


if __name__ == "__main__":
    main()
