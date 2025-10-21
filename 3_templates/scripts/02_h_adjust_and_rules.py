#!/usr/bin/env python3
# 3_templates/scripts/02_h_adjust_and_rules.py
from pathlib import Path
import json
import pandas as pd
import sys
import networkx as nx
from networkx.readwrite import json_graph

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
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting
except ModuleNotFoundError:
    _ensure_local_syntemp()
    from syntemp.SynITS.its_hadjuster import ITSHAdjuster
    from syntemp.SynRule.rules_extraction import RuleExtraction
    from syntemp.SynRule.rule_writing import RuleWriting

BASE   = Path(__file__).resolve().parent
DATA   = BASE.parent / "data"
IN_IR  = DATA / "its_raw.ndjson"
OUT_IH = DATA / "its_hfixed.ndjson"
RULE_DIR = DATA / "rules"
RULE_DIR.mkdir(parents=True, exist_ok=True)

def main():
    if not IN_IR.exists():
        raise SystemExit("Run 01_aam_ensemble_consensus.py first.")

    # 1) 로드
    records = []
    with IN_IR.open("r", encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            # Keep ITS as JSON data for H-adjuster and serialization
            # We'll convert to NetworkX graphs later when extracting rules
            records.append(rec)

    # 2) H-보정 - ITSHAdjuster는 인자 없이 인스턴스화
    hadj = ITSHAdjuster()
    try:
        # process_graph_data_parallel의 필수 인자들을 올바르게 전달
        fixed = hadj.process_graph_data_parallel(
            graph_data_list=records,
            column="ITS",
            n_jobs=4,
            verbose=0,
            max_hydrogen=4,
            get_priority_graph=True
        )
    except TypeError as e:
        print(f"[WARN] First attempt failed: {e}")
        # 폴백: 최소 필수 인자만으로 시도
        try:
            fixed = hadj.process_graph_data_parallel(
                records, "ITS", n_jobs=4, verbose=0
            )
        except Exception as e2:
            print(f"[ERROR] H-adjustment failed: {e2}")
            # 최후 폴백: H-보정 없이 진행
            fixed = records

    with OUT_IH.open("w", encoding="utf-8") as fo:
        for rec in fixed:
            # Ensure ITS data is JSON-serializable
            rec_copy = rec.copy()
            its_val = rec_copy.get("ITS") or rec_copy.get("ITSGraph")
            if its_val and isinstance(its_val, tuple):
                # Convert NetworkX graphs back to JSON format if needed
                try:
                    rec_copy["ITS"] = [
                        json_graph.node_link_data(its_val[0]),
                        json_graph.node_link_data(its_val[1]),
                        json_graph.node_link_data(its_val[2])
                    ]
                except:
                    pass  # Keep as is if conversion fails
            fo.write(json.dumps(rec_copy, ensure_ascii=False) + "\n")
    print(f"[OK] wrote {OUT_IH}")

    # 3) 룰 추출 - Use custom approach that generates rules from ITS graph directly
    rule_index = []
    print(f"[DEBUG] Processing {len(fixed)} records for rule extraction")
    
    for rec in fixed:
        rid = rec["reaction_id"]; rsmi = rec["reactions"]
        its = rec.get("ITS") or rec.get("ITSGraph")
        print(f"[DEBUG] Processing {rid}: ITS type={type(its)}, has_ITS={its is not None}")
        if its is None:
            print(f"[DEBUG] Skipping {rid}: no ITS data")
            continue
        try:
            # ITS should be a list of 3 JSON graph dictionaries
            if not (isinstance(its, list) and len(its) == 3):
                print(f"[WARN] Unexpected ITS format for {rid}: type={type(its)}, len={len(its) if hasattr(its, '__len__') else 'N/A'}")
                continue
            
            # Convert JSON graphs to NetworkX
            try:
                if isinstance(its[0], dict):
                    # Standard case: list of JSON graph dictionaries
                    reactants_graph = json_graph.node_link_graph(its[0])
                    products_graph = json_graph.node_link_graph(its[1])
                    its_graph = json_graph.node_link_graph(its[2])
                elif hasattr(its[0], 'nodes'):
                    # Already NetworkX graphs (shouldn't happen after our changes, but handle it)
                    reactants_graph, products_graph, its_graph = its
                else:
                    print(f"[WARN] Unexpected ITS element type for {rid}: {type(its[0])}")
                    continue
            except Exception as e:
                print(f"[WARN] Failed to convert JSON to NetworkX for {rid}: {e}")
                import traceback
                traceback.print_exc()
                continue
            
            # NEW APPROACH: Generate rule graphs directly from ITS graph
            # The ITS graph contains all information needed for the rule
            print(f"[DEBUG] Generating rule from ITS graph for {rid}: ITS={its_graph.number_of_nodes()} nodes, {its_graph.number_of_edges()} edges")
            
            # Create L graph (left/reactants): nodes from ITS, edges where order[0] > 0
            L = nx.Graph()
            for node, attrs in its_graph.nodes(data=True):
                # Get left-side element type from typesGH if available
                typesGH = attrs.get('typesGH')
                if typesGH and len(typesGH) >= 1:
                    left_type = typesGH[0]  # [element, aromatic, hcount, charge, neighbors]
                    node_attrs = {
                        'element': left_type[0],
                        'aromatic': left_type[1] if len(left_type) > 1 else False,
                        'hcount': left_type[2] if len(left_type) > 2 else 0,
                        'charge': left_type[3] if len(left_type) > 3 else 0,
                    }
                else:
                    node_attrs = {
                        'element': attrs.get('element', '*'),
                        'aromatic': attrs.get('aromatic', False),
                        'hcount': attrs.get('hcount', 0),
                        'charge': attrs.get('charge', 0),
                    }
                L.add_node(node, **node_attrs)
            
            for u, v, data in its_graph.edges(data=True):
                order = data.get('order', [0, 0])
                if isinstance(order, (list, tuple)) and len(order) >= 2 and order[0] > 0:
                    L.add_edge(u, v, order=order[0])
            
            # Create R graph (right/products): nodes from ITS, edges where order[1] > 0
            R = nx.Graph()
            for node, attrs in its_graph.nodes(data=True):
                # Get right-side element type from typesGH if available
                typesGH = attrs.get('typesGH')
                if typesGH and len(typesGH) >= 2:
                    right_type = typesGH[1]  # [element, aromatic, hcount, charge, neighbors]
                    node_attrs = {
                        'element': right_type[0],
                        'aromatic': right_type[1] if len(right_type) > 1 else False,
                        'hcount': right_type[2] if len(right_type) > 2 else 0,
                        'charge': right_type[3] if len(right_type) > 3 else 0,
                    }
                else:
                    node_attrs = {
                        'element': attrs.get('element', '*'),
                        'aromatic': attrs.get('aromatic', False),
                        'hcount': attrs.get('hcount', 0),
                        'charge': attrs.get('charge', 0),
                    }
                R.add_node(node, **node_attrs)
            
            for u, v, data in its_graph.edges(data=True):
                order = data.get('order', [0, 0])
                if isinstance(order, (list, tuple)) and len(order) >= 2 and order[1] > 0:
                    R.add_edge(u, v, order=order[1])
            
            # Create K graph (context): all nodes from ITS graph
            K = nx.Graph()
            for node, attrs in its_graph.nodes(data=True):
                node_attrs = {
                    'element': attrs.get('element', '*'),
                    'aromatic': attrs.get('aromatic', False),
                    'hcount': attrs.get('hcount', 0),
                    'charge': attrs.get('charge', 0),
                }
                K.add_node(node, **node_attrs)
            
            print(f"[DEBUG]   L: {L.number_of_nodes()} nodes, {L.number_of_edges()} edges")
            print(f"[DEBUG]   R: {R.number_of_nodes()} nodes, {R.number_of_edges()} edges")
            print(f"[DEBUG]   K: {K.number_of_nodes()} nodes, {K.number_of_edges()} edges")
            
            # Save the rule
            safe_rid = rid.replace("/", "_").replace("\\", "_").replace(":", "_").replace("*", "_").replace("?", "_").replace('"', "_").replace("<", "_").replace(">", "_").replace("|", "_")
            stem = f"{safe_rid}__r0"
            print(f"[DEBUG] Saving rule as {stem}")
            
            try:
                # Find nodes that changed between left (L) and right (R)
                changed_node_ids = RuleWriting.find_changed_nodes(L, R, attributes=["charge", "hcount", "element"])
                print(f"[DEBUG]   Changed nodes: {changed_node_ids}")
                gml = RuleWriting.RulesGrammar(L, R, K, stem, changed_node_ids)
                gml_len = len(gml)
                print(f"[DEBUG]   Generated GML: {gml_len} chars")
                (RULE_DIR / f"{stem}.gml").write_text(gml, encoding="utf-8")
                print(f"[DEBUG] Saved GML for {stem}")
            except Exception as e:
                print(f"[ERROR] GML save failed for {stem}: {e}")
                import traceback
                traceback.print_exc()
            
            try:
                # Save as simple JSON representation of the graphs
                rule_json = {
                    "rule_id": stem,
                    "left": json_graph.node_link_data(L),
                    "right": json_graph.node_link_data(R), 
                    "context": json_graph.node_link_data(K)
                }
                (RULE_DIR / f"{stem}.json").write_text(json.dumps(rule_json, ensure_ascii=False), encoding="utf-8")
                print(f"[DEBUG] Saved JSON for {stem}")
            except Exception as e:
                print(f"[ERROR] JSON save failed for {stem}: {e}")
            
            # Save metadata
            try:
                meta = {
                    "reaction_id": rid,
                    "reactions": rsmi,
                    "num_nodes": K.number_of_nodes(),
                    "num_left_edges": L.number_of_edges(),
                    "num_right_edges": R.number_of_edges(),
                    "changed_nodes": changed_node_ids
                }
                (RULE_DIR / f"{stem}.meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
            except Exception:
                pass
            
            rule_index.append({"reaction_id": rid, "safe_reaction_id": safe_rid, "file_stem": stem})
            print(f"[DEBUG] Added {stem} to rule index")
        except Exception as e:
            print(f"[ERROR] Rule extraction failed for {rid}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Write index; if extraction loop did not populate rule_index, build from files
    try:
        from pathlib import Path as _Path
        stems = [p.stem for p in RULE_DIR.glob("*.gml")]
        if not rule_index and stems:
            rule_index = [{"reaction_id": s.split("__r")[0], "file_stem": s} for s in stems]
    except Exception:
        pass
    pd.DataFrame(rule_index).to_csv(DATA / "rules_index.csv", index=False)
    print(f"[OK] rules saved under {RULE_DIR}")
    print(f"[OK] wrote {DATA / 'rules_index.csv'} (rows={len(rule_index)})")

if __name__ == "__main__":
    main()
