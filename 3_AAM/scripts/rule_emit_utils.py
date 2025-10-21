from __future__ import annotations

from pathlib import Path
import json
import networkx as nx


ORDER_TO_LABEL = {1: "-", 1.5: ":", 2: "=", 3: "#"}


def charge_to_string(charge: int) -> str:
    if charge > 0:
        return "+" if charge == 1 else f"{charge}+"
    if charge < 0:
        return "-" if charge == -1 else f"{-charge}-"
    return ""


def edge_label(edge_data: dict) -> str:
    if not isinstance(edge_data, dict):
        return "-"
    val = edge_data.get("order")
    if val is None:
        val = edge_data.get("bond_order", edge_data.get("label"))
    if isinstance(val, str) and val in ("-", ":", "=", "#"):
        return val
    if isinstance(val, (list, tuple)):
        val = val[0] if val else 1
    if isinstance(val, str):
        try:
            val = float(val)
        except Exception:
            return "-"
    if isinstance(val, (int, float)):
        if val in ORDER_TO_LABEL:
            return ORDER_TO_LABEL[val]
        if abs(val - 1.0) < 1e-6:
            return "-"
        if abs(val - 1.5) < 1e-6:
            return ":"
        if abs(val - 2.0) < 1e-6:
            return "="
        if abs(val - 3.0) < 1e-6:
            return "#"
    return "-"


def gml_section(graph: nx.Graph, section: str, changed_node_ids: list[int]) -> str:
    lines = [f"   {section} ["]
    if section == "context":
        for node, attrs in graph.nodes(data=True):
            if node in changed_node_ids:
                continue
            element = attrs.get("element", "*")
            charge = attrs.get("charge", 0)
            charge_str = charge_to_string(charge)
            lines.append(f'      node [ id {node} label "{element}{charge_str}" ]')
        for u, v, data in graph.edges(data=True):
            if u in changed_node_ids or v in changed_node_ids:
                continue
            label = edge_label(data)
            lines.append(f'      edge [ source {u} target {v} label "{label}" ]')
    else:
        for u, v, data in graph.edges(data=True):
            label = edge_label(data)
            lines.append(f'      edge [ source {u} target {v} label "{label}" ]')
        for node, attrs in graph.nodes(data=True):
            if node not in changed_node_ids:
                continue
            element = attrs.get("element", "*")
            charge = attrs.get("charge", 0)
            charge_str = charge_to_string(charge)
            lines.append(f'      node [ id {node} label "{element}{charge_str}" ]')
    lines.append("   ]")
    return "\n".join(lines)


def build_rule_gml(L: nx.Graph, R: nx.Graph, K: nx.Graph, rule_name: str, changed_node_ids: list[int]) -> str:
    parts = [
        "rule [",
        f'   ruleID "{rule_name}"',
        gml_section(L, "left", changed_node_ids),
        gml_section(K, "context", changed_node_ids),
        gml_section(R, "right", changed_node_ids),
        "]",
    ]
    return "\n".join(parts)


def sanitize_rule_stem(text: str) -> str:
    return (
        str(text)
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
        .replace("*", "_")
        .replace("?", "_")
        .replace('"', "_")
        .replace("<", "_")
        .replace(">", "_")
        .replace("|", "_")
    )


def write_rule_triplet(
    L: nx.Graph,
    R: nx.Graph,
    K: nx.Graph,
    stem: str,
    rules_dir: Path,
    changed_node_ids: list[int],
    rid: str | None = None,
    rsmi: str | None = None,
) -> None:
    rules_dir.mkdir(parents=True, exist_ok=True)
    gml = build_rule_gml(L, R, K, stem, changed_node_ids)
    (rules_dir / f"{stem}.gml").write_text(gml, encoding="utf-8")
    # JSON + meta
    from networkx.readwrite import json_graph  # local import

    rule_json = {
        "rule_id": stem,
        "left": json_graph.node_link_data(L),
        "right": json_graph.node_link_data(R),
        "context": json_graph.node_link_data(K),
    }
    (rules_dir / f"{stem}.json").write_text(json.dumps(rule_json, ensure_ascii=False), encoding="utf-8")
    if rid is not None:
        meta = {
            "reaction_id": rid,
            "reactions": rsmi,
            "num_nodes": K.number_of_nodes(),
            "num_left_edges": L.number_of_edges(),
            "num_right_edges": R.number_of_edges(),
            "changed_nodes": changed_node_ids,
        }
        (rules_dir / f"{stem}.meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

