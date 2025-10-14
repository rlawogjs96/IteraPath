import argparse
import csv
import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple

ALL_IPKS_SUBCLASS_REGEX = re.compile(r"(?i)\b(iterative\s*type\s*i|type\s*ii|type\s*iii)\b")
ITERATIVE_TYPE_I_REGEX = re.compile(r"(?i)\b(iterative\s*type\s*i)\b")

# Broad but practical set of bacterial genera commonly appearing in PKS BGCs
BACTERIAL_GENERA_REGEX = re.compile(
    r"(?i)\b("
    r"Streptomyces|Actinomadura|Actinoplanes|Amycolatopsis|Micromonospora|Nocardia|Kitasatospora|"
    r"Salinispora|Saccharopolyspora|Gordonia|Frankia|Marinispora|Verrucosispora|Streptosporangium|"
    r"Bacillus|Burkholderia|Pseudomonas|Myxococcus|Sorangium|Actinokineospora|Catellatospora|"
    r"Dactylosporangium|Lechevalieria|Rhodococcus|Dietzia|Corynebacterium|Arthrobacter|Agromyces|"
    r"Agrobacterium"
    r")\b"
)

def safe_get(d: Dict[str, Any], path: List[str], default: Any = None) -> Any:
    value: Any = d
    for key in path:
        if isinstance(value, dict) and key in value:
            value = value[key]
        else:
            return default
    return value

def iter_classes(entry: Dict[str, Any]) -> List[Dict[str, Any]]:
    classes = safe_get(entry, ["biosynthesis", "classes"], default=[])
    return classes if isinstance(classes, list) else []

def extract_compound_names(entry: Dict[str, Any]) -> List[str]:
    compounds = entry.get("compounds", [])
    names: List[str] = []
    for c in compounds:
        name = c.get("name")
        if name:
            names.append(name)
    # De-duplicate while preserving order
    seen = set()
    unique_names: List[str] = []
    for n in names:
        if n not in seen:
            seen.add(n)
            unique_names.append(n)
    return unique_names

def find_matching_pks_subclass(entry: Dict[str, Any], regex: re.Pattern) -> Optional[str]:
    for cls in iter_classes(entry):
        cls_name = cls.get("class")
        subclass = cls.get("subclass")
        if cls_name == "PKS" and isinstance(subclass, str) and regex.search(subclass):
            return subclass
    return None

def is_bacterial_taxon(entry: Dict[str, Any]) -> bool:
    org_name = safe_get(entry, ["taxonomy", "name"], default="")
    if not isinstance(org_name, str):
        return False
    return BACTERIAL_GENERA_REGEX.search(org_name) is not None

def gather_records_for_entry(entry: Dict[str, Any], file_name: str) -> Tuple[Optional[Dict[str, str]], Optional[Dict[str, str]]]:
    accession = entry.get("accession")
    organism = safe_get(entry, ["taxonomy", "name"], default="")
    compounds = extract_compound_names(entry)

    # All iPKS across Type I/II/III (Type II/III are mechanistically iterative)
    subclass_all = find_matching_pks_subclass(entry, ALL_IPKS_SUBCLASS_REGEX)
    record_all = None
    if subclass_all:
        record_all = {
            "Accession": accession or "",
            "Subclass": subclass_all,
            "Organism": organism or "",
            "Compounds": "; ".join(compounds),
            "File": file_name,
        }

    # Bacterial iterative Type I only
    subclass_iter1 = find_matching_pks_subclass(entry, ITERATIVE_TYPE_I_REGEX)
    record_bact_iter1 = None
    if subclass_iter1 and is_bacterial_taxon(entry):
        record_bact_iter1 = {
            "Accession": accession or "",
            "Subclass": subclass_iter1,
            "Organism": organism or "",
            "Compounds": "; ".join(compounds),
            "File": file_name,
        }

    return record_all, record_bact_iter1

def scan_directory(input_dir: str) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    all_ipks_rows: List[Dict[str, str]] = []
    bact_iter1_rows: List[Dict[str, str]] = []

    for entry in os.scandir(input_dir):
        if not entry.is_file():
            continue
        if not entry.name.lower().endswith(".json"):
            continue

        try:
            with open(entry.path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception:
            # Skip malformed or incompatible JSON files
            continue

        row_all, row_bact_iter1 = gather_records_for_entry(data, entry.name)
        if row_all:
            all_ipks_rows.append(row_all)
        if row_bact_iter1:
            bact_iter1_rows.append(row_bact_iter1)

    # Stable sort by accession
    all_ipks_rows.sort(key=lambda r: r.get("Accession", ""))
    bact_iter1_rows.sort(key=lambda r: r.get("Accession", ""))
    return all_ipks_rows, bact_iter1_rows

def write_csv(rows: List[Dict[str, str]], out_path: str) -> None:
    fieldnames = ["Accession", "Subclass", "Organism", "Compounds", "File"]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def main() -> None:
    parser = argparse.ArgumentParser(description="Export iPKS BGCs to CSV from MIBiG JSONs")
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing MIBiG JSON files",
    )
    parser.add_argument(
        "--out-all",
        default=None,
        help="Output CSV path for all iPKS (Iterative type I + Type II + Type III). Defaults to <input-dir>/ipks_bgcs.csv",
    )
    parser.add_argument(
        "--out-bact-iter1",
        default=None,
        help="Output CSV path for bacterial iterative Type I only. Defaults to <input-dir>/ipks_bacterial_iterative_typeI.csv",
    )

    args = parser.parse_args()
    input_dir = args.input_dir

    if not os.path.isdir(input_dir):
        raise SystemExit(f"Input directory not found: {input_dir}")

    out_all = args.out_all or os.path.join(input_dir, "ipks_bgcs.csv")
    out_bact_iter1 = args.out_bact_iter1 or os.path.join(input_dir, "ipks_bacterial_iterative_typeI.csv")

    all_ipks_rows, bact_iter1_rows = scan_directory(input_dir)
    write_csv(all_ipks_rows, out_all)
    write_csv(bact_iter1_rows, out_bact_iter1)

    print(f"Wrote {len(all_ipks_rows)} rows to: {out_all}")
    print(f"Wrote {len(bact_iter1_rows)} rows to: {out_bact_iter1}")

if __name__ == "__main__":
    main()