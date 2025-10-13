#!/usr/bin/env python3
"""
Add minimal metadata:
- reaction_id (stable)
- closure flag (PT_or_TE placeholder for CLOSURE step)
- domains_norm_no_acp (safety)
Input : 1_preprocessing/iPKS_rxn.reactions.validated.csv
Output: 1_preprocessing/iPKS_rxn.reactions.labeled.csv
"""

import re
from pathlib import Path
import pandas as pd

INCSV  = Path("../data/iPKS_rxn.reactions.validated.csv")
OUTCSV = Path("../data/iPKS_rxn.reactions.labeled.csv")

def strip_acp(s: str) -> str:
    s = (s or "").strip()
    toks = re.split(r"[,\s;/]+", s) if s else []
    toks = [t for t in toks if t and t.upper() != "ACP"]
    # dedup while preserving order
    seen, out = set(), []
    for t in toks:
        if t not in seen:
            out.append(t); seen.add(t)
    return ",".join(out)

def main():
    if not INCSV.exists():
        raise SystemExit(f"[ERROR] Missing {INCSV}")

    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)

    # reaction_id 생성
    required = {"bgc_id", "module_from", "module_to"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise SystemExit(f"[ERROR] Required columns missing: {missing}")

    df["reaction_id"] = df.apply(
        lambda r: f"{r['bgc_id']}_m{r['module_from']}_to_m{r['module_to']}", axis=1
    )

    # domains_norm 없으면 대체
    if "domains_norm" not in df.columns:
        if "domains" in df.columns:
            df["domains_norm"] = df["domains"].astype(str)
        else:
            df["domains_norm"] = ""

    # closure 컬럼 안전하게 보장 (여기서 .astype 사용은 Series에만)
    if "closure" not in df.columns:
        df["closure"] = ""
    else:
        df["closure"] = df["closure"].astype(str)

    # CLOSURE 단계에 placeholder 라벨 부여
    if "step_type" in df.columns:
        mask = (df["step_type"] == "CLOSURE") & (df["closure"].str.len() == 0)
        df.loc[mask, "closure"] = "PT_or_TE"

    # ACP 제거한 도메인
    df["domains_norm_no_acp"] = df["domains_norm"].map(strip_acp)

    df.to_csv(OUTCSV, index=False)
    print(f"[OK] wrote {OUTCSV} (rows={len(df)})")

if __name__ == "__main__":
    main()