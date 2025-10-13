#!/usr/bin/env python3
"""
Validate iPKS reaction pairs (minimal QC) and add lightweight tags.

Input : 1_preprocessing/iPKS_rxn.reactions.csv  (module_from/to, reactant_smiles, product_smiles, domains, at_substrate, ...)
Output: 1_preprocessing/iPKS_rxn.reactions.validated.csv
        1_preprocessing/iPKS_rxn.reactions_qc.md

What it checks:
  - basic field presence
  - domain tokens normalization (drops ACP)
  - thioester presence change (reactant/product contains "[S]"?)
  - step-type heuristics:
      * EXT step: domains has KS & AT (malonyl?) and both reactant/product have "[S]"
      * KR step: domains has KR  (optional)
      * DH step: domains has DH  (optional)
      * CLOSURE candidate: product has NO "[S]"   (PT/TE class; not deciding C2–C7 vs TE yet)
Notes:
  - RDKit 없는 환경도 동작 (원소/질량 밸런스는 이후 단계에서 상세 검증)
"""

import csv
import re
from pathlib import Path
import pandas as pd

INCSV  = Path("../data/iPKS_rxn.reactions.csv")
OUTCSV = Path("../data/iPKS_rxn.reactions.validated.csv")
OUTQC  = Path("../data/iPKS_rxn.reactions_qc.md")

def norm_domains(cell: str):
    if not isinstance(cell, str) or not cell.strip():
        return ""
    toks = re.split(r"[,\s;/]+", cell.strip())
    out = []
    for t in toks:
        u = t.upper()
        if u in {"", "ACP"}:  # ACP는 제외
            continue
        # 간단 별칭
        if u in {"THIOESTERASE"}: u = "TE"
        if u in {"PRODUCTTEMPLATE", "PTDOMAIN", "PTD"}: u = "PT"
        if u in {"KETOREDUCTASE"}: u = "KR"
        if u in {"DEHYDRATASE"}: u = "DH"
        if u in {"ENRED", "ENOYLREDUCTASE"}: u = "ER"
        if u in {"BETA-KETOACYLSYNTHASE", "KSQ"}: u = "KS"
        if u in {"ACYLTRANSFERASE"}: u = "AT"
        out.append(u)
    # dedup + 안정적 순서
    seen, uniq = set(), []
    for u in out:
        if u and u not in seen:
            uniq.append(u); seen.add(u)
    return ",".join(uniq)

def has_thioester(smiles: str) -> bool:
    # 매우 단순 휴리스틱: "[S]" 토큰 유무로 판단 (ACP/CoA 표기 관례를 그대로 활용)
    return isinstance(smiles, str) and "[S]" in smiles

def main():
    if not INCSV.exists():
        raise SystemExit(f"[ERROR] Missing {INCSV}")

    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)

    # 정규화
    for c in ["reactant_smiles","product_smiles","domains","at_substrate","kr_annotation","dh_annotation","er_annotation","te_annotation"]:
        if c not in df.columns: df[c] = ""
        df[c] = df[c].astype(str).str.strip()

    # 도메인 정리 (ACP 제거)
    df["domains_norm"] = df["domains"].map(norm_domains)

    # 티오에스터 플래그
    df["reactant_has_thio"] = df["reactant_smiles"].map(has_thioester)
    df["product_has_thio"]  = df["product_smiles"].map(has_thioester)

    # 스텝 타입 휴리스틱
    def infer_step(row):
        doms = set([d for d in row["domains_norm"].split(",") if d])
        r_thio = row["reactant_has_thio"]; p_thio = row["product_has_thio"]
        at = (("KS" in doms) and ("AT" in doms))
        kr = ("KR" in doms)
        dh = ("DH" in doms)
        # 규칙
        if not p_thio:         # 제품에 thioester가 없으면 방출/폐환 계열
            return "CLOSURE"
        if at and r_thio and p_thio:
            # 연장(EXT) 단계로 해석
            return "EXT"
        if kr and r_thio and p_thio:
            return "KR"
        if dh and r_thio and p_thio:
            return "DH"
        # 기타는 보류
        return "OTHER"

    df["step_type"] = df.apply(infer_step, axis=1)

    # 간단 QC 집계
    total = len(df)
    n_ext = (df["step_type"]=="EXT").sum()
    n_clo = (df["step_type"]=="CLOSURE").sum()
    n_kr  = (df["step_type"]=="KR").sum()
    n_dh  = (df["step_type"]=="DH").sum()
    n_other = (df["step_type"]=="OTHER").sum()

    # 도메인-스텝 상충 검출 (예: domains가 DH인데 CLOSURE로 나오는 경우 등)
    conflicts = []
    for i, r in df.iterrows():
        doms = set([d for d in r["domains_norm"].split(",") if d])
        st = r["step_type"]
        if st=="EXT" and not ({"KS","AT"} & doms):
            conflicts.append((i, "EXT but missing KS/AT in domains"))
        if st=="KR" and "KR" not in doms:
            conflicts.append((i, "KR but KR not in domains"))
        if st=="DH" and "DH" not in doms:
            conflicts.append((i, "DH but DH not in domains"))
        if st=="CLOSURE" and (("KS" in doms) or ("AT" in doms)):
            # 폐환인데 KS/AT만 있는 케이스는 의심: PT/TE 누락 가능성
            conflicts.append((i, "CLOSURE but domains lack PT/TE flag (suspect PT/TE)"))

    # 출력 저장
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTCSV, index=False, quoting=csv.QUOTE_MINIMAL)

    # QC 리포트
    lines = []
    lines.append("# Reactions QC Report\n")
    lines.append(f"- Total: {total}")
    lines.append(f"- EXT: {n_ext}, KR: {n_kr}, DH: {n_dh}, CLOSURE: {n_clo}, OTHER: {n_other}\n")
    if conflicts:
        lines.append("## Potential conflicts")
        for idx, msg in conflicts:
            row = df.iloc[idx]
            lines.append(f"- row {idx} | {row.get('bgc_id','')} {row.get('module_from','')}→{row.get('module_to','')}: {msg}")
    else:
        lines.append("## Potential conflicts\n- None")

    OUTQC.write_text("\n".join(lines), encoding="utf-8")
    print(f"[OK] wrote {OUTCSV}")
    print(f"[OK] wrote {OUTQC}")

if __name__ == "__main__":
    main()
