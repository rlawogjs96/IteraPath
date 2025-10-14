from pathlib import Path
import pandas as pd
import ast

IN_RULES = Path("../../3_templates/data/rules_index.csv")
IN_META  = Path("../../1_preprocessing/data/syntemp_input_meta.csv")
IN_BAL   = Path("../../2_balancing/data/balanced_only.csv")
OUT_PATH = Path("../data/rulemeta.csv")

def _safe_list(v):
    if pd.isna(v):
        return []
    if isinstance(v, list):
        return v
    s = str(v).strip()
    try:
        x = ast.literal_eval(s)
        if isinstance(x, list):
            return [str(i) for i in x]
    except Exception:
        pass
    return [p.strip() for p in s.split(",") if p.strip()]

def load_csv(p: Path, name: str) -> pd.DataFrame:
    if not p.exists():
        raise FileNotFoundError(f"{name} not found: {p}")
    df = pd.read_csv(p)
    # reaction_id 표준화
    if "reaction_id" in df.columns:
        df["reaction_id"] = df["reaction_id"].astype(str)
    return df

def first_nonnull(*vals):
    for v in vals:
        if v is None:
            continue
        s = str(v).strip()
        if s and s.lower() != "nan":
            return s
    return ""

def main():
    print("[START] Consolidating rule metadata...")
    print(f"[INFO] Loading rules_index from: {IN_RULES.resolve()}")
    print(f"[INFO] Loading input_meta  from: {IN_META.resolve()}")
    print(f"[INFO] Loading balanced    from: {IN_BAL.resolve()}")

    rules = load_csv(IN_RULES, "rules_index")
    meta  = load_csv(IN_META,  "input_meta")
    bal   = load_csv(IN_BAL,   "balanced_only")

    # 진단: 원본 행수/ID 세트
    rid_rules = set(rules["reaction_id"]) if "reaction_id" in rules.columns else set()
    rid_meta  = set(meta["reaction_id"]) if "reaction_id" in meta.columns else set()
    rid_bal   = set(bal["reaction_id"]) if "reaction_id" in bal.columns else set()

    print(f"[DIAG] rows: rules={len(rules)}, meta={len(meta)}, balanced={len(bal)}")
    print(f"[DIAG] unique reaction_id: rules={len(rid_rules)}, meta={len(rid_meta)}, balanced={len(rid_bal)}")

    # 기준: rules_index에 있는 reaction_id가 기준(빠지면 안 됨)
    # 병합
    base = rules.copy()

    meta_ren = meta.add_suffix("_im")
    base = base.merge(meta_ren, left_on="reaction_id", right_on="reaction_id_im", how="left")
    if "reaction_id_im" in base.columns:
        base = base.drop(columns=["reaction_id_im"])

    # balanced 컬럼 존재 여부 확인 후 병합
    bal_cols = ["reaction_id","reactant_smiles","product_smiles"]
    if "balanced" in bal.columns:
        bal_cols.append("balanced")
    base = base.merge(bal[bal_cols], on="reaction_id", how="left")

    # bgc_id 정규화: rules / meta / balanced 중 첫 유효값
    def unify_bgc(row):
        return first_nonnull(
            row.get("bgc_id"),
            row.get("bgc"),                     # 혹시 다른 이름일 경우
            row.get("bgc_id_im"),
        )

    # domains 정규화
    # 주의: rules_index.csv는 'domain' (단수), input_meta는 'domains' (복수)
    def unify_domains(row):
        d_rules = _safe_list(row.get("domain"))  # ← 'domain' (단수)
        d_meta  = _safe_list(row.get("domains_im"))
        dd = sorted(set([*d_rules, *d_meta]))
        return ",".join(dd) if dd else ""

    # substrate / at_substrate : at_substrate 우선
    def unify_substrate(row):
        at = str(row.get("at_substrate_im") or "").strip()
        if at:
            return at
        sub = str(row.get("substrate") or "").strip()
        return sub

    # step_type: meta 우선, 없으면 rule_type로 추정
    def unify_step_type(row):
        st = str(row.get("step_type_im") or "").strip()
        if st:
            return st
        rt = str(row.get("rule_type") or "").upper().strip()
        if rt in {"EXT","CLOSURE"}:
            return rt
        return "UNKNOWN"

    # closure: input_meta의 closure_im을 우선, 없으면 rules의 closure
    def unify_closure(row):
        c_im = str(row.get("closure_im") or "").strip()
        if c_im and c_im.lower() != "nan":
            return c_im
        c = str(row.get("closure") or "").strip()
        return c

    # origin: 없으면 iPKS_bacterial
    def unify_origin(row):
        o = str(row.get("origin") or "").strip()
        return o if o else "iPKS_bacterial"

    base["bgc_id"]       = base.apply(unify_bgc, axis=1)
    base["domains_norm"] = base.apply(unify_domains, axis=1)
    base["at_substrate"] = base.apply(unify_substrate, axis=1)
    base["step_type"]    = base.apply(unify_step_type, axis=1)
    base["closure"]      = base.apply(unify_closure, axis=1)
    base["origin"]       = base.apply(unify_origin, axis=1)

    # CLOSURE인데 closure가 비어 있으면 보정
    mask_clo = base["step_type"].str.upper() == "CLOSURE"
    base.loc[mask_clo & (base["closure"].eq("") | base["closure"].str.lower().eq("nan")), "closure"] = "PT_OR_TE"

    # 출력 컬럼 정리
    # module_from, module_to는 _im suffix가 붙어있을 수 있음
    out_cols = [
        "reaction_id","bgc_id",
        "module_from","module_from_im","module_to","module_to_im",
        "rule_type","step_type",
        "domains_norm","at_substrate","closure","origin",
        "reactant_smiles","product_smiles","balanced",
        "domain","substrate","domains_im","at_substrate_im","file_stem"
    ]
    # 존재하는 컬럼만 선택
    existing_cols = [c for c in out_cols if c in base.columns]
    out = base[existing_cols].copy()
    
    # module_from_im → module_from, module_to_im → module_to 변환
    if "module_from_im" in out.columns and "module_from" not in out.columns:
        out["module_from"] = out["module_from_im"]
    if "module_to_im" in out.columns and "module_to" not in out.columns:
        out["module_to"] = out["module_to_im"]
    
    # 최종 정리: _im suffix 컬럼 제거
    out = out[[c for c in out.columns if not c.endswith("_im")]].copy()

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, index=False)

    # 요약(결측 포함)
    st_counts = out["step_type"].value_counts(dropna=False).to_dict() if "step_type" in out.columns else {}
    print(f"[DONE] Wrote: {OUT_PATH.resolve()}")
    print(f"       Rows: {len(out)}")
    # bgc_id 결측도 집계
    if "bgc_id" in out.columns:
        print(f"       Unique BGCs: {out['bgc_id'].nunique(dropna=True)} (missing={int(out['bgc_id'].isna().sum())})")
    else:
        print("       Unique BGCs: NA")
    print("       Step types (incl. UNKNOWN):", st_counts)

    # 누락 진단 출력
    rid_out = set(out["reaction_id"])
    missing_from_meta = sorted(rid_rules - rid_meta)
    missing_from_bal  = sorted(rid_rules - rid_bal)
    missing_in_out    = sorted(rid_rules - rid_out)
    if missing_from_meta:
        print(f"[WARN] reaction_id present in rules_index but missing in input_meta: {missing_from_meta}")
    if missing_from_bal:
        print(f"[WARN] reaction_id present in rules_index but missing in balanced_only: {missing_from_bal}")
    if missing_in_out:
        print(f"[WARN] reaction_id present in rules_index but missing in OUTPUT: {missing_in_out}")

if __name__ == "__main__":
    main()
