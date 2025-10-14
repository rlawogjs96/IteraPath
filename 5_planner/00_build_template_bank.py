# 5_planner/00_build_template_bank.py
from pathlib import Path
import pandas as pd

# 입력 경로 (5_planner 루트 기준)
BALANCED = Path("../2_balancing/data/balanced_only.csv")
RULEMETA = Path("../4_rulemeta/data/rulemeta.csv")
OUTBANK  = Path("./data/template_bank.csv")

def load_tables():
    """입력 파일 로드"""
    bal = pd.read_csv(BALANCED)
    meta = pd.read_csv(RULEMETA)
    # reaction_id 문자열화
    for df in (bal, meta):
        if "reaction_id" in df.columns:
            df["reaction_id"] = df["reaction_id"].astype(str)
    return bal, meta

def build_bank(bal: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    BGC별 템플릿 뱅크 생성
    
    각 BGC에 대해:
    - 최종 코어 SMILES (CLOSURE 단계의 product)
    - 룰 시퀀스 (reaction_id 리스트)
    - 메타 요약 (n_ext, has_kr, has_dh, closure_type, at_substrate_set)
    """
    # balanced_only를 베이스로 시작
    ready = bal.copy()
    
    # rulemeta에서 step_type, closure 정보 병합
    if not meta.empty and "reaction_id" in meta.columns:
        meta_cols = ["reaction_id", "step_type", "closure"]
        available_cols = [c for c in meta_cols if c in meta.columns]
        ready = ready.merge(
            meta[available_cols].drop_duplicates("reaction_id"),
            on="reaction_id",
            how="left",
            suffixes=("_bal", "_meta")
        )
        
        # step_type 통합: _meta 우선 (_bal은 balanced_only에 없을 수 있음)
        if "step_type_meta" in ready.columns:
            ready["step_type"] = ready["step_type_meta"]
        elif "step_type_bal" in ready.columns:
            ready["step_type"] = ready["step_type_bal"]
    
    # step_type이 없으면 추론
    if "step_type" not in ready.columns:
        ready["step_type"] = "UNKNOWN"
    
    # 최종 코어 추출 함수
    def final_core(dfbgc: pd.DataFrame) -> str:
        # CLOSURE 단계의 product_smiles 우선
        if "step_type" in dfbgc.columns:
            clo = dfbgc[dfbgc["step_type"] == "CLOSURE"]
            if not clo.empty and "product_smiles" in clo.columns:
                sorted_clo = clo.sort_values(["module_from", "module_to"], na_position='last')
                return sorted_clo.iloc[-1]["product_smiles"]
        # fallback: BGC 내 마지막 product_smiles
        if "module_from" in dfbgc.columns and "module_to" in dfbgc.columns:
            sorted_bgc = dfbgc.sort_values(["module_from", "module_to"], na_position='last')
        else:
            sorted_bgc = dfbgc
        return sorted_bgc.iloc[-1]["product_smiles"]
    
    rows = []
    for bgc, dfb in ready.groupby("bgc_id"):
        dfb = dfb.sort_values(["module_from", "module_to"], na_position='last').reset_index(drop=True)
        
        # 최종 코어 SMILES
        core = final_core(dfb)
        
        # 룰 시퀀스
        rules_seq = dfb["reaction_id"].astype(str).tolist()
        
        # AT substrate set
        at_set = []
        if "at_substrate" in dfb.columns:
            at_set = sorted(set(dfb["at_substrate"].dropna().astype(str).tolist()))
            at_set = [x for x in at_set if x and x.lower() not in ['nan', 'none', '']]
        
        # EXT 개수 (디버깅 강화)
        n_ext = 0
        if "step_type" in dfb.columns:
            ext_mask = (dfb["step_type"] == "EXT")
            n_ext = int(ext_mask.sum())
            # 디버깅: step_type이 모두 NaN/UNKNOWN인 경우 대비
            if n_ext == 0:
                # rules_seq 길이에서 CLOSURE 1개 빼기로 추정
                total_steps = len(dfb)
                closure_count = int((dfb["step_type"] == "CLOSURE").sum())
                if closure_count > 0:
                    n_ext = total_steps - closure_count
                else:
                    # 최후의 fallback: 전체 스텝 - 1
                    n_ext = max(0, total_steps - 1)
        
        # KR/DH 존재 여부
        has_kr = False
        has_dh = False
        if "domains" in dfb.columns:
            has_kr = bool(dfb["domains"].fillna("").astype(str).str.contains("KR", case=False).any())
            has_dh = bool(dfb["domains"].fillna("").astype(str).str.contains("DH", case=False).any())
        
        # Closure type
        closure_type = "PT_OR_TE"
        if "closure" in dfb.columns:
            closure_vals = dfb["closure"].dropna().astype(str).unique().tolist()
            closure_vals = [x for x in closure_vals if x and x.lower() not in ['nan', 'none', '']]
            if closure_vals:
                closure_type = ";".join(sorted(closure_vals))
        
        rows.append({
            "bgc_id": bgc,
            "core_smiles": core,
            "rules_seq": "|".join(rules_seq),
            "n_ext": n_ext,
            "has_kr": has_kr,
            "has_dh": has_dh,
            "closure_type": closure_type,
            "at_substrate_set": "|".join(at_set) if at_set else "",
        })
    
    return pd.DataFrame(rows)

def main():
    print("[START] Building template bank...")
    OUTBANK.parent.mkdir(parents=True, exist_ok=True)
    
    bal, meta = load_tables()
    print(f"[INFO] Loaded: balanced_only={len(bal)} rows, rulemeta={len(meta)} rows")
    
    bank = build_bank(bal, meta)
    bank.to_csv(OUTBANK, index=False)
    
    print(f"[DONE] Template bank written: {OUTBANK}")
    print(f"       Rows: {len(bank)}")
    print(f"       BGCs: {', '.join(bank['bgc_id'].astype(str).tolist())}")

if __name__ == "__main__":
    main()
