#!/usr/bin/env python3
# 3_templates/scripts/01_aam_its_rules_complete.py
# AAM → ITS → H-adjustment → Rules 전체 파이프라인 (한 번에)

from pathlib import Path
import json
import pandas as pd
import warnings

# Suppress FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning, module="syntemp")

# SynAAM
from syntemp.SynAAM.atom_map_consensus import AAMConsensus
from syntemp.SynAAM.aam_postprocess import AMMPostprocessor
# SynITS
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynITS.its_hadjuster import ITSHAdjuster
# SynRule
from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynRule.rule_writing import RuleWriting

BASE = Path(__file__).resolve().parent
DATA = BASE.parent / "data"
META_PATH = BASE.parent.parent / "1_preprocessing" / "data" / "syntemp_input_meta.csv"

INCSV = DATA / "reactions.csv"
OUT_AAM = DATA / "aam.ndjson"
RULE_DIR = DATA / "rules"
RULE_DIR.mkdir(parents=True, exist_ok=True)

# 매퍼 설정 (의존성 문제로 RXNMapper만 사용)
MAPPERS = ["rxn_mapper"]  # Graphormer(chytorch 필요), LocalMapper(DGL C++ 필요) 제외

def infer_rule_metadata(reaction_id: str, meta_df: pd.DataFrame) -> dict:
    """유전자→룰 매핑용 메타데이터 추출"""
    meta = {"domain": [], "substrate": "unknown", "stereochem_tag": "unknown", "closure": "", "origin": ""}
    
    if meta_df is None or reaction_id not in meta_df["reaction_id"].values:
        return meta
    
    row = meta_df[meta_df["reaction_id"] == reaction_id].iloc[0]
    step_type = str(row.get("step_type", "")).upper()
    at_sub = str(row.get("at_substrate", "")).lower()
    closure = str(row.get("closure", ""))
    domains_norm = str(row.get("domains_norm", ""))
    
    # Domain 파싱
    if domains_norm:
        meta["domain"] = [d.strip() for d in domains_norm.split(",") if d.strip()]
    
    # Substrate (AT 특이성)
    if step_type == "EXT":
        if "mal" in at_sub:
            meta["substrate"] = "malonyl"
        elif "acet" in at_sub or at_sub == "ac":
            meta["substrate"] = "acetyl"
        else:
            meta["substrate"] = at_sub
    
    # Closure pattern
    if step_type == "CLOSURE":
        if "pt" in closure.lower():
            meta["closure"] = closure.upper()
        elif "te" in closure.lower():
            meta["closure"] = "TE_LAC"
        else:
            meta["closure"] = "PT_or_TE"
    
    # Stereochemistry
    if "KR" in meta["domain"]:
        meta["stereochem_tag"] = "KR_unknown"
    
    # Origin
    meta["origin"] = "iPKS_bacterial"
    
    return meta

def main():
    if not INCSV.exists():
        raise SystemExit("Run 00_make_reaction_smiles.py first.")
    
    # 메타데이터 로드
    meta_df = None
    if META_PATH.exists():
        meta_df = pd.read_csv(META_PATH, dtype=str, keep_default_na=False)
        print(f"[INFO] Loaded metadata for {len(meta_df)} reactions")
    
    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)
    
    # 1) AAM 앙상블
    batch = [{"reaction_id": rid, "reactions": rsmi} for rid, rsmi in df.values]
    print(f"[INFO] Running AAM with {len(batch)} reactions using: {MAPPERS}")
    
    cons = AAMConsensus(data=batch, mappers=MAPPERS)
    batch = cons.batch_consensus(batch, rsmi_column="reactions")
    print(f"[INFO] AAM consensus completed")
    
    # 2) AAM 후처리
    post = AMMPostprocessor()
    batch = post.parallel_postprocess(batch, MAPPERS, threshold=len(MAPPERS), n_jobs=4, verbose=0)
    print(f"[INFO] AAM postprocessing completed")
    
    # 3) ITS 추출
    correct, incorrect = ITSExtraction.parallel_process_smiles(
        batch,
        MAPPERS,
        n_jobs=4,
        verbose=0,
        id_column="reaction_id",
        check_method="RC",
        confident_mapper="rxn_mapper",
        sanitize=True
    )
    
    print(f"[INFO] ITS: {len(correct)} correct, {len(incorrect)} incorrect")
    
    # 4) H-보정 (correct만)
    if correct:
        hadj = ITSHAdjuster()
        try:
            fixed = hadj.process_graph_data_parallel(
                correct,
                column="ITSGraph",
                n_jobs=4,
                verbose=0,
                max_hydrogen=4,
                get_priority_graph=True
            )
            print(f"[INFO] H-adjustment completed for {len(fixed)} reactions")
        except Exception as e:
            print(f"[WARN] H-adjustment failed: {e}, using unfixed ITS")
            fixed = correct
    else:
        print("[WARN] No correct ITS found, cannot extract rules")
        fixed = []
    
    # 5) 룰 추출
    rule_index = []
    aam_records = []
    
    for rec in fixed:
        rid = rec.get("reaction_id", rec.get("R-id", "unknown"))
        
        # AAM 기록
        selected_mapper = None
        chosen_aam = None
        for mapper in MAPPERS:
            if mapper in rec and rec[mapper] and ">>" in str(rec[mapper]):
                chosen_aam = rec[mapper]
                selected_mapper = mapper
                break
        
        aam_records.append({
            "reaction_id": rid,
            "selected_mapper": selected_mapper or "none",
            "aam": chosen_aam,
            "equivariant": rec.get("equivariant", 0)
        })
        
        # ITSGraph 추출
        its_tuple = rec.get("ITSGraph")
        if its_tuple is None or not isinstance(its_tuple, tuple) or len(its_tuple) != 3:
            continue
        
        L, R, K = its_tuple  # Left (reactants), Right (products), K (context/ITS)
        
        # 메타데이터
        rule_meta = infer_rule_metadata(rid, meta_df)
        
        try:
            # 룰 추출 (올바른 API)
            rule_tuple = RuleExtraction.extract_reaction_rules(
                L, R, K,
                extend=True,
                n_knn=2  # iPKS용 반경 확장
            )
            
            if rule_tuple and len(rule_tuple) == 3:
                L_rule, R_rule, K_rule = rule_tuple
                stem = f"{rid}__r0"
                
                # 변경된 노드 찾기
                changed_nodes = RuleWriting.find_changed_nodes(L_rule, R_rule, attributes=["charge"])
                
                # GML 저장 (올바른 API)
                try:
                    gml = RuleWriting.RulesGrammar(L_rule, R_rule, K_rule, rid, changed_nodes)
                    (RULE_DIR / f"{stem}.gml").write_text(gml, encoding="utf-8")
                except Exception as e:
                    print(f"[WARN] GML save failed for {rid}: {e}")
                
                # 메타데이터만 JSON으로 저장 (SynTemp는 JSON rule 미지원)
                rule_meta_full = {
                    "reaction_id": rid,
                    "domain": rule_meta["domain"],
                    "substrate": rule_meta["substrate"],
                    "stereochem_tag": rule_meta["stereochem_tag"],
                    "closure": rule_meta["closure"],
                    "origin": rule_meta["origin"],
                    "rule_type": "EXT" if "AT" in rule_meta["domain"] else
                                 "CLOSURE" if rule_meta["closure"] else "RED",
                    "rule_file": f"{stem}.gml"
                }
                
                (RULE_DIR / f"{stem}.meta.json").write_text(
                    json.dumps(rule_meta_full, ensure_ascii=False, indent=2),
                    encoding="utf-8"
                )
                
                rule_index.append({
                    "reaction_id": rid,
                    "file_stem": stem,
                    **rule_meta
                })
                
        except Exception as e:
            print(f"[WARN] Rule extraction failed for {rid}: {e}")
    
    # 6) 저장
    with OUT_AAM.open("w", encoding="utf-8") as f:
        for rec in aam_records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    
    idx_df = pd.DataFrame(rule_index)
    idx_df.to_csv(DATA / "rules_index.csv", index=False)
    
    print(f"[OK] wrote {OUT_AAM}")
    print(f"[OK] wrote {DATA / 'rules_index.csv'} ({len(rule_index)} rules)")
    print(f"[OK] rules saved under {RULE_DIR}")
    print(f"[INFO] Metadata: domain, substrate, stereochem_tag, closure, origin")

if __name__ == "__main__":
    main()

