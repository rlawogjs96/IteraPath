# 5_planner/01_mcs_template_selector.py
from pathlib import Path
import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, DataStructs

# 입력/출력 경로 (5_planner 루트 기준)
BANK   = Path("./data/template_bank.csv")
TARGET = Path("../cases/orthosporin/target.smi")
OUTJS  = Path("./data/selector_out.json")

def load_target_smiles(p: Path) -> str:
    """타깃 SMILES 로드 (첫 번째 줄)"""
    text = p.read_text(encoding="utf-8").strip()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        raise ValueError(f"Empty target file: {p}")
    return lines[0]

def morgan_fp(mol, radius=2, nbits=2048):
    """Morgan fingerprint 생성"""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)

def score_pair(target_mol, core_mol):
    """
    Tanimoto similarity + MCS atoms 계산
    
    Returns:
        (tanimoto, mcs_atoms)
    """
    # Tanimoto similarity
    fp_t = morgan_fp(target_mol)
    fp_c = morgan_fp(core_mol)
    tanim = DataStructs.TanimotoSimilarity(fp_t, fp_c)
    
    # MCS (Maximum Common Substructure)
    try:
        mcs = rdFMCS.FindMCS([target_mol, core_mol], timeout=5)
        mcs_atoms = mcs.numAtoms if mcs else 0
    except Exception:
        mcs_atoms = 0
    
    return tanim, mcs_atoms

def main():
    print("[START] MCS template selector...")
    
    # 입력 로드
    bank = pd.read_csv(BANK)
    print(f"[INFO] Loaded template bank: {len(bank)} BGCs")
    
    target_smi = load_target_smiles(TARGET)
    print(f"[INFO] Target SMILES: {target_smi}")
    
    tmol = Chem.MolFromSmiles(target_smi)
    if tmol is None:
        raise ValueError(f"Invalid target SMILES: {target_smi}")
    
    # 각 BGC 코어와 타깃 비교
    results = []
    for _, row in bank.iterrows():
        cmol = Chem.MolFromSmiles(row["core_smiles"])
        if cmol is None:
            print(f"[WARN] Invalid core SMILES for {row['bgc_id']}: {row['core_smiles']}")
            continue
        
        tanim, mcs_atoms = score_pair(tmol, cmol)
        results.append({
            "bgc_id": str(row["bgc_id"]),
            "tanimoto": round(float(tanim), 6),
            "mcs_atoms": int(mcs_atoms),
            "core_smiles": row["core_smiles"],
            "rules_seq": row["rules_seq"],
            "n_ext": int(row["n_ext"]),
            "has_kr": bool(row["has_kr"]),
            "has_dh": bool(row["has_dh"]),
            "closure_type": str(row.get("closure_type", "")),
            "at_substrate_set": str(row.get("at_substrate_set", "")),
        })
    
    if not results:
        raise RuntimeError("No valid BGC cores found in template bank")
    
    # 결정적 정렬: n_ext>=5 최우선, 그 다음 (-tanimoto, -mcs_atoms, bgc_id)
    # Orthosporin 같은 hexaketide 타깃을 위해 n_ext>=5 후보 우선
    required_n_ext = 5
    results.sort(key=lambda x: (
        0 if x["n_ext"] >= required_n_ext else 1,  # n_ext>=5가 먼저
        -x["tanimoto"],
        -x["mcs_atoms"],
        x["bgc_id"]
    ))
    
    # 출력
    OUTJS.parent.mkdir(parents=True, exist_ok=True)
    output = {
        "target": target_smi,
        "ranked": results,
        "top1": results[0] if results else None
    }
    
    with open(OUTJS, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"[DONE] Selector out written: {OUTJS}")
    if results:
        top = results[0]
        print(f"       Top-1: {top['bgc_id']} | Tanimoto={top['tanimoto']:.4f} | MCS_atoms={top['mcs_atoms']}")

if __name__ == "__main__":
    main()
