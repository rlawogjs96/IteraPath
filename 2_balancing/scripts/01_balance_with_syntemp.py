#!/usr/bin/env python3
# 2_balancing/01_balance_with_syntemp.py
#
# Goal: Very-accurate balancing with SynTemp by augmenting true stoichiometry.
# We append small co-substrates/products so that BalanceReactionCheck sees a fully
# mass-balanced multi-molecular reaction, not "1→1" pairs.
#
# Assumptions (consistent with iterative PKS):
# - EXT (AT=malonyl):  acyl-S + malonyl-S  >> elongated acyl-S + CO2
# - KR:                ketone-S + H2       >> alcohol-S
# - DH:                beta-OH-S           >> enoyl-S + H2O   (i.e., dehydration produces water)
# - ER:                enoyl-S + H2        >> saturated-S
# - CLOSURE (TE/PT-like hydrolysis/lactonization, placeholder): 
#       acyl-S + H2O  >> cyclized/acid product + [SH]    (release of the thiol fragment)
#
# Notes:
# - We keep the dataset's thioester placeholder "[S]" (terminal sulfur atom).
# - Small molecules are encoded as: H2="[H][H]", H2O="O", CO2="O=C=O".
# - If EXT has extender other than malonyl, we leave it unhandled (fail) to stay precise.

from pathlib import Path
import pandas as pd
from syntemp.SynChemistry.balance_checker import BalanceReactionCheck

IN_STD   = Path("../data/syntemp_input.standardized.csv")   # from 00_standardize_with_syntemp.py
IN_META  = Path("../../1_preprocessing/data/syntemp_input_meta.csv")       # step_type, at_substrate, etc.
OUT_CSV  = Path("../data/balance_results.csv")
OUT_MD   = Path("../data/balance_summary.md")
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

# Small molecules (SMILES)
H2   = "[H][H]"
H2O  = "O"
CO2  = "O=C=O"

# BioPKS convention
MALONYL_THIO = "CC(=O)[S]"

def has_token(domstr: str, token: str) -> bool:
    if not isinstance(domstr, str): return False
    return token.upper() in {t.strip().upper() for t in domstr.replace(";",",").split(",") if t.strip()}

def augment_by_step(rs: str, ps: str, step: str, at_substrate: str, domains_norm: str) -> tuple[str,str,str]:
    s  = (step or "").upper()
    at = (at_substrate or "").lower()
    kr_here = has_token(domains_norm, "KR")
    dh_here = has_token(domains_norm, "DH")

    if s == "EXT":
        if at in ("mal","malonyl","malonyl-acp","malonyl_coa","malonylcoa"):
            # BioPKS: malonyl-S 축약 = CC(=O)[S]
            r_mix = f"{rs}.CC(=O)[S]"
            # KR 동시 존재 시 H2 소비
            if kr_here:
                r_mix = f"{r_mix}.[H][H]"
            # 생성물 쪽: 캐리어 방출은 [SH]로 (H까지 맞춤)
            p_mix = f"{ps}.[SH]"
            # DH 동시 존재 시 물 방출
            if dh_here:
                p_mix = f"{p_mix}.O"
            why = "EXT(mal): +acetyl-S; " + ("+H2(KR); " if kr_here else "") + ("+H2O(DH); " if dh_here else "") + "release [SH]"
            return r_mix, p_mix, why
        else:
            return rs, ps, f"EXT({at}) not handled"

    if s == "KR":
        return f"{rs}.[H][H]", ps, "KR: +H2 (reactant)"

    if s == "DH":
        return rs, f"{ps}.O", "DH: +H2O (product)"

    if s == "ER":
        return f"{rs}.[H][H]", ps, "ER: +H2 (reactant)"

    if s == "CLOSURE":
        # hydrolysis-like: 물 소비 + [SH] 방출
        return f"{rs}.O", f"{ps}.[SH]", "CLOSURE: +H2O (reactant), release [SH] (product)"

    return rs, ps, "OTHER/unknown: no augmentation"

def main():
    if not IN_STD.exists():
        raise SystemExit(f"[ERROR] Missing {IN_STD}")
    if not IN_META.exists():
        raise SystemExit(f"[ERROR] Missing {IN_META}")

    std  = pd.read_csv(IN_STD, dtype=str, keep_default_na=False)
    meta = pd.read_csv(IN_META, dtype=str, keep_default_na=False)

    need_std  = {"reaction_id", "reactant_smiles_std", "product_smiles_std"}
    need_meta = {"reaction_id", "step_type", "at_substrate"}
    miss_std  = need_std - set(std.columns)
    miss_meta = need_meta - set(meta.columns)
    if miss_std:
        raise SystemExit(f"[ERROR] Required columns missing in standardized: {miss_std}")
    if miss_meta:
        raise SystemExit(f"[ERROR] Required columns missing in meta: {miss_meta}")

    df = std.merge(
        meta[["reaction_id", "step_type", "at_substrate", "domains_norm"]],
        on="reaction_id",
        how="left"
    )

    # Build reaction strings "reactants>>products" with augmentation
    recs = []
    for _, r in df.iterrows():
        rid = r["reaction_id"]
        rs  = r["reactant_smiles_std"]
        ps  = r["product_smiles_std"]
        step= r.get("step_type", "")
        at  = r.get("at_substrate", "")

        r_mix, p_mix, why = augment_by_step(rs, ps, step, at, r.get("domains_norm",""))
        rsmi = f"{r_mix}>>{p_mix}"
        recs.append({"reaction_id": rid, "reactions": rsmi, "augmentation": why, "step_type": step, "at_substrate": at})

    # Use the official checker properly (batch)
    checker = BalanceReactionCheck(n_jobs=4, verbose=0)
    balanced, unbalanced = checker.dicts_balance_check(recs, rsmi_column="reactions")

    # Merge results back to a single table
    out_rows = []
    for row in balanced:
        out_rows.append({
            "reaction_id": row["reaction_id"],
            "reactions": row["reactions"],
            "balanced": True,
            "augmentation": row["augmentation"],
            "step_type": row.get("step_type",""),
            "at_substrate": row.get("at_substrate",""),
        })
    for row in unbalanced:
        out_rows.append({
            "reaction_id": row["reaction_id"],
            "reactions": row["reactions"],
            "balanced": False,
            "augmentation": row["augmentation"],
            "step_type": row.get("step_type",""),
            "at_substrate": row.get("at_substrate",""),
        })

    rep = pd.DataFrame(out_rows)
    rep = rep.sort_values(["balanced","reaction_id"], ascending=[False, True]).reset_index(drop=True)
    rep.to_csv(OUT_CSV, index=False)

    total = len(rep)
    n_ok  = int(rep["balanced"].sum())
    n_bad = total - n_ok
    OUT_MD.write_text(
        "\n".join([
            "# SynTemp Balance Summary (official checker, stoichiometry-augmented)",
            f"- Total: {total}",
            f"- Balanced: {n_ok}",
            f"- Unbalanced: {n_bad}",
        ]),
        encoding="utf-8"
    )
    print(f"[OK] wrote {OUT_CSV} and {OUT_MD}")

if __name__ == "__main__":
    main()
