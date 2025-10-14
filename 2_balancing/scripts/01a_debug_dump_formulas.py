#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from syntemp.SynChemistry.balance_checker import BalanceReactionCheck
from syntemp.utils.chemutils import get_combined_molecular_formula

IN_STD   = Path("../data/syntemp_input.standardized.csv")   # from 00_standardize_with_syntemp.py
IN_META  = Path("../../1_preprocessing/data/syntemp_input_meta.csv")       # step_type, at_substrate, etc.
OUT_CSV  = Path("../data/balance_debug_formulas.csv")

# BioPKS 축약 컨벤션: malonyl-S는 acetyl-S placeholder로 공급
MALONYL_THIO = "CC(=O)[S]"
H2   = "[H][H]"
H2O  = "O"
# Thiol released in CLOSURE/EXT is kept as [S] per BioPKS notation.
# Can change to [SH] if needed to verify H balance.

def augment_by_step(rs: str, ps: str, step: str, at_substrate: str):
    s  = (step or "").upper()
    at = (at_substrate or "").lower()

    if s == "EXT":
        if at in ("mal","malonyl","malonyl-acp","malonyl_coa","malonylcoa"):
            # BioPKS convention assumption: +acetyl-S, release carrier thiol
            r_mix = f"{rs}.{MALONYL_THIO}"
            p_mix = f"{ps}.[S]"           # ← Can change to "[SH]" to check H impact
            why   = "EXT(mal): +acetyl-S (BioPKS), +[S] on product"
            return r_mix, p_mix, why
        else:
            return rs, ps, f"EXT({at}) not handled"

    if s == "KR":
        return f"{rs}.{H2}", ps, "KR: +H2 (reactant)"

    if s == "DH":
        return rs, f"{ps}.{H2O}", "DH: +H2O (product)"

    if s == "ER":
        return f"{rs}.{H2}", ps, "ER: +H2 (reactant)"

    if s == "CLOSURE":
        return rs, f"{ps}.[S]", "CLOSURE: +[S] (thiol release)"

    return rs, ps, "OTHER/unknown"

def split_lr(rsmi: str):
    L, R = rsmi.split(">>")
    return L, R

def to_dict(formula: str):
    # "C4H6O2S2" → {"C":4,"H":6,"O":2,"S":2}
    out = {}
    num = ""
    sym = ""
    for ch in formula:
        if ch.isalpha():
            if sym:
                out[sym] = out.get(sym, 0) + (int(num) if num else 1)
            sym = ch
            num = ""
        else:
            num += ch
    if sym:
        out[sym] = out.get(sym, 0) + (int(num) if num else 1)
    return out

def dict_diff(a: dict, b: dict):
    keys = sorted(set(a)|set(b))
    return {k: b.get(k,0)-a.get(k,0) for k in keys if b.get(k,0)!=a.get(k,0)}

def main():
    if not (IN_STD.exists() and IN_META.exists()):
        raise SystemExit("Missing inputs. Run 00_standardize first.")

    std  = pd.read_csv(IN_STD, dtype=str, keep_default_na=False)
    meta = pd.read_csv(IN_META, dtype=str, keep_default_na=False)
    df   = std.merge(meta[["reaction_id","step_type","at_substrate"]], on="reaction_id", how="left")

    rows = []
    for _, r in df.iterrows():
        rid  = r["reaction_id"]
        rs   = r["reactant_smiles_std"]
        ps   = r["product_smiles_std"]
        step = r.get("step_type","")
        at   = r.get("at_substrate","")

        r_mix, p_mix, why = augment_by_step(rs, ps, step, at)
        rsmi = f"{r_mix}>>{p_mix}"

        # Official function only returns true/false, so we calculate formulas ourselves for recording
        L, R = split_lr(rsmi)
        fL = get_combined_molecular_formula(L)   # e.g.: "C4H6O2S2"
        fR = get_combined_molecular_formula(R)

        dL = to_dict(fL); dR = to_dict(fR)
        diff = dict_diff(dL, dR)

        ok = BalanceReactionCheck.rsmi_balance_check(rsmi)

        rows.append({
            "reaction_id": rid,
            "step_type": step,
            "at_substrate": at,
            "reactions": rsmi,
            "f_left": fL,
            "f_right": fR,
            "diff": diff,          # Which elements differ and by how much
            "balanced": bool(ok),
            "augmentation": why
        })

    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[OK] wrote {OUT_CSV} (rows={len(out)})")
    # Useful tip: check the first 5 diffs.

if __name__ == "__main__":
    main()
