# IteraPath

Mechanistic, deterministic rule-based design toolkit for iterative PKS (iPKS) cores with gene-aware outputs.

Select a template BGC via structure similarity (MCS/Tanimoto), then replay & plan the biosynthesis using SynTemp/SynKit-style reaction rules (AAM → ITS → Rule/GML). Supports target-specific constraints (e.g., `required_n_ext=5`, `pt_mode=c2c7`, ER disabled) to guarantee carbon accounting and closure logic.

## Why this exists

BioPKS-type brute-force approaches can be non-deterministic and greedy. Here we separate **selection** from **construction**:
- **MCS** = choose the template
- **Rules** = build the pathway

This yields reproducibility, explainability, and direct mapping to gene/domain specs.

---

## Key features

- **Template selection (MCS/Tanimoto)**: pick the closest iPKS core from a bank of known BGCs.
- **Deterministic planner**: gate-based replay with hard constraints (EXT count, PT/TE presence, ER off) and soft preferences.
- **Gene-aware outputs**: automatic gene/domain spec summarizing KS/AT/KR/DH/PT/TE usage and AT substrate.
- **Accounting & validation**: stoichiometric sanity (CO₂, H₂O, thioester release), AAM→ITS consistency, audit logs.
- **Extensible constraints**: target-specific `constraints.json` (no code edits per target).

---

## iPKS configurations

- **Iteration vs termination**: repeated EXT (KS–AT; AT=malonyl) steps → PT-guided aldol closure (e.g., C2–C7) → TE release.
- **Reductive patterning**: optional KR/DH usage per cycle (ER disabled for bacterial aromatic iPKS).
- **Closure mode**: PT vs TE; ability to prefer PT=C2–C7 for isocoumarin scaffolds.
- **Domain metadata**: each rule carries domain tags and substrates → gene-level design spec.

---

## Repository layout

```
.
├─ 1_preprocessing/
│  └─ data/                # iPKS_rxn.csv → syntemp_input.csv/meta (done)
├─ 2_balancing/
│  ├─ data/                # balance_results.csv, syntemp_ready.csv (done)
│  └─ scripts/             # 00_standardize_with_syntemp.py, 01_balance_with_syntemp.py
├─ 3_templates/
│  ├─ data/                # rules_index.csv, rules/*.gml (AAM→ITS→Rule)
│  └─ scripts/             # rule extraction pipeline (done)
├─ 4_rulemeta/
│  ├─ data/                # rulemeta.csv (consolidated; add pt_mode later)
│  └─ scripts/             # 00_consolidate_rule_metadata.py
├─ 5_planner/
│  ├─ data/                # template_bank.csv, selector_out.json
│  └─ scripts/
│     ├─ 00_build_template_bank.py
│     ├─ 01_mcs_template_selector.py
│     └─ 02_deterministic_planner.py
└─ cases/
   └─ orthosporin/
      ├─ target.smi        # target SMILES
      ├─ constraints.json  # per-target constraints (recommended)
      ├─ route.json        # selected rules sequence (output)
      ├─ gene_spec.json    # domain-level spec (output)
      └─ audit.md          # gate & balance audit (output)
```

---

## ⚙️ Installation

**Python 3.9–3.11**

Required packages:
- RDKit, joblib, pandas, numpy
- SynTemp (warning: emits deprecation notices), optional SynKit (future migration)

```bash
conda create -n ipks python=3.10 -y
conda activate ipks

# RDKit install per your platform (conda-forge suggested)
conda install -c conda-forge rdkit -y

pip install joblib pandas numpy
pip install syntemp   # ok for now; SynKit migration later
# pip install synkit  # (planned)
```

---

## Quickstart (Orthosporin demo)

### 1. Consolidate rule metadata

```bash
cd 4_rulemeta/scripts
python 00_consolidate_rule_metadata.py
# writes ../data/rulemeta.csv
```

### 2. Build template bank (BGC cores)

```bash
cd ../../5_planner/scripts
python 00_build_template_bank.py
# reads ../../4_rulemeta/data/rulemeta.csv
# writes ../data/template_bank.csv
```

### 3. Put your target

Create `cases/orthosporin/target.smi`:
```
CC(CC1=CC2=CC(=CC(=C2C(=O)O1)O)O)O
```

### 4. (Recommended) Add constraints

Create `cases/orthosporin/constraints.json`:
```json
{
  "starter": "acetyl",
  "at_substrate": "mal",
  "er_allowed": false,

  "required_n_ext": 5,
  "pt_mode_allowed": ["c2c7"],
  "require_pt": true,
  "require_te": true,

  "soft": {
    "prefer_bgc_ids": [],
    "penalty_if_wrong_pt": 5.0,
    "penalty_per_ext_mismatch": 2.0
  }
}
```

### 5. Select template by MCS

```bash
python 01_mcs_template_selector.py
# reads ../data/template_bank.csv, ../../cases/orthosporin/target.smi
# writes ../data/selector_out.json
```

### 6. Deterministic planning with gates

```bash
python 02_deterministic_planner.py
# reads ../data/selector_out.json, ../../4_rulemeta/data/rulemeta.csv
# writes outputs to ../../cases/orthosporin/
```

### Outputs

Files generated in `cases/orthosporin/`:

- **`route.json`**: ordered rule IDs (EXT×N + CLOSURE), intermediates
- **`gene_spec.json`**: required domains (KS/AT/(KR)/DH/PT/TE), AT substrate
- **`audit.md`**: gate checks (EXT count, PT present, ER disabled…), formula balance notes

---

## Specifics

**Reproducibility** (deterministic sorting/tie-breaks)

**Carbon accounting** via EXT gating (e.g., exactly 5 malonyl extensions for hexaketide)

**Closure gating** (PT present; TE release; ER off)

**Gene-aware spec**: KS–AT–(KR)–DH–ACP + PT/TE, AT substrate labels

### Current limitation (transparent)

PT mode (e.g., C2–C7 vs C4–C9) is **assumed/preferred** by constraints and documentation.  
Automatic `pt_mode` inference from GML reaction centers is planned (see Roadmap).

---

## Design logic

1. **Selector vs Constructor**: use MCS to pick the template BGC, then rules to mechanistically build—not the other way around.

2. **Target-specific constraints, not global hardcoding**: all "sensible defaults" (e.g., Orthosporin → EXT×5 + PT=C2–C7) live in `constraints.json`.

3. **Explaining by genes**: every rule carries domain/substrate metadata → export a gene/domain spec.

---

## Example (Orthosporin)

**Template picked**: a hexaketide iPKS BGC (e.g., BGC1000006) because `required_n_ext=5`

**Planned route**:
```
Acetyl + 5× Malonyl (AT=mal) 
  → pentaketide/hexaketide 
    → PT(C2–C7) closure 
      → TE release
```

**Gene spec**: KS, AT, (KR optional), DH, PT, TE, ER disabled

**Tailoring**: C3 side chain oxidation/alkylation via P450/MT → handled post-PKS (optional step)

---

## Roadmap

- [ ] **PT mode inference from GML** reaction center (C2–C7 vs C4–C9) → asserted metadata, not assumed
- [ ] **Multi-objective template scoring**: balance n_ext match, Tanimoto, MCS ratio, domain fit
- [ ] **SynKit migration**: replace deprecated SynTemp modules; keep the same interfaces
- [ ] **Tailoring rule library**: common P450/MT glycosyl/O-methyl extra steps

---

## Acknowledgments

Built with [SynTemp](https://github.com/TieuLongPhan/SynTemp) for reaction rule extraction.
