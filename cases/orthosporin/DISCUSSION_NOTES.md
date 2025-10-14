# Orthosporin Design: Discussion Points

## Current Status

### ✅ Completed
- **EXT×5 enforcement**: Code-level validation (Line 176-186 in planner)
- **PT/TE presence**: Gate checks both exist (Line 188-200)
- **Carbon accounting**: Hexaketide (acetyl + 5×malonyl) guaranteed
- **Gene specification**: KS-AT-DH-KR-PT-TE, AT=mal, ER=disabled
- **Reproducibility**: Deterministic template selection (n_ext>=5 priority)

### ⚠️ Open Questions

#### 1. PT Mode Policy

**Current**: "C2-C7 assumed" (trust template)

**Options**:
- **A. Keep as-is**: Document as "preferred for isocoumarin, inferred from template"
- **B. Manual assertion**: Add `pt_mode=c2c7` to rulemeta.csv (CLOSURE rows)
- **C. Auto-detection**: Parse GML files for C2-C7 bond formation (future work)

**Trade-off**: Transparency vs development time

---

#### 2. Template Selection Strategy

**Current**: `n_ext>=5` hard filter → Tanimoto/MCS

**Result**: 
- BGC1000006 selected (n_ext=5, **Tanimoto=0.13**, MCS=15)
- BGC1000000 rejected (n_ext=4, Tanimoto=0.56, MCS=14)

**Question**: Is carbon accounting > structure similarity acceptable?

**Options**:
- **A. Keep current**: Prioritize mechanistic correctness (n_ext)
- **B. Multi-objective**: Weighted score = α·n_ext + β·Tanimoto + γ·MCS
- **C. Constraint + score**: n_ext>=5 required, then maximize Tanimoto

**Recommended weights (if B)**:
```json
{
  "n_ext_match_bonus": 10.0,
  "tanimoto_weight": 5.0,
  "mcs_ratio_weight": 2.0,
  "n_ext_error_penalty": -1.0
}
```

---

#### 3. Constraint Management

**Current**: Hard-coded in planner (`force_n_ext=5`, `DISABLE_ER=True`)

**Question**: Move to external config?

**Options**:
- **A. Keep hard-coded**: Simple, but requires code changes for new targets
- **B. Load from constraints.json**: Flexible, but adds file dependency
- **C. CLI args + config**: `--required-n-ext 5 --pt-mode c2c7`

**Recommendation**: Start with (A), migrate to (B) when handling multiple targets

---

#### 4. Metadata Completeness

**Current rulemeta.csv**:
```csv
BGC1000006_m5_to_m6, CLOSURE, domains=DH, closure=PT_or_TE, pt_mode=(missing)
```

**Options**:
- **A. Add pt_mode column**: Manual tagging for CLOSURE reactions
- **B. Enhance domain field**: Change `DH` → `DH,PT` for clarity
- **C. Both**: Most transparent

**Minimal action**: Add `pt_mode` column, fill CLOSURE rows with `c2c7` or `unknown`

---

## Decisions Needed

| Question | Options | Impact |
|----------|---------|--------|
| PT mode verification | assumed / asserted / auto | Transparency |
| Template scoring | hard filter / multi-obj / hybrid | Similarity vs mechanism |
| Constraint location | code / json / cli | Flexibility |
| Metadata enhancement | pt_mode / domains / both | Completeness |

---

## Next Steps (After Agreement)

### Immediate (< 30 min)
1. Add `pt_mode` column to rulemeta.csv (manual tagging)
2. Update planner to read constraints.json (optional)
3. Adjust gene_spec wording: "assumed" → "asserted" (if metadata added)

### Short-term (< 2 hours)
1. Implement multi-objective scoring in selector (if agreed)
2. Add CLI args for constraints (if agreed)
3. Update audit.md template with policy notes

### Medium-term (future work)
1. GML parser for automatic pt_mode detection
2. Constraints validation framework
3. Multi-target benchmark suite

---

## Files for Reference

- Current design: `cases/orthosporin/route.json`, `gene_spec.json`, `audit.md`
- Proposed constraints: `cases/orthosporin/constraints.json`
- Code locations:
  - Template selection: `5_planner/01_mcs_template_selector.py` (Line 85-93)
  - Gate enforcement: `5_planner/02_deterministic_planner.py` (Line 176-200)
  - Metadata: `4_rulemeta/data/rulemeta.csv`

---

## Recommendation

**For publication/presentation**: Current state is sufficient
- Mechanism is correct (EXT×5, PT present)
- Reproducible (deterministic gates)
- Gene-aware (domain spec)

**For long-term robustness**: Choose policy now, implement incrementally
- PT mode: Start with "asserted" (manual metadata)
- Scoring: Keep n_ext priority, add Tanimoto weight later if needed
- Constraints: External config when handling multiple targets

