#!/usr/bin/env python3
import json

with open('gene_spec.json') as f:
    spec = json.load(f)

print("=" * 60)
print("ORTHOSPORIN BIOSYNTHETIC GENE CLUSTER DESIGN")
print("=" * 60)
print()

print("### iPKS Core Module (Iterative, 5 cycles)")
print("  Required domains:", ", ".join(spec['required_domains']))
print("  AT substrate:", spec['AT_substrates'][0], "(malonyl-CoA)")
print("  Closure mechanism:", spec['closure'])
print("  ER:", "disabled" if not spec['ER_allowed'] else "enabled")
print("  KR mode:", spec['KR_mode'], "(cycles 2, 4)")
print()

print("### Tailoring Enzymes (Post-PKS)")
for i, t in enumerate(spec['tailoring'], 1):
    print(f"  Step {i}: {t['name']}")
    print(f"          Domain: {t['domain']}")
    if 'substrate' in t:
        print(f"          Substrate: {t['substrate']}")
    print(f"          Change: {t['delta']}")
    print()

print("=" * 60)
print("RESULT: iPKS core (C11H8O3) + Tailoring (+C2H4O)")
print("        → Orthosporin (C13H12O5)")
print("=" * 60)

