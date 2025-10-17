from syntemp.SynAAM.rdt_wrapper import map_with_rdt
import importlib.resources
from pathlib import Path
import time

rdt_jar = str(importlib.resources.files("syntemp.SynAAM").joinpath("RDT_2.4.1.jar"))

# Test with different complexity levels
tests = [
    ("Simple", "CC>>CCO"),
    ("2-component", "CC.O>>CCO"),
    ("With thioester", "CC(=O)[S]>>CC(=O)O.[SH]"),
    ("Malonyl simple", "O=C(O)CC(=O)[S]>>O=C(O)CC(=O)O.[SH]"),
    ("Full EXT", "CC(=O)[S].O=C(O)CC(=O)[S]>>CC(=O)CC(=O)[S].[SH].O=C=O"),
]

for name, rxn in tests:
    working_dir = str(Path(f"C:/IteraPath/3_templates/data/rdt_test_{name.replace(' ', '_')}"))
    print(f"\n{'='*60}")
    print(f"Test: {name}")
    print(f"Reaction: {rxn}")
    try:
        result = map_with_rdt(rxn, rdt_jar, working_dir)
        if result != rxn:
            print(f"✓ MAPPED: {result[:80]}")
        else:
            print(f"✗ FAILED: Returned unmapped")
    except Exception as e:
        print(f"✗ ERROR: {e}")
    time.sleep(1)  # Avoid concurrent file access

