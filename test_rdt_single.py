from syntemp.SynAAM.rdt_wrapper import map_with_rdt
import importlib.resources
from pathlib import Path

# Get RDT JAR
rdt_jar = str(importlib.resources.files("syntemp.SynAAM").joinpath("RDT_2.4.1.jar"))
working_dir = str(Path("C:/IteraPath/3_templates/data/rdt_test"))

# Test first reaction with CORRECT malonyl-CoA
test_rxn = "CC(=O)[S].O=C(O)CC(=O)[S]>>CC(=O)CC(=O)[S].[SH].O=C=O"

print(f"Testing: {test_rxn}")
print(f"RDT JAR: {rdt_jar}")
print(f"Working dir: {working_dir}")
print("")

try:
    result = map_with_rdt(test_rxn, rdt_jar, working_dir)
    print(f"Result: {result}")
    print(f"Same as input? {result == test_rxn}")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

