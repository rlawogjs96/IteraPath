from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent / "3_templates/scripts"))
from rdt_mapper_fixed import map_with_rdt_fixed
import importlib.resources
import pandas as pd

rdt_jar = str(importlib.resources.files("syntemp.SynAAM").joinpath("RDT_2.4.1.jar"))
working_dir = "C:/IteraPath/3_templates/data/rdt_batch_debug"

# Test just 3 reactions
reactions_csv = Path("C:/IteraPath/3_templates/data/reactions.csv")
df = pd.read_csv(reactions_csv)

for idx in range(min(3, len(df))):
    rid = df.iloc[idx]["reaction_id"]
    rsmi = df.iloc[idx]["reactions"]
    
    print(f"\n{'='*60}")
    print(f"[{idx+1}] {rid}")
    print(f"Input: {rsmi[:60]}...")
    
    mapped = map_with_rdt_fixed(rsmi, rdt_jar, working_dir)
    
    print(f"Output: {mapped[:60]}...")
    print(f"Has '>>'? {'>>' in mapped}")
    print(f"Has ':'? {':' in mapped}")
    print(f"Different? {mapped != rsmi}")
    
    if mapped != rsmi and ">>" in mapped and ":" in mapped:
        print("✓ SUCCESS - Mapping worked!")
    else:
        print("✗ FAILED - Returned unmapped")

