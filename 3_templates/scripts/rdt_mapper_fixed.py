#!/usr/bin/env python3
"""
Fixed RDT mapper - removes shell redirection that causes Windows file locking
"""
import os
import shutil
from uuid import uuid4
import pandas as pd
import json
from pathlib import Path

def map_with_rdt_fixed(reaction_smiles: str, rdt_jar_path: str, working_dir: str) -> str:
    """Fixed version without shell redirection"""
    unique_id = uuid4()
    unique_dir = os.path.join(working_dir, f"RDT_{unique_id}")
    output_file = "ECBLAST_smiles_AAM.txt"
    
    original_cwd = os.getcwd()
    
    try:
        os.makedirs(unique_dir, exist_ok=True)
        os.chdir(unique_dir)
        
        # Run RDT WITHOUT shell redirection (RDT writes the file itself)
        # Use full Java path to avoid PATH issues in conda environments
        java_exe = r"C:\Program Files\Microsoft\jdk-17.0.16.8-hotspot\bin\java.exe"
        command = f'"{java_exe}" -jar "{rdt_jar_path}" -Q SMI -q "{reaction_smiles}" -c -j AAM -f TEXT'
        exit_code = os.system(command)
        
        # Check if RDT succeeded
        if not os.path.exists(output_file):
            print(f"[WARN] RDT did not create output file for: {reaction_smiles[:50]}...")
            return reaction_smiles
        
        # Read the mapped reaction (line 4, index 3)
        with open(output_file, "r") as f:
            lines = f.read().splitlines()
            if len(lines) < 4:
                print(f"[WARN] RDT output too short for: {reaction_smiles[:50]}...")
                return reaction_smiles
            mapped = lines[3]
        
        # Validate mapping
        if ">>" in mapped and ":" in mapped:
            return mapped
        else:
            print(f"[WARN] RDT mapping invalid for: {reaction_smiles[:50]}...")
            return reaction_smiles
            
    except Exception as e:
        print(f"[ERROR] RDT failed: {e}")
        return reaction_smiles
    finally:
        os.chdir(original_cwd)
        # Clean up temp directory
        if os.path.exists(unique_dir):
            shutil.rmtree(unique_dir, ignore_errors=True)


def main():
    import importlib.resources
    
    # Setup
    BASE = Path(__file__).resolve().parent
    DATA = BASE.parent / "data"
    INCSV = DATA / "reactions.csv"
    OUT_DIR = DATA / "individual_mappers"
    
    rdt_jar = str(importlib.resources.files("syntemp.SynAAM").joinpath("RDT_2.4.1.jar"))
    working_dir = str(DATA / "rdt_fixed_temp")
    
    # Read reactions
    print(f"[INFO] Reading reactions from {INCSV}")
    df = pd.read_csv(INCSV, dtype=str, keep_default_na=False)
    df = df.dropna(subset=["reaction_id","reactions"]).copy()
    df = df[df["reactions"].str.contains(">>")].drop_duplicates("reaction_id")
    
    print(f"[INFO] Processing {len(df)} reactions with fixed RDT wrapper")
    
    # Map reactions
    results = []
    for idx, (rid, rsmi) in enumerate(df.values, 1):
        print(f"[{idx}/{len(df)}] {rid}...", end=" ")
        mapped = map_with_rdt_fixed(rsmi, rdt_jar, working_dir)
        
        success = mapped != rsmi and ">>" in mapped and ":" in mapped
        results.append({
            "reaction_id": rid,
            "original_reaction": rsmi,
            "mapped_reaction": mapped,
            "success": success
        })
        print("OK" if success else "FAIL")
    
    # Save results
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    output_csv = OUT_DIR / "rdt_fixed_mappings.csv"
    output_json = OUT_DIR / "rdt_fixed_mappings.json"
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    
    with output_json.open("w", encoding="utf-8") as f:
        for result in results:
            f.write(json.dumps(result, ensure_ascii=False) + "\n")
    
    success_count = sum(1 for r in results if r["success"])
    print(f"\n[OK] RDT (fixed): {success_count}/{len(results)} successful")
    print(f"[OK] Saved to: {output_csv}")
    print(f"[OK] Saved to: {output_json}")


if __name__ == "__main__":
    main()

