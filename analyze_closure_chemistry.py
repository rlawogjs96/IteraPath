"""
Analyze CLOSURE reactions to understand the chemistry behind H2O production
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd
import re

df = pd.read_csv('C:/IteraPath/2_balancing/data/balance_results.csv')
closures = df[df['step_type'] == 'CLOSURE'].copy()

print("="*80)
print("CLOSURE REACTION CHEMISTRY ANALYSIS")
print("="*80)

for _, row in closures.iterrows():
    rxn_id = row['reaction_id']
    rxn = row['reactions']
    
    # Split reaction
    left, right = rxn.split('>>')
    reactants = left.split('.')
    products = right.split('.')
    
    # Get the chain (first reactant, excluding H2O)
    chain = [r for r in reactants if r != 'O'][0]
    # Get the aromatic product (first product, excluding [SH] and O)
    aromatic = [p for p in products if p not in ['[SH]', 'O']][0]
    
    # Count H2O on products
    h2o_count = products.count('O')
    
    # Analyze structures
    mol_chain = Chem.MolFromSmiles(chain)
    mol_aromatic = Chem.MolFromSmiles(aromatic)
    
    if mol_chain and mol_aromatic:
        # Get molecular formulas
        formula_chain = rdMolDescriptors.CalcMolFormula(mol_chain)
        formula_aromatic = rdMolDescriptors.CalcMolFormula(mol_aromatic)
        
        def parse_formula(f):
            elements = {}
            for match in re.finditer(r'([A-Z][a-z]?)(\d*)', f):
                elem = match.group(1)
                count = int(match.group(2)) if match.group(2) else 1
                elements[elem] = elements.get(elem, 0) + count
            return elements
        
        elem_chain = parse_formula(formula_chain)
        elem_aromatic = parse_formula(formula_aromatic)
        
        # Calculate H difference (how many H atoms lost)
        h_chain = elem_chain.get('H', 0)
        h_aromatic = elem_aromatic.get('H', 0)
        h_lost = h_chain - h_aromatic
        
        # Each dehydration removes 2H + 1O = H2O
        # Lactonization removes 2H + 1O = H2O (from -COOH + -OH -> ester + H2O)
        dehydrations_needed = h_lost // 2
        
        print(f"\n{rxn_id}:")
        print(f"  Chain: {chain}")
        print(f"    Formula: {formula_chain}")
        print(f"    H atoms: {h_chain}")
        print(f"  ")
        print(f"  Aromatic: {aromatic}")
        print(f"    Formula: {formula_aromatic}")
        print(f"    H atoms: {h_aromatic}")
        print(f"  ")
        print(f"  Transformation:")
        print(f"    H atoms lost: {h_lost}")
        print(f"    Dehydrations needed: {dehydrations_needed} (each removes 2H as H2O)")
        print(f"  ")
        print(f"  H2O balance:")
        print(f"    Found in products: {h2o_count} x H2O")
        print(f"    Expected: {dehydrations_needed} x H2O")
        print(f"  ")
        if h2o_count == dehydrations_needed:
            print(f"  CHECK: CONSISTENT")
        else:
            print(f"  CHECK: DISCREPANCY (got {h2o_count}, expected {dehydrations_needed})")

print("\n" + "="*80)
print("INTERPRETATION:")
print("="*80)
print("Each new C=C double bond formed requires one dehydration (removal of H2O).")
print("Aromatic ring formation requires multiple dehydrations to create the π-system.")
print("Lactonization (ester ring formation) releases one additional H2O.")
print("="*80)

