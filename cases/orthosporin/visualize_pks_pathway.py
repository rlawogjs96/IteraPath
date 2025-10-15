#!/usr/bin/env python3
"""
PKS Domain Architecture + Biosynthetic Pathway Visualization
Creates a comprehensive diagram showing:
1. PKS module domain architecture (KS, AT, DH, KR, ACP)
2. Intermediate structures at each step
3. Substrate additions (Mal-CoA, Ace-CoA)
4. Tailoring reactions
"""
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io
import numpy as np

# Paths
BASE_DIR = Path(__file__).resolve().parent
ROUTE_FILE = BASE_DIR / "route.json"
AAM_FILE = BASE_DIR.parent.parent / "3_templates" / "data" / "aam.ndjson"

# Color scheme (matching the reference image)
COLORS = {
    'KS': '#C85A3E',      # Brownish-red (active)
    'AT': '#C85A3E',      # Brownish-red (active)
    'DH': '#FFFFFF',      # White
    'KR': '#FFFFFF',      # White
    'ACP': '#FFFFFF',     # White
    'PT': '#C85A3E',      # Brownish-red (active for cyclization)
    'TE': '#FFFFFF',      # White
}

DOMAIN_LABELS = {
    'KS': 'KS',
    'AT': 'AT',
    'DH': 'DH',
    'KR': 'KR',
    'ACP': 'ACP',
    'PT': 'PT',
    'TE': 'TE',
}


def load_route_data():
    """Load route.json"""
    with ROUTE_FILE.open('r') as f:
        return json.load(f)


def load_intermediate_smiles():
    """Load intermediate structures from aam.ndjson"""
    intermediates = {}
    
    with AAM_FILE.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            
            if not rec['reaction_id'].startswith('BGC1000006'):
                continue
            
            # Parse mapped SMARTS to get product structure
            rxn_smarts = rec['aam']
            parts = rxn_smarts.split('>>')
            if len(parts) == 2:
                product_smarts = parts[1]
                
                # Remove atom mapping numbers to get clean SMILES
                mol = Chem.MolFromSmarts(product_smarts)
                if mol:
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    product_smiles = Chem.MolToSmiles(mol)
                    
                    # Store with step number
                    step_num = rec['reaction_id'].split('_')[-1]
                    intermediates[step_num] = product_smiles
    
    # Add starter unit (m0) - Acetyl-CoA
    intermediates['m0'] = 'CC(=O)[S]'
    
    return intermediates


def draw_molecule(smiles, img_size=(300, 200)):
    """Draw molecule using RDKit"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2
    opts.padding = 0.1
    
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    png_data = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_data))
    
    return img


def draw_domain_circle(ax, x, y, domain_type, radius=0.3):
    """Draw a single domain as a circle"""
    color = COLORS.get(domain_type, '#FFFFFF')
    
    circle = patches.Circle(
        (x, y), radius,
        facecolor=color,
        edgecolor='black',
        linewidth=2,
        zorder=2
    )
    ax.add_patch(circle)
    
    # Add label
    ax.text(
        x, y, DOMAIN_LABELS.get(domain_type, domain_type),
        ha='center', va='center',
        fontsize=11, fontweight='bold',
        zorder=3
    )


def draw_module(ax, x_start, y, domains, module_idx):
    """Draw a complete PKS module"""
    domain_spacing = 0.8
    x = x_start
    
    for i, domain in enumerate(domains):
        draw_domain_circle(ax, x, y, domain)
        
        # Draw connection to next domain
        if i < len(domains) - 1:
            ax.plot([x + 0.35, x + domain_spacing - 0.35], [y, y],
                   'k-', linewidth=1.5, zorder=1)
        
        x += domain_spacing
    
    return x  # Return end position


def add_substrate_label(ax, x, y, substrate):
    """Add substrate label (Mal-CoA, Ace-CoA)"""
    ax.text(x, y, substrate, ha='center', va='bottom',
           fontsize=9, style='italic', color='#555')


def add_reaction_arrow(ax, x1, y1, x2, y2):
    """Add arrow between modules"""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
               arrowprops=dict(arrowstyle='->', lw=2, color='black'))


def create_pks_pathway_diagram():
    """Create the complete PKS pathway diagram"""
    print("=" * 70)
    print("Creating PKS Pathway Diagram")
    print("=" * 70)
    
    # Load data
    route_data = load_route_data()
    intermediates = load_intermediate_smiles()
    
    print(f"\nLoaded {len(route_data['steps'])} reaction steps")
    print(f"Loaded {len(intermediates)} intermediate structures")
    
    # Create figure with single large canvas
    fig, ax = plt.subplots(1, 1, figsize=(18, 22))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 24)
    ax.axis('off')
    
    # Title
    ax.text(8, 23.5, 'BGC1000006: Polyketide Biosynthesis Pathway',
           ha='center', fontsize=18, fontweight='bold')
    ax.text(8, 23, 'Tailoring Method',
           ha='center', fontsize=16, fontweight='bold')
    
    # Current y position
    y = 21.5
    y_step = 3.5
    
    # Process each step
    for step_idx, step in enumerate(route_data['steps']):
        step_type = step['step_type']
        domains_str = step['domains']
        domains = [d.strip() for d in domains_str.split(',')]
        at_substrate = step.get('at_substrate', '')
        reaction_id = step['reaction_id']
        
        step_num_from = reaction_id.split('_')[-3]
        step_num_to = reaction_id.split('_')[-1]
        
        print(f"\n  Step {step_idx + 1}: {step_num_from} → {step_num_to}")
        print(f"    Domains: {domains}")
        print(f"    Type: {step_type}")
        
        # Draw reactant structure (left side)
        if step_num_from in intermediates:
            smiles_from = intermediates[step_num_from]
            mol_img = draw_molecule(smiles_from, img_size=(280, 180))
            
            if mol_img:
                img_array = np.array(mol_img)
                extent = [0.5, 2.8, y - 1.2, y + 0.3]
                ax.imshow(img_array, extent=extent, aspect='auto', zorder=1)
        
        # Draw module (center)
        x_module_start = 4.5
        end_x = draw_module(ax, x_module_start, y - 0.5, domains, step_idx)
        
        # Add substrate label
        if at_substrate == 'mal':
            ax.text(x_module_start + 1.2, y + 0.5, 'Mal-CoA',
                   ha='center', va='bottom', fontsize=10, style='italic')
        elif at_substrate == 'ace':
            ax.text(x_module_start + 1.2, y + 0.5, 'Ace-CoA',
                   ha='center', va='bottom', fontsize=10, style='italic')
        
        # Add PT/TE label for cyclization
        if step_type == 'CLOSURE':
            mechanism = step.get('mechanism', '')
            if 'PT' in mechanism:
                ax.text(x_module_start + 0.5, y + 0.5, 'PT + TE',
                       ha='center', va='bottom',
                       fontsize=10, style='italic', color='#555')
        
        # Add arrow
        ax.annotate('', xy=(11.5, y - 0.5), xytext=(end_x + 0.5, y - 0.5),
                   arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))
        
        # Draw product structure (right side)
        if step_num_to in intermediates:
            smiles_to = intermediates[step_num_to]
            mol_img = draw_molecule(smiles_to, img_size=(280, 180))
            
            if mol_img:
                img_array = np.array(mol_img)
                extent = [12, 14.3, y - 1.2, y + 0.3]
                ax.imshow(img_array, extent=extent, aspect='auto', zorder=1)
        
        # Move to next row
        y -= y_step
    
    plt.tight_layout()
    
    # Save
    output_path = BASE_DIR / 'BGC1000006_pks_pathway_diagram.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n[SAVED] {output_path}")
    
    plt.show()


if __name__ == "__main__":
    create_pks_pathway_diagram()

