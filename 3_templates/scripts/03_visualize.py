#!/usr/bin/env python3
"""
Complete visualization of BGC1000006 all reaction steps
Creates a single figure with all 6 reactions showing atom mapping
"""
import json
import matplotlib.pyplot as plt
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io
import numpy as np

# Paths
BASE = Path(__file__).resolve().parent
DATA = BASE.parent / "data"
AAM_FILE = DATA / "aam.ndjson"
OUTPUT_DIR = BASE.parent / "data"


def visualize_reaction_with_aam(mapped_smarts, title, img_size=(2000, 800)):
    """
    Visualize a reaction with atom mapping numbers using RDKit.
    Shows element:number format (e.g., C:1, O:2) for better clarity.
    """
    # Parse reaction
    rxn = rdChemReactions.ReactionFromSmarts(mapped_smarts, useSmiles=True)
    if rxn is None:
        return None
    
    # Preprocess
    rdChemReactions.PreprocessReaction(rxn)
    
    # Add atom map numbers as labels in "Element:Number" format
    for mol in list(rxn.GetReactants()) + list(rxn.GetProducts()):
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                # Format as "Element:Number" (e.g., "C:1", "O:2")
                symbol = atom.GetSymbol()
                atom.SetProp("atomLabel", f"{symbol}:{map_num}")
    
    # Create drawer
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    
    # Configure options
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3.5
    opts.padding = 0.02  # Minimal padding for tight layout
    opts.addStereoAnnotation = True
    opts.annotationFontScale = 1.5  # Increase font size for atom labels
    
    # Draw reaction with color highlighting (will be yellow by default)
    drawer.DrawReaction(rxn, highlightByReactant=True)
    drawer.FinishDrawing()
    
    # Get PNG data and convert to PIL Image
    png_data = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_data))
    
    # Convert yellow highlights to grapefruit/coral color
    # Yellow RGB: (255, 255, 0) → Grapefruit/Coral RGB: (255, 120, 100)
    img_array = np.array(img)
    
    # Define color ranges for yellow highlighting (with some tolerance)
    # Target yellow range in RDKit highlighting
    yellow_min = np.array([240, 240, 0])
    yellow_max = np.array([255, 255, 100])
    
    # Grapefruit/coral target color
    grapefruit_color = np.array([255, 120, 100])
    
    # Create mask for yellow pixels
    mask = np.all((img_array >= yellow_min) & (img_array <= yellow_max), axis=-1)
    
    # Replace yellow with grapefruit
    img_array[mask] = grapefruit_color
    
    # Convert back to PIL Image
    img = Image.fromarray(img_array)
    
    return img


def main():
    """Create complete visualization for BGC1000006"""
    
    print("=" * 70)
    print("BGC1000006 Complete Pathway Visualization")
    print("=" * 70)
    
    # Read AAM data
    print(f"\nReading data from: {AAM_FILE}")
    reactions = []
    
    with AAM_FILE.open('r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                rec = json.loads(line)
                if rec['reaction_id'].startswith('BGC1000006'):
                    reactions.append(rec)
    
    print(f"Found {len(reactions)} reactions for BGC1000006:")
    for r in sorted(reactions, key=lambda x: x['reaction_id']):
        print(f"  - {r['reaction_id']}")
    
    # Create figure with tight spacing
    print("\nCreating visualization...")
    fig, axes = plt.subplots(3, 2, figsize=(26, 30))
    
    fig.suptitle(
        'BGC1000006: Complete Polyketide Biosynthesis Pathway\n'
        'All 6 reaction steps with atom-to-atom mapping (Element:Number)',
        fontsize=28, fontweight='bold', y=0.997
    )
    
    # Flatten axes
    axes = axes.flatten()
    
    # Visualize each reaction
    for idx, rec in enumerate(sorted(reactions, key=lambda x: x['reaction_id'])):
        reaction_id = rec['reaction_id']
        mapped_smarts = rec['aam']
        unmapped_rxn = rec['reactions']
        
        print(f"  Processing {reaction_id}...")
        
        # Create visualization with larger size for bigger fonts
        img = visualize_reaction_with_aam(mapped_smarts, reaction_id, img_size=(2000, 800))
        
        if img is None:
            print(f"    [ERROR] Failed to create image")
            continue
        
        # Display in subplot
        ax = axes[idx]
        ax.imshow(img)
        ax.axis('off')
        
        # Add title with much larger font
        step_name = reaction_id.replace('BGC1000006_', '')
        ax.set_title(
            f'{step_name}\n{unmapped_rxn}',
            fontsize=16, fontweight='bold', pad=12, family='monospace'
        )
    
    # Adjust layout - very tight margins
    plt.tight_layout(rect=[0, 0, 1, 0.997], h_pad=0.2, w_pad=0.3)
    
    # Save
    output_path = OUTPUT_DIR / 'BGC1000006_complete_pathway_aam.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\n[SAVED] {output_path}")
    
    # Also save high-res version
    output_path_hires = OUTPUT_DIR / 'BGC1000006_complete_pathway_aam_hires.png'
    plt.savefig(output_path_hires, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[SAVED] {output_path_hires} (high resolution)")
    
    plt.close()
    
    print("\n" + "=" * 70)
    print("Visualization complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
