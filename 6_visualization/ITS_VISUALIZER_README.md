# ITS Visualizer Documentation

A comprehensive toolkit for visualizing Imaginary Transition State (ITS) graphs, reaction centers, and reaction templates inspired by the SynTemp framework.

## Overview

The `ITSVisualizer` class provides visualization capabilities similar to those shown in the SynTemp paper (Figure 1), allowing you to create publication-quality diagrams of:

- ITS graphs with color-coded bond changes
- Reaction centers and extended reaction centers
- Complete reaction visualizations (Reactants → ITS → Products)
- Hierarchical clustering of reaction templates

## Quick Start

```python
from its_visualizer import ITSVisualizer, visualize_its
import networkx as nx

# Create your ITS graph
its_graph = nx.Graph()
its_graph.add_nodes_from([
    (0, {'element': 'C'}),
    (1, {'element': 'O'}),
    (2, {'element': 'C'}),
])
its_graph.add_edges_from([
    (0, 1, {'standard_order': 1}),   # Breaking bond (red)
    (1, 2, {'standard_order': -1}),  # Forming bond (green)
])

# Quick visualization
fig = visualize_its(its_graph, title="My ITS Graph")
```

## Key Features

### 1. ITS Graph Visualization

```python
visualizer = ITSVisualizer(seed=42)

fig = visualizer.visualize_its_graph(
    its_graph,
    title="Diels-Alder ITS",
    figsize=(10, 8),
    show_atom_labels=True,      # Show element symbols
    show_atom_numbers=True,     # Show atom mapping numbers
    show_legend=True,            # Display edge color legend
    node_size=700,
    layout="spring"              # or "kamada_kawai", "circular", "planar"
)
```

**Edge Color Coding:**
- 🔴 **Red** = Breaking bonds
- 🟢 **Green** = Forming bonds
- 🔵 **Blue** = Changing bonds (e.g., single → double)
- ⚫ **Black** = Unchanged bonds

### 2. Reaction Center Highlighting

```python
# Define the reaction center (atoms involved in bond changes)
reaction_center = {0, 1, 4, 5}

# Visualize with highlighting
fig = visualizer.visualize_reaction_center(
    its_graph,
    reaction_center_nodes=reaction_center,
    radius=0,              # 0 = only RC, 1+ = extended context
    highlight_color="#FFD700",
    title="Reaction Center (r=0)"
)
```

### 3. Complete Reaction Visualization

```python
from its_visualizer import visualize_reaction

fig = visualize_reaction(
    reactants_graph,
    products_graph,
    its_graph,
    figsize=(18, 6),
    show_atom_numbers=True
)
```

Creates a three-panel view: **Reactants | ITS Graph | Products**

### 4. Custom Styling

```python
# Custom edge colors
custom_colors = {
    "unchanged": "#666666",
    "breaking": "#FF4444",
    "forming": "#44FF44",
    "changing": "#4444FF",
}

visualizer = ITSVisualizer(edge_colors=custom_colors, seed=42)
```

## Edge Attributes

The visualizer supports two edge attribute formats:

### Format 1: `standard_order` (from ChemicalGraphVisualizer)

```python
# standard_order > 0: breaking bond
# standard_order < 0: forming bond  
# standard_order = 0: unchanged bond

its_graph.add_edge(0, 1, standard_order=1)   # Breaking
its_graph.add_edge(1, 2, standard_order=-1)  # Forming
its_graph.add_edge(2, 3, standard_order=0)   # Unchanged
```

### Format 2: `order` tuple (SynTemp style)

```python
# order = (reactant_order, product_order)
# 0 = no bond, 1 = single, 2 = double, 3 = triple

its_graph.add_edge(0, 1, order=(2, 1))  # Double → Single (changing)
its_graph.add_edge(1, 2, order=(0, 1))  # No bond → Single (forming)
its_graph.add_edge(2, 3, order=(1, 0))  # Single → No bond (breaking)
its_graph.add_edge(3, 4, order=(1, 1))  # Single → Single (unchanged)
```

## Node Attributes

```python
# Required
its_graph.add_node(0, element='C')

# Optional
its_graph.add_node(0, atom_map=10)  # For atom numbering
```

## Integration with Existing Code

Works seamlessly with your existing `ChemicalGraphVisualizer`:

```python
from chemical_graph_visualizer import ChemicalGraphVisualizer
from its_visualizer import ITSVisualizer

# Use both visualizers
chem_viz = ChemicalGraphVisualizer(seed=42)
its_viz = ITSVisualizer(seed=42)

# Create comparison
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
chem_viz.its_vis(graph, ax=axes[0])
its_viz.visualize_its_graph(graph, ax=axes[1], show_legend=True)
```

## Layout Algorithms

- **`spring`** (default) - Force-directed layout, good for most cases
- **`kamada_kawai`** - Minimizes energy, better for symmetric molecules
- **`circular`** - Nodes arranged in a circle
- **`planar`** - For planar graphs only

```python
fig = visualizer.visualize_its_graph(
    its_graph,
    layout="kamada_kawai",
    # Pass kwargs to layout algorithm
    scale=2.0,
    center=(0, 0)
)
```

## Extended Reaction Centers

Visualize context around the reaction center at different radii:

```python
reaction_center = {2, 3, 4}

# r=0: Only atoms in reaction center
fig0 = visualizer.visualize_reaction_center(its_graph, reaction_center, radius=0)

# r=1: RC + immediate neighbors
fig1 = visualizer.visualize_reaction_center(its_graph, reaction_center, radius=1)

# r=2: RC + neighbors within 2 bonds
fig2 = visualizer.visualize_reaction_center(its_graph, reaction_center, radius=2)
```

## Saving Figures

```python
# High-resolution for publications
fig = visualize_its(its_graph, title="My Reaction")
fig.savefig("its_graph.png", dpi=300, bbox_inches='tight')
fig.savefig("its_graph.pdf", bbox_inches='tight')  # Vector format

# Save multiple reactions
for rxn_id, its_g in reaction_data.items():
    fig = visualize_its(its_g, title=f"Reaction {rxn_id}")
    fig.savefig(f"figures/rxn_{rxn_id}.png", dpi=300)
    plt.close(fig)
```

## Advanced Example: iPKS Reaction Rules

```python
import pickle
from its_visualizer import ITSVisualizer

# Load your iPKS ITS graphs
with open('ipks_its_graphs.pkl', 'rb') as f:
    ipks_data = pickle.load(f)

visualizer = ITSVisualizer(seed=42)

# Visualize all reactions
for rule_id, data in ipks_data.items():
    # Extract graphs
    reactants = data['reactants_graph']
    products = data['products_graph']
    its = data['its_graph']
    rc_nodes = data['reaction_center']
    
    # Create multi-panel figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Panel A: Full ITS
    visualizer.visualize_its_graph(
        its, ax=axes[0, 0],
        title=f"Rule {rule_id}: ITS Graph"
    )
    
    # Panel B: Reaction center
    visualizer.visualize_reaction_center(
        its, rc_nodes, radius=0,
        ax=axes[0, 1],
        title="Reaction Center"
    )
    
    # Panel C: Extended (r=1)
    visualizer.visualize_reaction_center(
        its, rc_nodes, radius=1,
        ax=axes[1, 0],
        title="Extended Context (r=1)"
    )
    
    # Panel D: Complete reaction
    visualizer.visualize_reaction_triple(
        reactants, products, its,
        ax=axes[1, 1]
    )
    
    plt.tight_layout()
    plt.savefig(f"ipks_rules/rule_{rule_id}.png", dpi=300)
    plt.close()
```

## Comparison with SynTemp Paper

This visualizer recreates the visualization style from SynTemp (Phan et al., 2025):

| Figure Element | Method |
|----------------|--------|
| Panel A: Reaction with AAM | `visualize_reaction_triple()` |
| Panel B: Complete ITS | `visualize_its_graph()` |
| Panel C: Reaction Center | `visualize_reaction_center(radius=0)` |
| Panel D: Clustering | `visualize_hierarchical_clustering()` |

## API Reference

### `ITSVisualizer`

Main class for ITS visualization.

**Methods:**
- `visualize_its_graph()` - Visualize ITS with colored edges
- `visualize_reaction_center()` - Highlight reaction center
- `visualize_reaction_triple()` - Show reactants → ITS → products
- `visualize_hierarchical_clustering()` - Show template clustering

### Convenience Functions

- `visualize_its(its_graph, **kwargs)` - Quick ITS visualization
- `visualize_reaction(reactants, products, its, **kwargs)` - Quick reaction view

## Requirements

- `networkx` - Graph data structure
- `matplotlib` - Plotting
- `numpy` - Numerical operations
- Optional: `pygraphviz` for better hierarchical layouts

## Tips

1. **For publication figures:**
   - Use `dpi=300` or higher
   - Save as PDF for vector graphics
   - Set `figsize` appropriately (e.g., `(10, 8)` for single column)

2. **For large graphs:**
   - Increase `node_size` and `figsize`
   - Try different layout algorithms
   - Consider showing only reaction center

3. **For presentations:**
   - Use larger `node_size` (800-1000)
   - Increase font sizes
   - Use high contrast colors

## Troubleshooting

**Issue:** Nodes overlapping
- Try different layout: `layout="kamada_kawai"`
- Increase `figsize`
- Adjust layout parameters

**Issue:** Colors not showing correctly
- Check edge attributes (`order` or `standard_order`)
- Verify attribute format matches documentation
- Use custom colors if needed

**Issue:** Slow rendering
- Use simpler layouts (`circular` instead of `spring`)
- Reduce `node_size`
- Visualize reaction center only, not full ITS

## Citation

If you use this visualizer in your research, please cite:

- **SynTemp Paper:** Phan, T.-L., et al. (2025). "SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases." *J. Chem. Inf. Model.*, 65, 2882-2896.

## License

Part of the IteraPath project.

