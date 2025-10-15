"""
ITS (Imaginary Transition State) Graph Visualizer

This module provides visualization tools for ITS graphs, reaction centers,
and extended reaction centers inspired by the SynTemp framework.
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from typing import Optional, Dict, List, Tuple, Set
import numpy as np
from IPython.display import display, SVG
import io


# Chemical element color scheme (CPK colors)
ELEMENT_COLORS = {
    "H": "#FFFFFF",  # White
    "C": "#909090",  # Gray
    "N": "#3050F8",  # Blue
    "O": "#FF0D0D",  # Red
    "F": "#90E050",  # Green
    "Cl": "#1FF01F",  # Green
    "Br": "#A62929",  # Dark Red/Brown
    "I": "#940094",  # Purple
    "P": "#FF8000",  # Orange
    "S": "#FFFF30",  # Yellow
}

# Edge colors for ITS graphs
EDGE_COLORS = {
    "unchanged": "#000000",  # Black - bonds that don't change
    "breaking": "#FF0000",   # Red - bonds being broken
    "forming": "#00FF00",    # Green - bonds being formed
    "changing": "#0000FF",   # Blue - bonds changing order
}


class ITSVisualizer:
    """
    Visualizer for Imaginary Transition State (ITS) graphs.
    
    Provides methods to visualize:
    - ITS graphs with colored edges showing bond changes
    - Reaction centers
    - Extended reaction centers with different radii
    - Hierarchical clustering of reaction templates
    """
    
    def __init__(
        self,
        element_colors: Optional[Dict[str, str]] = None,
        edge_colors: Optional[Dict[str, str]] = None,
        seed: Optional[int] = 42,
    ):
        """
        Initialize the ITS visualizer.
        
        Parameters:
        -----------
        element_colors : dict, optional
            Mapping of element symbols to hex color codes
        edge_colors : dict, optional
            Mapping of edge types to hex color codes
        seed : int, optional
            Random seed for reproducible layouts
        """
        self.element_colors = element_colors or ELEMENT_COLORS
        self.edge_colors = edge_colors or EDGE_COLORS
        self.seed = seed
    
    def visualize_its_graph(
        self,
        G: nx.Graph,
        node_size: int = 500,
        figsize: Tuple[int, int] = (10, 8),
        show_atom_labels: bool = True,
        show_atom_numbers: bool = False,
        show_legend: bool = True,
        title: str = "ITS Graph",
        layout: str = "spring",
        ax: Optional[plt.Axes] = None,
        **layout_kwargs
    ) -> plt.Figure:
        """
        Visualize an ITS graph with colored edges indicating bond changes.
        
        Parameters:
        -----------
        G : nx.Graph
            ITS graph to visualize
        node_size : int
            Size of nodes in the visualization
        figsize : tuple
            Figure size (width, height)
        show_atom_labels : bool
            Show element symbols on nodes
        show_atom_numbers : bool
            Show atom numbering on nodes
        show_legend : bool
            Display legend for edge colors
        title : str
            Graph title
        layout : str
            Layout algorithm ('spring', 'kamada_kawai', 'circular', 'planar')
        ax : plt.Axes, optional
            Matplotlib axes to draw on
        **layout_kwargs
            Additional arguments for layout algorithm
        
        Returns:
        --------
        fig : matplotlib.Figure
            The figure object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        
        # Compute layout
        pos = self._compute_layout(G, layout, **layout_kwargs)
        
        # Get node colors based on elements
        node_colors = [
            self.element_colors.get(G.nodes[node].get("element", "C"), "#909090")
            for node in G.nodes()
        ]
        
        # Classify edges by type
        edge_lists = self._classify_edges(G)
        
        # Draw edges with different colors
        for edge_type, edges in edge_lists.items():
            if edges:
                nx.draw_networkx_edges(
                    G, pos,
                    edgelist=edges,
                    edge_color=self.edge_colors.get(edge_type, "#000000"),
                    width=2.5,
                    ax=ax,
                    alpha=0.8
                )
        
        # Draw nodes
        nx.draw_networkx_nodes(
            G, pos,
            node_color=node_colors,
            node_size=node_size,
            edgecolors='black',
            linewidths=1.5,
            ax=ax
        )
        
        # Draw labels
        if show_atom_labels or show_atom_numbers:
            labels = self._create_node_labels(G, show_atom_labels, show_atom_numbers)
            nx.draw_networkx_labels(
                G, pos,
                labels=labels,
                font_size=10,
                font_weight='bold',
                ax=ax
            )
        
        # Add legend
        if show_legend:
            self._add_edge_legend(ax)
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
        plt.tight_layout()
        
        return fig
    
    def visualize_reaction_center(
        self,
        its_graph: nx.Graph,
        reaction_center_nodes: Set,
        radius: int = 0,
        figsize: Tuple[int, int] = (10, 8),
        highlight_color: str = "#FFD700",
        title: str = "Reaction Center",
        ax: Optional[plt.Axes] = None,
    ) -> plt.Figure:
        """
        Visualize the reaction center within an ITS graph.
        
        Parameters:
        -----------
        its_graph : nx.Graph
            The full ITS graph
        reaction_center_nodes : set
            Nodes in the reaction center
        radius : int
            Extension radius around reaction center
        figsize : tuple
            Figure size
        highlight_color : str
            Color to highlight reaction center nodes
        title : str
            Graph title
        ax : plt.Axes, optional
            Matplotlib axes
        
        Returns:
        --------
        fig : matplotlib.Figure
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        
        # Extend reaction center by radius
        extended_nodes = self._extend_nodes_by_radius(
            its_graph, reaction_center_nodes, radius
        )
        
        # Compute layout
        pos = self._compute_layout(its_graph, "spring")
        
        # Classify nodes
        core_nodes = list(reaction_center_nodes)
        extended_only = list(extended_nodes - reaction_center_nodes)
        other_nodes = list(set(its_graph.nodes()) - extended_nodes)
        
        # Get node colors
        def get_node_color(node, node_set):
            if node in node_set:
                return highlight_color
            element = its_graph.nodes[node].get("element", "C")
            return self.element_colors.get(element, "#909090")
        
        # Draw edges
        edge_lists = self._classify_edges(its_graph)
        for edge_type, edges in edge_lists.items():
            if edges:
                # Highlight edges in reaction center
                rc_edges = [(u, v) for u, v in edges if u in reaction_center_nodes and v in reaction_center_nodes]
                other_edges = [(u, v) for u, v in edges if (u, v) not in rc_edges]
                
                if rc_edges:
                    nx.draw_networkx_edges(
                        its_graph, pos,
                        edgelist=rc_edges,
                        edge_color=self.edge_colors.get(edge_type, "#000000"),
                        width=4.0,
                        ax=ax,
                        alpha=1.0
                    )
                if other_edges:
                    nx.draw_networkx_edges(
                        its_graph, pos,
                        edgelist=other_edges,
                        edge_color=self.edge_colors.get(edge_type, "#000000"),
                        width=1.5,
                        ax=ax,
                        alpha=0.4
                    )
        
        # Draw nodes in layers
        if core_nodes:
            nx.draw_networkx_nodes(
                its_graph, pos,
                nodelist=core_nodes,
                node_color=[get_node_color(n, reaction_center_nodes) for n in core_nodes],
                node_size=600,
                edgecolors='red',
                linewidths=3.0,
                ax=ax,
                label="Reaction Center"
            )
        
        if extended_only:
            nx.draw_networkx_nodes(
                its_graph, pos,
                nodelist=extended_only,
                node_color=[self.element_colors.get(its_graph.nodes[n].get("element", "C"), "#909090") for n in extended_only],
                node_size=500,
                edgecolors='orange',
                linewidths=2.0,
                ax=ax,
                label=f"Extended (r={radius})" if radius > 0 else None
            )
        
        if other_nodes:
            nx.draw_networkx_nodes(
                its_graph, pos,
                nodelist=other_nodes,
                node_color=[self.element_colors.get(its_graph.nodes[n].get("element", "C"), "#909090") for n in other_nodes],
                node_size=400,
                edgecolors='black',
                linewidths=1.0,
                ax=ax,
                alpha=0.3
            )
        
        # Draw labels
        labels = {node: its_graph.nodes[node].get("element", "") for node in its_graph.nodes()}
        nx.draw_networkx_labels(
            its_graph, pos,
            labels=labels,
            font_size=10,
            font_weight='bold',
            ax=ax
        )
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
        if radius > 0 or core_nodes:
            ax.legend(loc='upper right')
        plt.tight_layout()
        
        return fig
    
    def visualize_reaction_triple(
        self,
        reactants_graph: nx.Graph,
        products_graph: nx.Graph,
        its_graph: nx.Graph,
        figsize: Tuple[int, int] = (18, 6),
        show_atom_numbers: bool = False,
        titles: Optional[Tuple[str, str, str]] = None,
    ) -> plt.Figure:
        """
        Visualize reactants, ITS graph, and products side by side.
        
        Parameters:
        -----------
        reactants_graph : nx.Graph
            Graph of reactants
        products_graph : nx.Graph
            Graph of products
        its_graph : nx.Graph
            ITS graph
        figsize : tuple
            Figure size
        show_atom_numbers : bool
            Show atom mapping numbers
        titles : tuple of str, optional
            Custom titles for (reactants, ITS, products)
        
        Returns:
        --------
        fig : matplotlib.Figure
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        if titles is None:
            titles = ("Reactants", "ITS Graph", "Products")
        
        # Visualize reactants
        self._visualize_simple_graph(
            reactants_graph,
            ax=axes[0],
            title=titles[0],
            show_atom_numbers=show_atom_numbers
        )
        
        # Visualize ITS
        self.visualize_its_graph(
            its_graph,
            ax=axes[1],
            title=titles[1],
            show_atom_numbers=show_atom_numbers,
            show_legend=False
        )
        
        # Visualize products
        self._visualize_simple_graph(
            products_graph,
            ax=axes[2],
            title=titles[2],
            show_atom_numbers=show_atom_numbers
        )
        
        # Add legend to the right side
        self._add_edge_legend(axes[1], loc='upper right')
        
        plt.tight_layout()
        return fig
    
    def visualize_hierarchical_clustering(
        self,
        cluster_tree: Dict,
        max_depth: int = 3,
        figsize: Tuple[int, int] = (14, 10),
        node_size: int = 300,
    ) -> plt.Figure:
        """
        Visualize hierarchical clustering of reaction templates.
        
        Parameters:
        -----------
        cluster_tree : dict
            Hierarchical tree structure of clusters
        max_depth : int
            Maximum depth to visualize
        figsize : tuple
            Figure size
        node_size : int
            Size of nodes in tree
        
        Returns:
        --------
        fig : matplotlib.Figure
        """
        fig, ax = plt.subplots(figsize=figsize)
        
        # Build tree graph
        tree_graph = self._build_tree_graph(cluster_tree, max_depth)
        
        # Use hierarchical layout
        pos = nx.nx_agraph.graphviz_layout(tree_graph, prog='dot') if hasattr(nx, 'nx_agraph') else \
              self._hierarchy_pos(tree_graph)
        
        # Draw tree
        nx.draw_networkx_edges(tree_graph, pos, ax=ax, arrows=True, arrowsize=10)
        nx.draw_networkx_nodes(
            tree_graph, pos,
            node_color='lightblue',
            node_size=node_size,
            edgecolors='black',
            linewidths=2,
            ax=ax
        )
        
        # Draw labels
        labels = nx.get_node_attributes(tree_graph, 'label')
        nx.draw_networkx_labels(tree_graph, pos, labels=labels, font_size=8, ax=ax)
        
        ax.set_title("Hierarchical Clustering of Reaction Templates", fontsize=14, fontweight='bold')
        ax.axis('off')
        plt.tight_layout()
        
        return fig
    
    # Helper methods
    
    def _compute_layout(self, G: nx.Graph, layout: str = "spring", **kwargs) -> Dict:
        """Compute node positions using specified layout algorithm."""
        if self.seed is not None and 'seed' not in kwargs:
            kwargs['seed'] = self.seed
        
        if layout == "spring":
            return nx.spring_layout(G, **kwargs)
        elif layout == "kamada_kawai":
            return nx.kamada_kawai_layout(G, **kwargs)
        elif layout == "circular":
            return nx.circular_layout(G, **kwargs)
        elif layout == "planar":
            return nx.planar_layout(G, **kwargs)
        else:
            return nx.spring_layout(G, **kwargs)
    
    def _classify_edges(self, G: nx.Graph) -> Dict[str, List]:
        """
        Classify edges by their type (unchanged, breaking, forming, changing).
        
        Expected edge attributes:
        - 'order': tuple (reactant_order, product_order) or single value
        - 'standard_order': numeric value indicating change
        """
        edge_lists = {
            "unchanged": [],
            "breaking": [],
            "forming": [],
            "changing": []
        }
        
        for u, v, data in G.edges(data=True):
            # Check for standard_order attribute (from your existing code)
            if 'standard_order' in data:
                order = data['standard_order']
                if order == 0:
                    edge_lists["unchanged"].append((u, v))
                elif order > 0:
                    edge_lists["breaking"].append((u, v))
                elif order < 0:
                    edge_lists["forming"].append((u, v))
            # Check for order tuple
            elif 'order' in data:
                order = data['order']
                if isinstance(order, tuple) and len(order) == 2:
                    r_order, p_order = order
                    if r_order == p_order and r_order > 0:
                        edge_lists["unchanged"].append((u, v))
                    elif r_order > 0 and p_order == 0:
                        edge_lists["breaking"].append((u, v))
                    elif r_order == 0 and p_order > 0:
                        edge_lists["forming"].append((u, v))
                    elif r_order != p_order and r_order > 0 and p_order > 0:
                        edge_lists["changing"].append((u, v))
                else:
                    edge_lists["unchanged"].append((u, v))
            else:
                # Default to unchanged if no order info
                edge_lists["unchanged"].append((u, v))
        
        return edge_lists
    
    def _create_node_labels(
        self,
        G: nx.Graph,
        show_elements: bool = True,
        show_numbers: bool = False
    ) -> Dict:
        """Create node labels with element symbols and/or numbers."""
        labels = {}
        for node in G.nodes():
            label_parts = []
            
            if show_elements:
                element = G.nodes[node].get("element", "")
                if element:
                    label_parts.append(element)
            
            if show_numbers:
                # Check for atom mapping number
                atom_num = G.nodes[node].get("atom_map", node)
                label_parts.append(str(atom_num))
            
            labels[node] = "\n".join(label_parts) if label_parts else ""
        
        return labels
    
    def _add_edge_legend(self, ax: plt.Axes, loc: str = 'upper right') -> None:
        """Add legend explaining edge colors."""
        legend_elements = [
            mpatches.Patch(color=self.edge_colors["breaking"], label='Breaking bonds'),
            mpatches.Patch(color=self.edge_colors["forming"], label='Forming bonds'),
            mpatches.Patch(color=self.edge_colors["changing"], label='Changing bonds'),
            mpatches.Patch(color=self.edge_colors["unchanged"], label='Unchanged bonds'),
        ]
        ax.legend(handles=legend_elements, loc=loc, fontsize=9)
    
    def _visualize_simple_graph(
        self,
        G: nx.Graph,
        ax: plt.Axes,
        title: str = "",
        show_atom_numbers: bool = False
    ) -> None:
        """Visualize a simple molecular graph (reactants or products)."""
        pos = self._compute_layout(G, "spring")
        
        node_colors = [
            self.element_colors.get(G.nodes[node].get("element", "C"), "#909090")
            for node in G.nodes()
        ]
        
        nx.draw_networkx_edges(G, pos, ax=ax, width=2.0)
        nx.draw_networkx_nodes(
            G, pos,
            node_color=node_colors,
            node_size=500,
            edgecolors='black',
            linewidths=1.5,
            ax=ax
        )
        
        labels = self._create_node_labels(G, show_elements=True, show_numbers=show_atom_numbers)
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_weight='bold', ax=ax)
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
    
    def _extend_nodes_by_radius(
        self,
        G: nx.Graph,
        core_nodes: Set,
        radius: int
    ) -> Set:
        """Extend a set of nodes by a given radius in the graph."""
        extended = set(core_nodes)
        current_layer = set(core_nodes)
        
        for _ in range(radius):
            next_layer = set()
            for node in current_layer:
                next_layer.update(G.neighbors(node))
            extended.update(next_layer)
            current_layer = next_layer
        
        return extended
    
    def _build_tree_graph(self, cluster_tree: Dict, max_depth: int) -> nx.DiGraph:
        """Build a directed graph from cluster tree structure."""
        tree = nx.DiGraph()
        
        def add_nodes(node_data, parent=None, depth=0):
            if depth > max_depth:
                return
            
            node_id = node_data.get('id', id(node_data))
            label = node_data.get('label', f"Node {node_id}")
            tree.add_node(node_id, label=label)
            
            if parent is not None:
                tree.add_edge(parent, node_id)
            
            for child in node_data.get('children', []):
                add_nodes(child, node_id, depth + 1)
        
        add_nodes(cluster_tree)
        return tree
    
    def _hierarchy_pos(
        self,
        G: nx.DiGraph,
        root=None,
        width=1.0,
        vert_gap=0.2,
        vert_loc=0,
        xcenter=0.5
    ) -> Dict:
        """
        Create hierarchical layout for tree.
        
        Based on Joel Spolsky's algorithm:
        https://rachel53461.wordpress.com/2014/04/20/algorithm-for-drawing-trees/
        """
        if not nx.is_tree(G):
            raise TypeError('Cannot use hierarchy_pos on a graph that is not a tree')
        
        if root is None:
            root = next(iter(nx.topological_sort(G)))
        
        def _hierarchy_pos_recursive(
            G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5, pos=None, parent=None
        ):
            if pos is None:
                pos = {root: (xcenter, vert_loc)}
            else:
                pos[root] = (xcenter, vert_loc)
            
            children = list(G.neighbors(root))
            if not isinstance(G, nx.DiGraph) and parent is not None:
                children.remove(parent)
            
            if len(children) != 0:
                dx = width / len(children)
                nextx = xcenter - width / 2 - dx / 2
                for child in children:
                    nextx += dx
                    pos = _hierarchy_pos_recursive(
                        G, child, width=dx, vert_gap=vert_gap,
                        vert_loc=vert_loc - vert_gap, xcenter=nextx,
                        pos=pos, parent=root
                    )
            return pos
        
        return _hierarchy_pos_recursive(G, root, width, vert_gap, vert_loc, xcenter)


# Convenience functions

def visualize_its(
    its_graph: nx.Graph,
    title: str = "ITS Graph",
    figsize: Tuple[int, int] = (10, 8),
    **kwargs
) -> plt.Figure:
    """
    Quick function to visualize an ITS graph.
    
    Parameters:
    -----------
    its_graph : nx.Graph
        ITS graph to visualize
    title : str
        Graph title
    figsize : tuple
        Figure size
    **kwargs
        Additional arguments for ITSVisualizer.visualize_its_graph
    
    Returns:
    --------
    fig : matplotlib.Figure
    """
    visualizer = ITSVisualizer()
    return visualizer.visualize_its_graph(its_graph, title=title, figsize=figsize, **kwargs)


def visualize_reaction(
    reactants: nx.Graph,
    products: nx.Graph,
    its_graph: nx.Graph,
    figsize: Tuple[int, int] = (18, 6),
    **kwargs
) -> plt.Figure:
    """
    Quick function to visualize a complete reaction.
    
    Parameters:
    -----------
    reactants : nx.Graph
        Reactants graph
    products : nx.Graph
        Products graph
    its_graph : nx.Graph
        ITS graph
    figsize : tuple
        Figure size
    **kwargs
        Additional arguments
    
    Returns:
    --------
    fig : matplotlib.Figure
    """
    visualizer = ITSVisualizer()
    return visualizer.visualize_reaction_triple(
        reactants, products, its_graph, figsize=figsize, **kwargs
    )

