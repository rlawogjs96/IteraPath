import warnings

warnings.warn(
    "The 'SynVis' subpackage is deprecated and will be removed in future releases. "
    "Please migrate to the 'synkit' package as soon as possible,"
    + " which offers enhanced functionality. "
    "You can install it directly using pip: `pip install synkit`.",
    FutureWarning,
)

# Import visualization classes
from .its_visualizer import ITSVisualizer, visualize_its, visualize_reaction
from .chemical_graph_visualizer import ChemicalGraphVisualizer
from .chemical_reaction_visualizer import ChemicalReactionVisualizer

__all__ = [
    "ITSVisualizer",
    "visualize_its",
    "visualize_reaction",
    "ChemicalGraphVisualizer",
    "ChemicalReactionVisualizer",
]