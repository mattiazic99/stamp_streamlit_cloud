"""Integration with the STRING protein-protein interaction network.

Computes average inverse distance from non-switching genes to the
cumulative set of switching genes per age bracket.
"""
import networkx as nx
import pandas as pd


def load_string_network(path: str, min_score: int = 900) -> nx.Graph:
    """Load STRING PPI links above a confidence threshold."""
    raise NotImplementedError("TODO")


def average_inverse_distance(
    graph: nx.Graph,
    cumulative_switching_genes: list[str],
) -> float:
    """Mean inverse shortest-path distance from other genes to the set."""
    raise NotImplementedError("TODO")


def ppi_heatmap_data(
    graph: nx.Graph,
    sets_by_tissue: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    """Compute the (tissue x age) heatmap of average inverse distances."""
    raise NotImplementedError("TODO")
