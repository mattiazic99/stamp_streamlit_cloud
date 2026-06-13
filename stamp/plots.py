"""Plot generation: bar charts, heatmaps, network diagrams."""
import pandas as pd


def bar_switching_per_age(
    sets_by_bracket: dict[str, list[str]], output_path: str
) -> None:
    """Bar chart of switching gene counts across age brackets."""
    raise NotImplementedError("TODO")


def heatmap_jaccard(similarity: pd.DataFrame, output_path: str) -> None:
    """Render the Jaccard tissue-similarity heatmap."""
    raise NotImplementedError("TODO")


def heatmap_ppi(ppi_data: pd.DataFrame, output_path: str) -> None:
    """Render the PPI average-inverse-distance heatmap."""
    raise NotImplementedError("TODO")


def network_migration(migrations: pd.DataFrame, output_path: str) -> None:
    """Render the tissue-to-tissue migration network."""
    raise NotImplementedError("TODO")
