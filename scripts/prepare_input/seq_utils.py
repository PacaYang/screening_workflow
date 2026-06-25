###################################
# Shared helpers for parsing protein sequence cells
###################################
"""
Utilities shared by the structure-prediction input generators.

A target's ``sequence`` cell in ``sequences.csv`` may hold several protein
chains as a single comma-separated string (the B2_gen2/MCP convention), e.g.::

    MAGAL...ASA,TVFFW...AKA

``parse_chains`` turns that cell into an ordered list of per-chain sequences.
"""

import string


def parse_chains(cell):
    """Split a comma-separated sequence cell into a list of chain sequences.

    Whitespace around each chain is stripped and empty fragments are dropped,
    so a single-chain cell returns a one-element list.

    Args:
        cell: the raw ``sequence`` value (str).

    Returns:
        list[str]: ordered, non-empty chain sequences.
    """
    if cell is None:
        return []
    return [chunk.strip() for chunk in str(cell).split(",") if chunk.strip()]


def chain_ids(n):
    """Return ``n`` protein chain ids: A, B, C, ... (A-Z, then AA, AB, ...)."""
    letters = string.ascii_uppercase
    ids = []
    for i in range(n):
        if i < len(letters):
            ids.append(letters[i])
        else:
            # Two-letter fallback for very large complexes (rare).
            ids.append(letters[i // len(letters) - 1] + letters[i % len(letters)])
    return ids


def ligand_ids(n):
    """Return ``n`` ligand chain ids starting at Z and counting down: Z, Y, X...

    The screened compound is always id ``Z`` so downstream scorers that key on
    chain ``Z`` keep working regardless of how many ligands are co-folded.
    """
    letters = string.ascii_uppercase  # ...XYZ
    # Z, Y, X, W, ... by walking the alphabet backwards.
    return [letters[len(letters) - 1 - i] for i in range(n)]
