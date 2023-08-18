"""Module used to remove all unspecified items from data during community visualisation.

"""

from tqdm import tqdm
from typing import Dict, List

def remove_unspecified_A(atom_data: Dict[int, Dict],
                         specific_residues: List[int],
                         specific_atoms: List[int],
                         specific_communities: List[int],
                         community_matrix: List[int]) -> None:

    """Remove all unspecified atoms from atom data.

    This function deletes all unspecified residues, atoms, and communities from the `Protein` object's atom data.
    Only the atoms that remain will proceed to be coloured according to their communities.
    This function is implemented within the `Protein.visualise_A()` and `Protein.visualise_A_multi()` methods.

    Parameters
    ----------
    atom_data : Dict[int, Dict]
        The atom_data of the `Protein` object.
        Typically stored as `Protein.results.atom_data`.

    specific_residues : List[int]
        The specific residues to be coloured, specified by res_num.

    specific_atoms : List[int]
        The specific atoms to be coloured, specified by atom_id.

    specific_communities : List[int]
        The specific communities to be coloured, specified by community number.

    community_matrix : List[int]
        The community matrix for the specific Markov timescale.
        Used to determine which community number each atom belongs to.

    """
    # remove unspecified residues
    if len(specific_residues) > 0:
        for atom_id in tqdm(list(atom_data),
                            desc="Removing unspecified residues",
                            unit="atoms"):
            if (int(atom_data[atom_id]["res_num"]) not in specific_residues):
                del atom_data[atom_id]

    # remove unspecified atoms
    if len(specific_atoms) > 0:
        for atom_id in tqdm(list(atom_data),
                            desc="Removing unspecified atoms",
                            unit="atoms"):
            if atom_id not in specific_atoms:
                del atom_data[atom_id]

    # remove unspecified communities
    if len(specific_communities) > 0:
        for atom_id in tqdm(list(atom_data),
                            desc="Removing unspecified communities",
                            unit="atoms"):
            if community_matrix[atom_id] not in specific_communities:
                del atom_data[atom_id]

