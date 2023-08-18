"""Module containing all functions used to colour atoms and bonds.

Atoms/Bonds are coloured according to their Markov communities at each Markov timescale.
The colouring is achieved using PyMol's API commands.
This functions in this module are implemented within the `Protein.visualise_A()` and `Protein.visualise_A()` methods.

"""

from pymol import cmd
from tqdm import tqdm
from random import random

def delete_previous_A() -> None:

    """Resets all atoms colours to grey.

    Used to reset atom colours between Markov timescales.
    """
    cmd.color("grey")

def allocate_colours(total_communities: int) -> None:

    """Generates a random RGB colour for each Markov community number.

    Parameters
    ----------
    total_communities : int
        The total number of Markov communities present in the specific Markov timescale.

    """
    for i in tqdm(range(total_communities),
                  desc="Allocating colours",
                  unit="colours"):
        cmd.set_color(str(i), [random(), random(), random()])

def colour_atom(atom_community: int,
                atom_PDBnum: int) -> None:

    """Colours an atom according to the atom's Markov community.

    Parameters
    ----------
    atom_community : int
        The Markov community that the atom belongs to.

    atom_PDBnum : int
        The atom's PDB number, equivalent to the atom's (atom_id + 1).

    """
    # colour according to community
    cmd.color(str(atom_community), f"id {atom_PDBnum}")

