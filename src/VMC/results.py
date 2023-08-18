"""Module defining the `Results` class.

Each `Protein` object has to have its data loaded into in its own `Results` attribute (stored as `Protein.results`).
Data must be loaded into `Protein.results` before atom visualisation can be run.

"""

from tqdm import tqdm
from typing import Any, Dict, Optional
import pickle
import csv

class Results():

    """The `Results` class is used to load and store:
        1. the results of Markov Stability analysis
        2. the atom data
        3. the bond data
        of the `Protein` object.

    Attributes
    ----------
    pdb_code : str
        The PDB code of the `Protein` object.

    matrix_type : str {'A'}
        The type of matrix on which Markov Stability computations were performed:
            'A' - adjacency matrix

    constructor : str {'linearized',
                       'continuous_combinatorial',
                       'continuous_normalised',
                       'signed_modularity',
                       'signed_combinatorial',
                       'directed',
                       'linearized_directed'}
        The constructor used to perform Markov Stability computations within PyGenStability [1]_.

    datetime : str
        The datetime specifying the Markov Stability computation.

    Markov_results : Optional[Dict[str, Any]]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    bond_data : Optional[Dict[int, Dict]]
        The bond_data of the `Protein` object.
        Typically stored in `Protein.results.bond_data`.

    atom_data : Optional[Dict[int, Dict]]
        The atom_data of the `Protein` object.
        Typically stored in `Protein.results.atom_data`.

    """
    def __init__(self,
                 pdb_code: str,
                 matrix_type: str,
                 constructor: str,
                 datetime: str,
                 noncovalent: bool) -> None:

        """Initialisation of the `Results` object.

        Stores the parameters used in the Markov Stability analysis.
        Loads the protein's Markov Stability results.
        Loads the protein's atom_data.
        Loads the protein's bond_data and formats it to facilitate colouring.

        """
        # store the parameters used in MS analysis
        self.pdb_code = pdb_code
        self.matrix_type = matrix_type
        self.constructor = constructor
        self.datetime = datetime
        self.noncovalent = noncovalent

        # store data used to colour bonds
        self.Markov_results = self.load_Markov_results()
        if self.Markov_results is None:
            return
        self.bond_data = self.load_bond_data()
        if self.bond_data is None:
            return
        self.atom_data = self.load_atom_data()
        if self.atom_data is None:
            return

    def load_Markov_results(self) -> Optional[Dict[str, Any]]:

        """Unpickles and returns the results of the Markov Stability analysis.

        Returns
        -------
        Optional[Dict[str, Any]]
            The results of the Markov Stability analysis.
            If the specified results.pkl file cannot be found, returns None.

        """
        # NOTE: requires specific file naming convention
        if self.noncovalent:
            results_dir = f"./pygenstability/results_{self.pdb_code}_{self.matrix_type}_noncovalent_{self.constructor}_{self.datetime}.pkl"
        else:
            results_dir = f"./pygenstability/results_{self.pdb_code}_{self.matrix_type}_{self.constructor}_{self.datetime}.pkl"

        try:
            with open(results_dir, "rb") as results_file:
                print("Markov Stability results loaded successfully")
                return pickle.load(results_file)

        except FileNotFoundError:
            print(f"Markov results file: '{results_dir}' cannot be found")
            print("Please load results again")
            return None

    def load_bond_data(self) -> Optional[Dict[int, Dict]]:

        """Reads, formats, and returns the bond data of the `Protein` object.

        Returns
        -------
        Optional[Dict[int, Dict]]
            The bond data of the protein.
            If the bonds.csv file cannot be found, returns None.

        """
        # NOTE: requires specific file naming convention
        bonds_dir = f"./bagpype/{self.pdb_code}_bonds.csv"

        try:
            with open(bonds_dir, 'r') as bonds_file:

                # count total bonds for tqdm progress bar
                total_bonds = line_counter(bonds_dir)

                bond_data_raw = {}
                bond_reader = csv.DictReader(bonds_file)

                # if only noncovalent bonds are to be coloured
                if self.noncovalent:
                    noncovalent_bond_id = 0
                    for row in tqdm(bond_reader,
                                    desc = "Loading noncovalent bond data",
                                    total = total_bonds,
                                    unit = "bonds"):
                        if ("COVALENT" not in row["bond_type"]) and ("DISULFIDE" not in row["bond_type"]):
                            bond_data_raw[noncovalent_bond_id] = row
                            noncovalent_bond_id += 1
                # if all bonds are to be coloured
                else:
                    for row in tqdm(bond_reader,
                                    desc = "Loading all bond data",
                                    total = total_bonds,
                                    unit = "bonds"):
                        bond_data_raw[int(row["bond_id"])] = row

                format_bond_data(bond_data_raw)
                print("Bond data loaded successfully")
                return bond_data_raw

        except FileNotFoundError:
            print(f"Bond data file: '{bonds_dir}' cannot be found")
            print("Please load results again")
            return None

    def load_atom_data(self):

        """Reads, formats, and returns the atom data of the `Protein` object.

        Returns
        -------
        Optional[Dict[int, Dict]]
            The atom data of the protein.
            If the atoms.csv file cannot be found, returns None.

        """

        # NOTE: requires specific file naming convention
        atoms_dir = f"./bagpype/{self.pdb_code}_atoms.csv"

        try:
            with open(atoms_dir, 'r') as atoms_file:

                # count total atoms for tqdm progress bar
                total_atoms = line_counter(atoms_dir)

                atom_data_raw = {}
                atom_reader = csv.DictReader(atoms_file)

                for row in tqdm(atom_reader,
                                desc = "Loading all atom data",
                                total = total_atoms,
                                unit = "atoms"):
                    atom_data_raw[int(row["id"])] = row

                print("Atom data loaded successfully")
                return atom_data_raw

        except FileNotFoundError:
            print(f"Atom data file: '{atoms_dir}' cannot be found")
            print("Please load results again")
            return None

def format_bond_data(bond_data_raw: Dict[int, Dict]) -> None:

    """Formats the bond data of the protein into a format that facilitates colouring.

    Implemented within the `load_bond_data()` method.

    Parameters
    ----------
    bond_data_raw : Dict[int, Dict]
       Unformatted bond data loaded from the bonds.csv file

    """
    for bond in tqdm(bond_data_raw.values(),
                     desc = "Formatting bond data",
                     unit = "bonds"):

        bond["atom1_id"] = int(bond["atom1_id"]) + 1 # +1 to match PDBnum
        bond["atom2_id"] = int(bond["atom2_id"]) + 1 # +1 to match PDBnum

def line_counter(csv_dir: str) -> int:

    """Counts the total number of lines in the .csv file, excluding the header.

    Used to calculate the total of `tqdm` progress bars in the `load_bond_data()` and `load_atom_data` methods.

    Parameters
    ----------
    csv_dir : str
        The relative directory of the .csv file.

    Returns
    -------
    int
        The total number of lines in the .csv file, excluding the header.

    """
    with open(csv_dir, "r") as csv_file:
        return (len(csv_file.readlines()) - 1) # -1 to account for header

"""
Notes
-----
.. [1]
@article{pygenstability,
         author = {Arnaudon, Alexis and Schindler, Dominik J. and Peach, Robert L. and Gosztolai, Adam and Hodges, Maxwell and Schaub, Michael T. and Barahona, Mauricio},
         title = {PyGenStability: Multiscale community detection with generalized Markov Stability},
         publisher = {arXiv},
         year = {2023},
         doi = {10.48550/ARXIV.2303.05385},
         url = {https://arxiv.org/abs/2303.05385}
}
"""
