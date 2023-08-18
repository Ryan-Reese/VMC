import os
import cv2 as cv
from tqdm import tqdm
from typing import List, Optional, Tuple, Union
from pymol import cmd
from random import random
from typing import Any, Dict, List, Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.axes import Axes
import pickle
import csv
from copy import deepcopy

class Protein:

    def __init__(self, pdb_code: str) -> None:

        self.pdb_code = pdb_code
        self.results = None

    def load_PDB(self,
                 representation: str = "sticks",
                 pdb_file_suffix: str = "_stripped_H.pdb") -> None:

        # PDB filename is determined by BagPype parsing options
        cmd.load(f"./PDBs/{self.pdb_code}" + pdb_file_suffix)
        cmd.show_as(representation)

    def load_results(self,
                     matrix_type: str,
                     constructor: str,
                     datetime: str,
                     noncovalent: bool = False) -> None:

        temp = Results(self.pdb_code,
                       matrix_type,
                       constructor,
                       datetime,
                       noncovalent)

        if (temp.Markov_results == None):
            return
        if (temp.bond_data == None):
            return
        if (temp.atom_data == None):
            return

        self.results = temp

    def image(self, *params: Any) -> None:

        # construct filename separated by underscores
        file_header = f"./pymol_images/{self.pdb_code}"
        filename = [file_header]
        for param in params:
            filename.append(str(param))
        filename = "_".join(filename)

        # centre on specific Protein object and capture entire object
        cmd.zoom(complete=1)

        # image options (can be modified to your liking)
        cmd.png(filename=filename,
                width=2400,
                height=2400,
                dpi=150,
                ray=1)
        print(f"Image saved as {filename}.png")

    def compile(self,
                prefix: Optional[str] = None,
                video_name: Optional[str] = None,
                fps: int = 2) -> None:

        if prefix == None:
            compile_images(prefix=self.pdb_code,
                           video_name=video_name,
                           fps=fps)
        else:
            compile_images(prefix=prefix,
                           video_name=video_name,
                           fps=fps)

    def plot(self) -> None:
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        create_plot(self.results.Markov_results, scale_axis=True)

    def reset(self) -> None:
        cmd.delete("all")

    def visualise_A(self,
                    scale: int,
                    specific_residues: List[int] = [],
                    specific_atoms: List[int] = [],
                    specific_communities: List[int] = [],
                    show_entire_community: bool = False,
                    image: bool = False) -> None:

        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        assert self.results.atom_data != None

        scale_check = verify_scale(self.results.Markov_results, scale)
        if not scale_check[0]:
            print(f"Time scale is not in range {scale_check[1]}, please provide a valid scale")
            return

        # get data from Markov results file
        total_communities = self.results.Markov_results['number_of_communities'][scale]
        community_matrix = self.results.Markov_results['community_id'][scale]

        # create copy of bond data to remove unspecified
        atom_data = deepcopy(self.results.atom_data)

        if not show_entire_community:

            # delete atom colouring from previous run
            delete_previous_A()
            # allocate colours equal to total number of communities
            allocate_colours(total_communities)
            # remove unspecified residues/atoms/communities
            remove_unspecified_A(atom_data,
                                 specific_residues,
                                 specific_atoms,
                                 specific_communities,
                                 community_matrix)

            # colour every atom remaining in the atom_data dictionary
            for atom_id in tqdm(atom_data.keys(),
                                desc="Colouring atoms",
                                unit="atoms"):

                atom_community = community_matrix[atom_id]
                atom_PDBnum = atom_data[atom_id]["PDBnum"]

                colour_atom(atom_community, atom_PDBnum)

            # save as png image
            if image:
                self.image(self.results.matrix_type,
                           self.results.constructor,
                           self.results.datetime,
                           f"scale_{scale}")

        else:
            # create set of communities to be coloured
            communities_to_colour = set()

            # remove unspecified residues/atoms/communities
            remove_unspecified_A(atom_data,
                                 specific_residues,
                                 specific_atoms,
                                 specific_communities,
                                 community_matrix)

            for atom_id in tqdm(atom_data.keys(),
                                desc="Scanning for communities",
                                unit="atoms"):

                # add bond community to set
                atom_community = community_matrix[atom_id]
                communities_to_colour.add(atom_community)

            # colour entire communities
            self.visualise_A(scale,
                             specific_residues=[],
                             specific_atoms=[],
                             specific_communities=list(communities_to_colour),
                             show_entire_community=False,
                             image=image)

    def visualise_multi_A(self,
                          scales: Optional[Union[List[int], range]] = None,
                          specific_residues: List[int] = [],
                          specific_atoms: List[int] = [],
                          specific_communities: List[int] = [],
                          show_entire_community: bool = False,
                          image: bool = True) -> None:

        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        assert self.results.atom_data != None

        if scales == None:
            scales = range(self.results.Markov_results['run_params']['n_scale'])
            print("Visualising all scales")
        else:
            scale_check = verify_scale_multi(self.results.Markov_results, scales)
            if not scale_check[0]:
                print(f"At least one of the time scales is not in range {scale_check[1]}, please provide a valid range of scales")
                return
            print("Time scales have been verified")

        # get data from Markov results file
        max_communities = max(self.results.Markov_results["number_of_communities"])

        # allocate colours equal to the maximum number of communities
        allocate_colours(max_communities)

        for scale in tqdm(scales,
                          desc="Capturing scales",
                          unit="scales"):

            # get community matrix for the specific scale
            scale_community_matrix = self.results.Markov_results['community_id'][scale]
            # create copy of atom data to remove unspecified
            scale_atom_data = deepcopy(self.results.atom_data)
            # remove unspecified residues/atoms/communities
            remove_unspecified_A(scale_atom_data,
                                 specific_residues,
                                 specific_atoms,
                                 specific_communities,
                                 scale_community_matrix)

            if not show_entire_community:

                # delete atom colouring from previous scale
                delete_previous_A()

                # colour every atom remaining in atom_data dictionary
                for atom_id in tqdm(scale_atom_data.keys(),
                                    desc="Colouring atoms",
                                    unit="atoms"):

                    atom_community = scale_community_matrix[atom_id]
                    atom_PDBnum = scale_atom_data[atom_id]["PDBnum"]

                    colour_atom(atom_community, atom_PDBnum)

                # save as png image
                if image:
                    self.image(self.results.matrix_type,
                               self.results.constructor,
                               self.results.datetime,
                               f"scale_{scale}")

            else:
                # create set of communities to be coloured
                scale_communities_to_colour = set()

                for atom_id in tqdm(scale_atom_data.keys(),
                                    desc="Scanning for communities",
                                    unit="atoms"):

                    # add atom community to set
                    atom_community = scale_community_matrix[atom_id]
                    scale_communities_to_colour.add(atom_community)

                self.visualise_A(scale,
                                 specific_residues=[],
                                 specific_atoms=[],
                                 specific_communities=list(scale_communities_to_colour),
                                 show_entire_community=False,
                                 image=image)

    def visualise_M(self,
                    scale: int,
                    # TODO: change specific residue input to only resi_num
                    specific_residues: List[int] = [],
                    specific_bonds: List[int] = [],
                    specific_communities: List[int] = [],
                    show_entire_community: bool = False,
                    image: bool = False) -> None:

        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        assert self.results.bond_data != None

        scale_check = verify_scale(self.results.Markov_results, scale)
        if not scale_check[0]:
            print(f"Time scale is not in range {scale_check[1]}, please provide a valid scale")
            return
        print("Time scale has been verified")

        # get data from Markov results file
        total_communities = self.results.Markov_results['number_of_communities'][scale]
        community_matrix = self.results.Markov_results['community_id'][scale]
        total_bonds = len(community_matrix)

        # create copy of bond data to remove unspecified
        bond_data = deepcopy(self.results.bond_data)

        if not show_entire_community:

            # delete bonds from previous run
            delete_previous_M(total_bonds)
            # allocate colours equal to total number of communities
            allocate_colours(total_communities)
            # remove unspecified residues/bonds/communities
            remove_unspecified_M(bond_data,
                                 specific_residues,
                                 specific_bonds,
                                 specific_communities,
                                 community_matrix)

            # colour every bond remaining in bond_data dictionary
            for bond_id in tqdm(bond_data.keys(),
                                desc="Colouring bonds",
                                unit="bonds"):

                bond_community = community_matrix[bond_id]
                atom1_id = bond_data[bond_id]["atom1_id"]
                atom2_id = bond_data[bond_id]["atom2_id"]

                colour_bond(bond_community, bond_id, atom1_id, atom2_id)

            # hide distance labels in PyMol
            cmd.hide('labels')

            # save as png image
            if image:
                if self.results.noncovalent:
                    self.image(self.results.matrix_type,
                               "noncovalent",
                               self.results.constructor,
                               self.results.datetime,
                               f"scale_{scale}")
                else:
                    self.image(self.results.matrix_type,
                               self.results.constructor,
                               self.results.datetime,
                               f"scale_{scale}")

        else:
            # create set of communities to be coloured
            communities_to_colour = set()

            # remove unspecified residues/bonds/communities
            remove_unspecified_M(bond_data,
                               specific_residues,
                               specific_bonds,
                               specific_communities,
                               community_matrix)

            for bond_id in tqdm(bond_data.keys(),
                                desc="Scanning for communities",
                                unit="bonds"):

                # add bond community to set
                bond_community = community_matrix[bond_id]
                communities_to_colour.add(bond_community)

            # colour entire communities
            self.visualise_M(scale,
                             specific_residues=[],
                             specific_bonds=[],
                             specific_communities=list(communities_to_colour),
                             show_entire_community=False,
                             image=image)

    def visualise_multi_M(self,
                          scales: Optional[Union[List[int], range]] = None,
                          specific_residues: List[int] = [],
                          specific_bonds: List[int] = [],
                          specific_communities: List[int] = [],
                          show_entire_community: bool = False,
                          image: bool = True) -> None:

        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        assert self.results.bond_data != None

        if scales == None:
            scales = range(self.results.Markov_results['run_params']['n_scale'])
            print("Visualising all scales")
        else:
            scale_check = verify_scale_multi(self.results.Markov_results, scales)
            if not scale_check[0]:
                print(f"At least one of the time scales is not in range {scale_check[1]}, please provide a valid range of scales")
                return
            print("Time scales have been verified")

        # get data from Markov results file
        total_bonds = len(self.results.Markov_results["community_id"][0])
        max_communities = max(self.results.Markov_results["number_of_communities"])

        # allocate colours equal to the maximum number of communities
        allocate_colours(max_communities)

        for scale in tqdm(scales,
                          desc="Capturing scales",
                          unit="scales"):

            # get community matrix for the specific scale
            scale_community_matrix = self.results.Markov_results['community_id'][scale]
            # create copy of bond data to remove unspecified
            scale_bond_data = deepcopy(self.results.bond_data)
            # remove unspecified residues/bonds/communities
            remove_unspecified_M(scale_bond_data,
                                 specific_residues,
                                 specific_bonds,
                                 specific_communities,
                                 scale_community_matrix)

            if not show_entire_community:

                # delete bonds from previous scale
                delete_previous_M(total_bonds)

                # colour every bond remaining in bond_data dictionary
                for bond_id in tqdm(scale_bond_data.keys(),
                                    desc="Colouring bonds",
                                    unit="bonds"):

                    bond_community = scale_community_matrix[bond_id]
                    atom1_id = scale_bond_data[bond_id]["atom1_id"]
                    atom2_id = scale_bond_data[bond_id]["atom2_id"]

                    colour_bond(bond_community, bond_id, atom1_id, atom2_id)

                # hide distance labels in PyMol
                cmd.hide('labels')

                # save as png image
                if image:
                    if self.results.noncovalent:
                        self.image(self.results.matrix_type,
                                   "noncovalent",
                                   self.results.constructor,
                                   self.results.datetime,
                                   f"scale_{scale}")
                    else:
                        self.image(self.results.matrix_type,
                                   self.results.constructor,
                                   self.results.datetime,
                                   f"scale_{scale}")
            else:
                # create set of communities to be coloured
                scale_communities_to_colour = set()

                for bond_id in tqdm(scale_bond_data.keys(),
                                    desc="Scanning for communities",
                                    unit="bonds"):

                    # add bond community to set
                    bond_community = scale_community_matrix[bond_id]
                    scale_communities_to_colour.add(bond_community)

                self.visualise_M(scale,
                                 specific_residues=[],
                                 specific_bonds=[],
                                 specific_communities=list(scale_communities_to_colour),
                                 show_entire_community=False,
                                 image=image)

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

    matrix_type : str, {'A', 'M'}
        The type of matrix on which Markov Stability computations were performed:
            'A' - adjacency matrix
            'M' - bond-to-bond propensity matrix

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

    noncovalent : bool
        Only relevant if matrix_type == 'M'.
        If True, indicates that Markov Stability computations were only performed on noncovalent bonds of the protein.
        If False, indicates that Markov Stability computations were perfored on all bonds of the protein.

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

        """Initialisation of the `Results` class.

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
def _get_scales(results: Dict[str, Any],
                scale_axis: bool = True) -> Union[np.ndarray, List[float]]:

    """Returns the Markov timescales of the Markov Stability analysis.

    Used to determine the values displayed on the scale axis of the graph.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scale_axis : bool
        If True, display the actual Markov timescales on the scale axis.
        If False, simply number the Markov timescales on the scale axis.

    Returns
    -------
    Union[np.ndarray, List[float]]
        The Markov timescales displayed on the scale axis.

    """
    if not scale_axis:
        return np.arange(len(results["scales"]))
    else:
        if results["run_params"]["log_scale"]:
            # returns the base-10 logarithm of the scales array
            return np.log10(results["scales"])
        else:
            return results["scales"]

def _plot_number_comm(results: Dict[str, Any],
                      ax: Axes,
                      scales: Union[np.ndarray, List]):

    """Plot the number of Markov communities at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot the number of communities.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the scale axis.

    """
    # plot the number of the communities
    ax.plot(scales,
            results["number_of_communities"],
            linestyle="solid",
            linewidth=1.0)
    # label and set x/y-axes
    ax.set_ylabel("# communities")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["number_of_communities"]) * 1.1))

def _plot_stability(results: Dict[str, Any],
                    ax: Axes,
                    scales: Union[np.ndarray, List]):

    """Plot the value of Markov Stability at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot Markov Stability.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the scale axis.

    """
    # plot Markov stability
    ax.plot(scales,
            results["stability"],
            linestyle="solid",
            linewidth=1.0)
    # label and set x/y-axes
    ax.set_ylabel("Stability")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["stability"]) * 1.1))

def _plot_NVI(results: Dict[str, Any],
              ax: Axes,
              scales: Union[np.ndarray, List],
              scale_axis: bool = True):

    """Plot the variation of information at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot the variation of information.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the scale axis.

    """
    # plot variation of information
    ax.plot(scales,
            results["NVI"],
            linestyle="solid",
            linewidth=1.0)
    # label and set x/y-axes
    if not scale_axis:
        ax.set_xlabel(r"time scale")
    else:
        if results["run_params"]["log_scale"]:
            ax.set_xlabel(r"$log_{10}(t)$")
        else:
            ax.set_xlabel(r"t")
    ax.set_ylabel(r"NVI")
    ax.axhline(1, ls="--", lw=0.5, c="C2")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["NVI"]) * 1.1))

def plot_optimal_scales(results: Dict[str, Any],
                        axes: List[Axes],
                        scales: Union[np.ndarray, List]):

    """Plot the optimal Markov timescales found during the Markov Stability analysis.

    Determination of optimal Markov timescales is detailed in PyGenStability [1]_.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    axes : Axes
        The axes where the number of communities, Markov Stability, and NVI have been plot.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the scale axis.

    """
    optimal_scales = [int(scale) for scale in results["selected_partitions"]]

    for ax in axes:
        for scale in optimal_scales:
            ax.axvline(scales[scale],
                       linestyle="dashed",
                       linewidth=0.5)
            # label scale number and number of communities in subplot 0
            if ax == axes[0]:
                n_comms = results["number_of_communities"][scale]
                annotation = f"scale={scale}"+"\n"+f"communities={n_comms}"
                ax.annotate(text=annotation,
                            xy=(scales[scale], n_comms),
                            fontsize="xx-small")

def create_plot(Markov_results: Dict[str, Any],
                scale_axis: bool = True):

    """Creates a plot of the Markov Stability results.

    The main function of the module, implemented within the 'Protein.plot()` method.

    Parameters
    ----------
    Markov_results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    """
    # get Markov time scales
    scales = _get_scales(Markov_results, scale_axis=scale_axis)

    # initialisation
    fig = plt.figure()
    gs = gridspec.GridSpec(3,1, wspace=0.0, hspace=0.1)  # border space can be changed
    axes = []

    # plot stability in subplot 0
    ax0 = fig.add_subplot(gs[0])
    axes.append(ax0)
    _plot_number_comm(Markov_results, ax0, scales)

    # plot stability in subplot 1
    ax1 = fig.add_subplot(gs[1])
    axes.append(ax1)
    _plot_stability(Markov_results, ax1, scales)

    # plot NVI in subplot 2
    ax2 = fig.add_subplot(gs[2])
    axes.append(ax2)
    _plot_NVI(Markov_results, ax2, scales, scale_axis=scale_axis)

    # plot optimal scales and label communities
    plot_optimal_scales(Markov_results, axes, scales)

    # post-processing options
    ax0.set_xticks([])
    ax1.set_xticks([])
    fig.align_labels()

    # show plot
    plt.show()

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
def init_options() -> None:

    # sets background colour to grey
    # assists in visualising the community colours
    cmd.bg_color("grey")  # Can be changed

    # gives ray-traced images an opaque (default: grey) background
    cmd.set("ray_opaque_background", "on")

    # puts PyMol in maximum-quality mode
    cmd.util.performance(0)

    # Disabling undo greatly reduces memory cost.
    # PyMol `cmd.set` operations apparently have a very high memory overhead.
    # From rough tests, this option will increase the rate of atom/bond colouring by ~100x.
    cmd.undo_disable()
def verify_scale(results: Dict[str, Any], scale: int) -> Tuple[bool, int]:

    """Verifies that the Markov time scale to be visualised is present in the Markov results.

    Assists the `Protein.visualise_A/M()` methods.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scale : int
        The Markov time scale to be visualised.

    Returns
    -------
    Tuple[bool, int]
        A tuple containing:
            - A boolean indicating whether the Markov time scale to be visualised is within `results`
            - An integer indicating the total number of Markov time scales in `results`.
    """
    n_scale = results['run_params']['n_scale']
    return (scale in range(n_scale), n_scale)

def verify_scale_multi(results: Dict[str, Any],
                       scales: Union[List[int], range]) -> Tuple[bool, int]:

    """Verifies that the Markov time scales to be visualised are present in the Markov results.

    Assists the `Protein.visualise_A/M_multi()` methods.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scales : Union[List[int], range]
        The list or range of Markov time scales to be visualised.

    Returns
    -------
    Tuple[bool, int]
        A tuple containing:
            - A boolean indicating whether the Markov time scales to be visualised are all within `results`
            - An integer indicating the total number of Markov time scales in `results`.

    """
    n_scale = results['run_params']['n_scale']
    for scale in tqdm(scales,
                      desc="Verifying scales",
                      unit="scales"):
        if scale not in range(n_scale):
            return (False, n_scale)
    return (True, n_scale)

def delete_previous_A() -> None:

    """Resets all atoms colours to grey.

    Used to reset atom colours between Markov timescales.
    """
    cmd.color("grey")

def delete_previous_M(total_bonds: int) -> None:

    """Deletes all previously coloured bonds generated for another Markov timescale.

    Used to reset bond colours between Markov timescales.

    Parameters
    ----------
    total_bonds : int
        The total number of bonds in the `Protein` object's bond data, corresponding to the maximum number of possible Markov communities.

    """
    for i in tqdm(range(total_bonds),
                  desc="Deleting previous communities",
                  unit="communities"):
        cmd.delete(f"community-{i}*")

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

def colour_bond(bond_community: int,
                bond_id: int,
                atom1_id: int,
                atom2_id: int) -> None:

    """Colours a bond between two atoms according to the bond's Markov community.

    Visual options for the bond can be modified by altering the PyMol API commands in this function.

    Parameters
    ----------
    bond_community : int
        The Markov community that the bond belongs to.

    bond_id : int
        The bond's bond_id.

    atom1_id : int
        The first atom in the bond's atom_id.

    atom2_id : int
        The second atom in the bond's atom_id.

    """
    # set bond name
    bond_name = f"community-{bond_community}-{bond_id}"
    # create distance object between the two atoms
    cmd.distance(bond_name, f"id {atom1_id}", f"id {atom2_id}")
    # colour the distance object according to the Markov community
    cmd.color(str(bond_community), bond_name)

    # visual options for bond
    # the bottom two API commands can be modified
    cmd.set("dash_gap", 0, bond_name)
    cmd.set("dash_radius", 0.15, bond_name)

def remove_unspecified_A(atom_data: Dict[int, Dict],
                         specific_residues: List[int],
                         specific_atoms: List[int],
                         specific_communities: List[int],
                         community_matrix: List[int]) -> None:

    """Remove all unspecified atoms from atom data.

    This function deletes all unspecified residues, atoms, and communities from the Protein's atom data.
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

def remove_unspecified_M(bond_data: Dict[int, Dict],
                         specific_residues: List[int],
                         specific_bonds: List[int],
                         specific_communities: List[int],
                         community_matrix: List[int]) -> None:

    """Remove all unspecified bonds from bond data.

    This function deletes all unspecified residues, bonds, and communities from the Protein's bond data.
    Only the bonds that remain will proceed to be coloured according to their communities.
    This function is implemented within the `Protein.visualise_M()` and `Protein.visualise_M_multi()` methods.

    Parameters
    ----------
    bond_data : Dict[int, Dict]
        The bond_data of the `Protein` object.
        Typically stored as `Protein.results.bond_data`.

    specific_residues : List[int]
        The specific residues to be coloured, specified by res_num.

    specific_bonds : List[int]
        The specific bonds to be coloured, specified by bond_id.

    specific_communities : List[int]
        The specific communities to be coloured, specified by community number.

    community_matrix : List[int]
        The community matrix for the specific Markov timescale.
        Used to determine which community number each bond belongs to.

    """
    # remove unspecified residues
    if len(specific_residues) > 0:
        for bond_id in tqdm(list(bond_data),
                            desc="Removing unspecified residues",
                            unit="bonds"):
            if (int(bond_data[bond_id]["atom1_res_num"]) not in specific_residues) and (int(bond_data[bond_id]["atom2_res_num"]) not in specific_residues):
                del bond_data[bond_id]

    # remove unspecified bonds
    if len(specific_bonds) > 0:
        for bond_id in tqdm(list(bond_data),
                            desc="Removing unspecified bonds",
                            unit="bonds"):
            if bond_id not in specific_bonds:
                del bond_data[bond_id]

    # remove unspecified communities
    if len(specific_communities) > 0:
        for bond_id in tqdm(list(bond_data),
                            desc="Removing unspecified communities",
                            unit="bonds"):
            if community_matrix[bond_id] not in specific_communities:
                del bond_data[bond_id]

def get_images(prefix: Optional[str],
               image_dir: str) -> Optional[List[str]]:

    """Returns the list of images to be compiled.

    Parameters
    ----------
    prefix : Optional[str]
        If none, all images in `image_dir` are returned.
        Otherwise, only images beginning with `prefix` are returned.

    image_dir : str
       The relative directory containing the images to be compiled.

    Returns
    -------
    Optional[List[str]]
        If no images are found in `image_dir` matching the (optional) `prefix`, None is returned.
        Otherwise, the list of image names to be compiled is returned.

    """
    # sort images in list by their scale number
    if prefix == None:
        images = sorted([image for image in os.listdir(image_dir)],
                        key=image_scale_number)

        if len(images) == 0:
            print(f"No images were found in {image_dir}, please try again")
            return None

    else:
        images = sorted([image for image in os.listdir(image_dir)
                         if image.startswith(prefix)],
                        key=image_scale_number)

        if len(images) == 0:
            print(f"No images were found in {image_dir} with 'prefix' {prefix}, please try again")
            return None

    return images

def image_scale_number(image_name: str) -> Union[int, str]:

    """Returns the scale number of the image from its name (if present).

    Used to sort images by their scale number when retrieving image names in `get_images()`.

    Parameters
    ----------
    image_name : str
        The name of the image.

    Returns
    -------
    Union[int, str]
        If the image contains a scale number, return the scale number as an integer.
        Otherwise, return the original image name as a string.

    """
    params = image_name.split("_")

    # tries to find scale number in image name
    try:
        scale_location = params.index("scale")
        if params[scale_location + 1].endswith(".png"):
            return int(params[scale_location + 1][:-4])
        else:
            return int(params[scale_location + 1])

    # if no scale number is found, returns image name
    except ValueError:
        return image_name

def get_dimensions(images: List[str], image_dir: str) -> Tuple[int, int]:

    """Returns the height and width of the image.

    Used to determine the dimensions of the final video.
    Only the first image in the list of images is passed to this function,
    as all images are assumed to have the same dimensions.

    Parameters
    ----------
    images : List[str]
        List of images to be compiled.

    image_dir : str
        The relative directory containing the images to be compiled.

    Returns
    -------
    Tuple[int, int]
        The height and width of the image.

    """
    frame = cv.imread(os.path.join(image_dir, images[0]))
    height, width = frame.shape[0], frame.shape[1]
    return height, width

def get_video_dir(video_name: Optional[str],
                  video_dir: str,
                  images: List[str]) -> str:

    """Returns the relative directory and name of the video to be saved.

    If `video_name` is specified, the specified name will be used.
    Otherwise, the video name is based off the name of the first image in the list of images.

    Parameters
    ----------
    video_name : Optional[str]
        Optional custom name for the video.

    video_dir : str
        Relative directory in which to save the video.

    images : List[str]
        List of images to be compiled.

    Returns
    -------
    str
        Relative directory and name of the video.

    """
    # if video name is specified, simply ensure that it ends in .mp4v
    if video_name != None:
        if not video_name.endswith(".mp4v"):
            video_name += ".mp4v"
        return f"{video_dir}/{video_name}"

    else:
        # if no video name is specified, base the video name off the first image name
        image_name = images[0]

        # tries to remove scale parameter from image name if present
        try:
            scale_location = image_name.index("scale")
            video_name = f"{image_name[:scale_location].rstrip('_')}.mp4v"
            return f"{video_dir}/{video_name}"

        # if no scale parameter is found, returns image name
        except ValueError:
            if image_name.endswith(".png"):
                video_name = f"{image_name[:-4]}.mp4v"
            else:
                video_name = f"{image_name}.mp4v"
            return f"{video_dir}/{video_name}"

def compile_images(prefix: Optional[str] = None,
                   video_name: Optional[str] = None,
                   fps: int = 2) -> None:

    """Compiles images into a .mp4v video file.

    The main function within this module, implemented within the `Protein.compile()` method.

    Parameters
    ----------
    prefix : Optional[str]
        If none, all images in `image_dir` are compiled.
        Otherwise, only images beginning with `prefix` are compiled.

    video_name : Optional[str]
        If specified, the video will be saved with this name.
        Otherwise, the video name is based on the name of the first image in the list of images.

    fps : int
        Frames per second of the video.
        A higher value corresponds to a faster playback speed of the video.

    """
    # specify relative image and video directories
    ##############################################
    image_dir = "./pymol_images"
    video_dir = "./pymol_videos"
    ##############################################

    # store all image files in a list
    images = get_images(prefix, image_dir)

    # checks that image list is not empty
    if images == None:
        return

    # get dimensions of photos
    height, width = get_dimensions(images, image_dir)

    # specify video name and location
    full_video_dir = get_video_dir(video_name, video_dir, images)

    # specify video options
    video = cv.VideoWriter(full_video_dir,
                           fourcc=0,
                           fps=fps,
                           frameSize=(width,height),
                           isColor=True)

    # write images to video
    for image_file in tqdm(images,
                           desc="Writing images to video",
                           unit="images"):

        image = cv.imread(os.path.join(image_dir, image_file))
        cv.putText(img=image,
                   text=image_file,
                   org=(100,2300),  # may have to be adjusted if cmd.png options are changed
                   fontFace=cv.FONT_HERSHEY_SIMPLEX,  # can be modified
                   fontScale=2,  # can be modified
                   color=(255,255,255),  # can be modified
                   thickness=3,  # can be modified
                   lineType=cv.LINE_AA,  # can be modified
                   bottomLeftOrigin=False)
        video.write(image)

    video.release()
    cv.destroyAllWindows()
    print(f"Animation '{full_video_dir}' has been generated successfully")

init_options()
