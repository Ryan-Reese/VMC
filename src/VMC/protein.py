"""The main module of VMC, defining the `Protein` class.

The standard usage of VMC always begins with initialising a `Protein` object.
The `Protein` object allows the user to access the operations within VMC (executed through methods of the `Protein` object).

The currently supported operations are:
    Protein.load_PDB()
    Protein.load_results()
    Protein.image()
    Protein.compile()
    Protein.plot()
    Protein.visualise_A()
    Protein.visualise_multi_A()

"""

# External packages
from pymol import cmd
from tqdm import tqdm
from copy import deepcopy
from typing import Any, List, Optional, Union

# Internal modules
from VMC.compile import compile_images
from VMC.cleaning import remove_unspecified_A
from VMC.colouring import allocate_colours, colour_atom, delete_previous_A
from VMC.helpers import verify_scale, verify_scale_multi
from VMC.plot import create_plot
from VMC.results import Results

class Protein:

    """This main object of VMC, used to access and run the package's operations.

    Attributes
    ----------
    pdb_code : str
        The 4-character PDB code of the protein.

    results : Results
        The `Results` attribute of the `Protein` object.
        Data must be loaded into `Protein.results` before atom visualisation can be run.

    """
    def __init__(self, pdb_code: str) -> None:

        """Initialisation of the `Protein` object.

        Stores the PDB code of the protein and initialises the `Protein.results` attribute to None.

        """
        self.pdb_code = pdb_code
        self.results = None

    def load_PDB(self,
                 representation: str = "sticks",
                 pdb_file_suffix: str = "_stripped_H.pdb") -> None:

        """Loads the structure of the protein into PyMol using its PDB file.

        Typically the first operation to be run on the `Protein` object.

        Parameters
        ----------
        representation : str {'sticks', 'lines', 'cartoon', 'ribbon'}
            The PyMol representation that determines how PyMol draws the protein.
            For adjacency matrices, 'sticks' or 'lines' is recommended.

        pdb_file_suffix : str
            The suffix of the PDB file loaded into PyMol.
            Determined by the parsing options used in BagPype [1]_.

        """
        cmd.load(f"./PDBs/{self.pdb_code}" + pdb_file_suffix)
        cmd.show_as(representation)

    def load_results(self,
                     matrix_type: str,
                     constructor: str,
                     datetime: str,
                     noncovalent: bool = False) -> None:

        """Attempts to load the data of the protein into the `Protein.results` attribute.

         `Protein.results` is only reassigned if Markov Stability analysis, atom data, and bond data are all successfully loaded.
         Details of input parameters can be found in the `results.py` module.

        """
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

        """Saves a .png image of the protein within PyMol.

        Images are saved by default in the ./pymol_images/ directory.
        Imaging options can be modified to your liking within the `cmd.png` API command.

        Parameters
        ----------
        *params : Tuple[Any]
            Additional parameters to be appended to the filename of the .png image.

        """
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
                dpi=200,
                ray=1)
        print(f"Image saved as {filename}.png")

    def compile(self,
                prefix: Optional[str] = None,
                video_name: Optional[str] = None,
                fps: int = 2) -> None:

        """Compiles images in the ./pymol_images/ directory into a .mp4v video file.

        Videos are saved by default in the ./pymol_videos/ directory.
        Details of input parameters and other options can be found in the `compile.py` module.

        """
        if prefix == None:
            compile_images(prefix=self.pdb_code,
                           video_name=video_name,
                           fps=fps)
        else:
            compile_images(prefix=prefix,
                           video_name=video_name,
                           fps=fps)

    def plot(self) -> None:

        """Plots the results of the Markov Stability analysis using MatPlotLib.

        Details can be found in the `plotting.py` module.

        """
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        create_plot(self.results.Markov_results, scale_axis=True)

    def visualise_A(self,
                    scale: int,
                    specific_residues: List[int] = [],
                    specific_atoms: List[int] = [],
                    specific_communities: List[int] = [],
                    show_entire_community: bool = False,
                    image: bool = False) -> None:

        """The main operation of VMC, colouring atoms based on their communities at a specific Markov timescale.

        This operation is only to be used when Markov Stability analysis has been run on adjacency (A) matrices.
        If communities at multiple Markov timescales would like to be visualised, please use the `visualise_multi_A()` method instead.

        Parameters
        ----------
        scale : int
            The Markov timescale whose atom communities are to be visualised.

        specific_residues : List[int]
            The specific residues to be coloured, specified by res_num.
            An empty list corresponds to colouring all residues.

        specific_atoms : List[int]
            The specific atoms to be coloured, specified by atom_id.
            An empty list corresponds to colouring all atoms.

        specific_communities : List[int]
            The specific communities to be coloured, specified by community number.
            An empty list corresponds to colouring all communities.

        show_entire_community : bool
            If True, entire communities will be coloured if at least one residue/atom belonging to the community is specified.

        image : bool
            If True, will save a .png image of the protein within PyMol following community visualisation.

        """
        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None
        assert self.results.atom_data != None

        scale_check = verify_scale(self.results.Markov_results, scale)
        if not scale_check[0]:
            print(f"Markov timescale is not in range {scale_check[1]}, please provide a valid scale")
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

        """Colour atoms based on their communities at multiple Markov timescales.

        This operation is only to be used when Markov Stability analysis has been run on adjacency (A) matrices.
        If atom communities at a singular Markov timescale would like to be visualised, please use the `visualise_A()` method instead.

        Parameters
        ----------
        scales : Optional[Union[List[int], range]]
            The Markov timescales whose atom communities are to be visualised.
            Scales be provided as either a list of integers or a range object.

        Please refer to the `visualise_A()` method for details of the remaining parameters.

        """
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
                print(f"At least one of the Markov timescales is not in range {scale_check[1]}, please provide a valid range of scales")
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

"""
Notes
-----
.. [1]
Florian Song, Mauricio Barahona & Sophia N. Yaliraki.
BagPype: A Python package for the construction of atomistic, energy-weighted graphs from biomolecular structures.
(2022).
"""
