from pymol import cmd
from tqdm import tqdm
from copy import deepcopy
from typing import Any, Optional, List, Dict, Tuple, Union
import pickle
import csv
from random import random
import cv2 as cv
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.axes import Axes
from typing import Any, Dict, List, Union
from tqdm import tqdm

def _get_scales(results: Dict[str, Any], scale_axis: bool = True):
    """Get the scale vector."""

    if not scale_axis:
        return np.arange(len(results["scales"]))

    if results["run_params"]["log_scale"]:
        # returns the base-10 logarithm of the scales array
        return np.log10(results["scales"])

    else:
        return results["scales"]

def _plot_number_comm(results: Dict[str, Any],
                      ax: Axes,
                      scales: Union[np.ndarray, List]):

    """Plot number of communities."""
    ax.plot(scales,
            results["number_of_communities"],
            linestyle="solid",
            linewidth=1.0)

    ax.set_ylabel("# communities")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["number_of_communities"]) * 1.1))

def _plot_stability(results: Dict[str, Any],
                    ax: Axes,
                    scales: Union[np.ndarray, List]):

    """Plot stability."""
    ax.plot(scales,
            results["stability"],
            linestyle="solid",
            linewidth=1.0)

    ax.set_ylabel("Stability")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["stability"]) * 1.1))

def _plot_NVI(results: Dict[str, Any],
              ax: Axes,
              scales: Union[np.ndarray, List]):

    """Plot variation information."""
    ax.plot(scales,
            results["NVI"],
            linestyle="solid",
            linewidth=1.0)

    ax.set_xlabel(r"$log_{10}(t)$")
    ax.set_ylabel(r"NVI")
    ax.axhline(1, ls="--", lw=0.5, c="C2")
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["NVI"]) * 1.1))

def _plot_optimal_scales(results: Dict[str, Any],
                         axes: List[Axes],
                         scales: Union[np.ndarray, List]):

    """Plot optimal Markov scales"""
    for ax in axes:
        for scale in results["selected_partitions"]:
            ax.axvline(scales[scale],
                       linestyle="dashed",
                       linewidth=0.5)
            if ax == axes[0]: # label number of communities
                n_comms = results["number_of_communities"][scale]
                annotation = f"scale={scale}"+"\n"+f"communities={n_comms}"
                ax.annotate(text=annotation,
                            xy=(scales[scale], n_comms),
                            fontsize="small")

def create_plot(Markov_results: Dict[str, Any]):

    # get Markov time scales
    scales = _get_scales(Markov_results, scale_axis=True)

    # create figure
    fig = plt.figure()
    gs = gridspec.GridSpec(3,1, wspace=0.0, hspace=0.1)
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
    _plot_NVI(Markov_results, ax2, scales)

    # plot optimal scales and label communities
    _plot_optimal_scales(Markov_results, axes, scales)

    # post-processing options
    ax0.set_xticks([])
    ax1.set_xticks([])
    fig.align_labels()
    # show plot
    plt.show()

def get_images(prefix: Optional[str],
               image_dir: str) -> Optional[List[str]]:

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
            print(f"No images were found in {image_dir} with prefix {prefix}, please try again")
            return None

    return images

def image_scale_number(image_name: str) -> Union[int, str]:
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

def get_dimensions(images: List[str], image_dir: str):
    # assumes that all images are of the same size
    # gets dimensions of the first image in the list
    frame = cv.imread(os.path.join(image_dir, images[0]))
    height, width = frame.shape[0], frame.shape[1]
    return height, width

def get_video_dir(video_name: Optional[str],
                  video_dir: str,
                  images: List[str]) -> str:

    # if video name is specified, simply ensure that it ends in .mp4v
    if video_name != None:
        if not video_name.endswith(".mp4v"):
            video_name += ".mp4"
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
            video_name = f"{image_name}.mp4v"
            if image_name.endswith(".png"):
                video_name = f"{image_name[:-4]}.mp4v"
            else:
                video_name = f"{image_name}.mp4v"
            return f"{video_dir}/{video_name}"

def compile_images(prefix: Optional[str] = None,
                   video_name: Optional[str] = None,
                   fps: int = 2) -> None:

    # specify image and video directories
    image_dir = "./PyMol_images"
    video_dir = "./PyMol_videos"

    # store all image files in a list
    images = get_images(prefix, image_dir)
    # checks that image list is not empty
    if images == None:
        return

    # get dimensions of photos
    height, width = get_dimensions(images, image_dir)

    # specify video options
    video_dir = get_video_dir(video_name, video_dir, images)

    video = cv.VideoWriter(video_dir,
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
    print(f"Animation '{video_name}' has been generated successfully")

def verify_scale_multi(results: Dict[str, Any],
                       scales: Union[List[int], range]) -> Tuple[bool, int]:
    n_scale = results['run_params']['n_scale']
    for scale in tqdm(scales,
                      desc="Verifying scales",
                      unit="scales"):
        if scale not in range(n_scale):
            return (False, n_scale)
    return (True, n_scale)

class Results():

    def __init__(self,
                 pdb_code: str,
                 matrix_type: str,
                 constructor: str,
                 datetime: str,
                 noncovalent: bool) -> None:

        # store the parameters used in MS analysis
        self.pdb_code = pdb_code
        self.matrix_type = matrix_type
        self.constructor = constructor
        self.datetime = datetime
        self.noncovalent = noncovalent

        # store data used to colour bonds
        self.Markov_results = self.load_Markov_results()
        self.bond_data = self.load_bond_data()

    def load_Markov_results(self) -> Optional[Dict[str, Any]]:

        if self.noncovalent:
            results_dir = f"./PyGenStability_results/results_{self.pdb_code}_{self.matrix_type}_noncovalent_{self.constructor}_{self.datetime}.pkl"
        else:
            results_dir = f"./PyGenStability_results/results_{self.pdb_code}_{self.matrix_type}_{self.constructor}_{self.datetime}.pkl"

        try:
            with open(results_dir, "rb") as results_file:
                print("Markov Stability results loaded successfully")
                return pickle.load(results_file)

        except FileNotFoundError:
            print("Markov results file: '{results_dir}' does not exist")
            print("Please try to load results again")
            return None

    def load_bond_data(self) -> Optional[Dict[int, Dict]]:

        bonds_dir = f"./bin/{self.pdb_code}_bonds.csv"

        # TODO: add total for tqdm
        try:
            with open(bonds_dir, 'r') as bonds_file:
                bond_data_raw = {}
                bond_reader = csv.DictReader(bonds_file)

                if self.noncovalent:
                    noncovalent_bond_id = 0
                    for row in tqdm(bond_reader,
                                    desc = "Loading noncovalent bond data",
                                    unit = "bonds"):
                        if ("COVALENT" not in row["bond_type"]) and ("DISULFIDE" not in row["bond_type"]):
                            bond_data_raw[noncovalent_bond_id] = row
                            noncovalent_bond_id += 1
                else:
                    for row in tqdm(bond_reader,
                                    desc = "Loading all bond data",
                                    unit = "bonds"):
                        bond_data_raw[int(row["bond_id"])] = row

                format_bond_data(bond_data_raw)
                print("Bond data loaded successfully")
                return bond_data_raw

        except FileNotFoundError:
            print(f"Bond data file: '{bonds_dir}' does not exist")
            print("Please try to load results again")
            return None

def format_bond_data(bond_data_raw: Dict[int, Dict]) -> None:
    for bond in tqdm(bond_data_raw.values(),
                          desc = "Formatting bond data",
                          total = len(bond_data_raw),
                          unit = "bonds"):
        bond["atom1_id"] = int(bond["atom1_id"]) + 1 # +1 to match id column from PDB file
        bond["atom2_id"] = int(bond["atom2_id"]) + 1 # +1 to match id column from PDB file
        bond["atom1_res"] = (int(bond["atom1_res_num"]), bond["atom1_chain"])
        bond["atom2_res"] = (int(bond["atom2_res_num"]), bond["atom2_chain"])

def verify_scale(results: Dict[str, Any], scale: int) -> Tuple[bool, int]:
    n_scale = results['run_params']['n_scale']
    return (scale in range(n_scale), n_scale)

def delete_previous(total_bonds: int) -> None:
    for i in tqdm(range(total_bonds),
                  desc="Deleting previous communities",
                  unit="communities"):
        cmd.delete(f"community-{i}*")

def allocate_colours(total_communities: int) -> None:
    for i in tqdm(range(total_communities),
                  desc="Allocating colours",
                  unit="colours"):
        cmd.set_color(str(i), [random(), random(), random()])

def colour_bond(bond_community: int,
                bond_id: int,
                atom1_id: int,
                atom2_id: int) -> None:
    # set bond name
    bond_name = f"community-{bond_community}-{bond_id}"
    # create distance object between the two atoms
    cmd.distance(bond_name, f"id {atom1_id}", f"id {atom2_id}")
    # colour according to community
    cmd.color(str(bond_community), bond_name)
    # visual options for bond
    cmd.set("dash_gap", 0, bond_name)
    cmd.set("dash_radius", 0.15, bond_name)

def remove_unspecified(bond_data: Dict[int, Dict],
                       specific_residues: List[Tuple[int, str]],
                       specific_bonds: List[int],
                       specific_communities: List[int],
                       community_matrix: List[int]) -> None:

    # remove unspecified residues
    if len(specific_residues) > 0:
        for bond_id in tqdm(list(bond_data),
                            desc="Removing unspecified residues",
                            unit="bonds"):
            if (bond_data[bond_id]["atom1_res"] not in specific_residues) and (bond_data[bond_id]["atom2_res"] not in specific_residues):
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

class Protein:

    def __init__(self, pdb_code: str) -> None:

        self.pdb_code = pdb_code
        self.results = None

    def load_PDB(self,
                 representation: str = "sticks") -> None:

        cmd.load(f"./PDBs/{self.pdb_code}_stripped_H.pdb")
        cmd.show_as(representation)

    def load_results(self,
                     matrix_type: str,
                     constructor: str,
                     datetime: str,
                     noncovalent: bool) -> None:

        temp = Results(self.pdb_code,
                       matrix_type,
                       constructor,
                       datetime,
                       noncovalent)

        if (temp.Markov_results != None) and (temp.bond_data != None):
            self.results = temp

    def image(self, *params: Any) -> None:

        # construct filename separated by underscores
        file_header = f"./PyMol_images/{self.pdb_code}"
        filename = [file_header]
        for param in params:
            filename.append(str(param))
        filename = "_".join(filename)

        # image options (can be modified)
        cmd.png(filename=filename,
            width=2400,
            height=2400,
            dpi=150,
            ray=1)

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
        assert self.results.Markov_results != None, "Markov results have not been loaded"
        create_plot(self.results.Markov_results)

    def reset(self) -> None:
        cmd.delete("all")

    def visualise(self,
                  scale: int,
                  specific_residues: List[Tuple[int, str]] = [],
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


        # get data from Markov results file
        total_communities = self.results.Markov_results['number_of_communities'][scale]
        community_matrix = self.results.Markov_results['community_id'][scale]
        total_bonds = len(community_matrix)

        # create copy of bond data to remove unspecified
        bond_data = deepcopy(self.results.bond_data)

        if not show_entire_community:

            # delete bonds from previous run
            delete_previous(total_bonds)
            # allocate colours equal to total number of communities
            allocate_colours(total_communities)
            # remove unspecified residues/bonds/communities
            remove_unspecified(bond_data,
                               specific_residues,
                               specific_bonds,
                               specific_communities,
                               community_matrix)

            # colour every bond remaining in bond_data dictionary
            for bond_id in tqdm(bond_data.keys(),
                                desc="Colouring bonds",
                                total=len(bond_data),
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
            remove_unspecified(bond_data,
                               specific_residues,
                               specific_bonds,
                               specific_communities,
                               community_matrix)

            for bond_id in tqdm(bond_data.keys(),
                                desc="Scanning for communities",
                                total=len(bond_data),
                                unit="bonds"):

                # add bond community to set
                bond_community = community_matrix[bond_id]
                communities_to_colour.add(bond_community)

            # colour entire communities
            self.visualise(scale,
                           specific_residues=[],
                           specific_bonds=[],
                           specific_communities=list(communities_to_colour),
                           show_entire_community=False,
                           image=image)

    def visualise_multi(self,
                        scales: Optional[Union[List[int], range]] = None,
                        specific_residues: List[Tuple[int, str]] = [],
                        specific_bonds: List[int] = [],
                        specific_communities: List[int] = [],
                        show_entire_community: bool = False,
                        image: bool = True) -> None:

        # preliminary checks
        if self.results == None:
            print("Results have not been loaded, please load results first")
            return
        assert self.results.Markov_results != None, "No Markov results"
        assert self.results.bond_data != None, "No bond data"

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
            remove_unspecified(scale_bond_data,
                               specific_residues,
                               specific_bonds,
                               specific_communities,
                               scale_community_matrix)

            if not show_entire_community:

                # delete bonds from previous scale
                delete_previous(total_bonds)

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

                self.visualise_multi(scales,
                                     specific_residues=[],
                                     specific_bonds=[],
                                     specific_communities=list(scale_communities_to_colour),
                                     show_entire_community=False,
                                     image=image)

def init_options() -> None:
    cmd.bg_color("grey")
    cmd.set("ray_opaque_background", "on")
    cmd.util.performance(0)

    """
    Disabling undo greatly reduces memory cost.
    cmd.set operations apparently have a very high memory overhead.
    This option will increase visualisation by ~100x.
    """
    cmd.undo_disable()

def main():
    init_options()
    cmd.util.cbag()

main()
