import bagpype
import networkx as nx
import pygenstability
from datetime import datetime
from pygenstability.plotting import plot_scan
from scipy.io import loadmat, savemat
from scipy.sparse import csr_matrix
from typing import Any, Dict, Optional

def parse_pdb_to_prot(pdb_code: str) -> bagpype.molecules.Protein:

    # create `Protein` object in BagPype
    myprot = bagpype.molecules.Protein()

    # specify relative path to pdb file
    # if download=False, then pdb file must already be stored locally
    pdb_file = f'./PDBs/{pdb_code}.pdb'
    parser = bagpype.parsing.PDBParser(pdb_file, download=True)

    """
    Parse information from the PDB file into `Protein` object.
    All parsing options can be modified.
    Documentation can be found in BagPype [1]_.
    """

    parser.parse(myprot,
                 model=1,
                 strip={"res_name": ["HOH"], "chain": ["B", "C", "D", "E"]},
                 strip_HETATM=True,
                 trim_H="all",
                 add_H=True,
                 MakeMultimer_number=None)

    return myprot

def construct_atomistic_graph(myprot: bagpype.molecules.Protein,
                              pdb_code: str) -> None:

    # specify relative paths to save atoms/bonds csv files
    atoms_file = f'./bagpype/{pdb_code}_atoms.csv'
    bonds_file = f'./bagpype/{pdb_code}_bonds.csv'

    ggenerator = bagpype.construction.Graph_constructor()
    ggenerator.construct_graph(myprot,
                               atoms_file_name=atoms_file,
                               bonds_file_name=bonds_file)

def verify_connectness(graph: nx.Graph) -> bool:
    if nx.number_connected_components(graph) > 1:
        return False
    return True

def calculate_A(myprot: bagpype.molecules.Protein) -> csr_matrix:

    # passing nodelist=atoms ensures that the row/column order corresponds atom ids
    atoms = [atom.id for atom in myprot.atoms]
    return nx.adjacency_matrix(myprot.graph, nodelist=atoms)

def save_A(A: csr_matrix,
           pdb_code: str) -> None:

    savemat(f'./adjacency_matrices/A_{pdb_code}.mat', {'A': A})

def create_A(pdb_code: str) -> None:

    # parse PDB file to Bagpype protein object
    myprot = parse_pdb_to_prot(pdb_code)

    # construct atomistic graph of protein using Bagpype
    construct_atomistic_graph(myprot, pdb_code)

    # ensure that the input (adjacency) matrix is connected
    if not verify_connectness(myprot.graph):
        print("""The atomistic graph has more than 1 connected component.
              This is incompatible with perturbation calculations.""")
        return

    # calculate A
    A = calculate_A(myprot)

    # save A as .mat file according to PDB code
    save_A(A, pdb_code)

def load_A(pdb_code: str) -> Optional[csr_matrix]:

    mat_dir = f'./adjacency_matrices/A_{pdb_code}.mat'

    try:
        mat_file = loadmat(mat_dir)
        A = mat_file["A"]

    except FileNotFoundError:
        print(f"Adjacency matrix file '{mat_dir}' does not exist")
        return None

    return A

def get_datetime() -> str:

    return datetime.now().strftime("%d%m%y-%H:%M")

def get_constructor() -> Optional[str]:

    constructors = [constructor.partition("_")[2] for constructor
                    in dir(pygenstability.constructors)
                    if constructor.startswith('constructor')]

    constructor = input("Input constructor:\nRecommended: linearized, continuous_combinatorial")

    # if the constructor is not recognised
    if constructor not in constructors:
        print(f"Constructor {constructor} not recognised.")
        print(f"Constructor must be one of: {constructors}")
        return None

    return constructor

def get_results(A: csr_matrix,
                pdb_code: str,
                constructor: str,
                datetime: str) -> Dict[str, Any]:

    """
    Compute graph clustering across scales with Markov Stability.
    All run parameters can be modified.
    Documentation can be found in PyGenStability [2]_.
    """

    return pygenstability.run(graph=A,
                              constructor=constructor,
                              min_scale=-2.0,
                              max_scale=2.0,
                              n_scale=100,
                              log_scale=True,
                              scales=None,
                              n_tries=100,
                              with_NVI=True,
                              n_NVI=100,
                              with_postprocessing=True,
                              with_ttprime=True,
                              with_spectral_gap=False,
                              exp_comp_mode="spectral",
                              result_file=f"./pygenstability/results_{pdb_code}_A_{constructor}_{datetime}.pkl",
                              n_workers=20,
                              tqdm_disable=False,
                              with_optimal_scales=True,
                              optimal_scales_kwargs=None,
                              method="leiden",
                              constructor_kwargs=None)

def plot_results(results: Dict[str, Any],
                 pdb_code: str,
                 constructor: str,
                 datetime: str) -> None:

    plot_scan(results,
              scale_axis=True,
              figure_name=f'./PyGenStability/results_{pdb_code}_A_{constructor}_{datetime}.pdf',
              use_plotly=False,
              live=True,
              plotly_filename=f'./PyGenStability/results_{pdb_code}_A_{constructor}_{datetime}.html')

def run_stab_A(pdb_code: str) -> None:

    # load adjacency matrix from .mat file
    A = load_A(pdb_code)
    # if the adjacency matrix cannot be found
    if A is None:
        return

    # input constructor
    constructor = get_constructor()
    # if the constructor is not recognised
    if constructor is None:
        return

    # get datetime
    datetime = get_datetime()

    # run Markov stability computations
    results = get_results(A, pdb_code, constructor, datetime)

    # save plot of results as pdf file
    plot_results(results, pdb_code, constructor, datetime)

def main():

    # take PDB code from user input
    pdb_code = input("PDB Code: ")

    # create A matrix
    create_A(pdb_code)

    # run Markov Stability
    run_stab_A(pdb_code)

main()

"""
Notes
-----
.. [1]
Florian Song, Mauricio Barahona & Sophia N. Yaliraki.
BagPype: A Python package for the construction of atomistic, energy-weighted graphs from biomolecular structures.
(2022).

.. [2]
Arnaudon, Alexis, Dominik J. Schindler, Robert L. Peach, Adam Gosztolai, Maxwell Hodges, Michael T. Schaub, and Mauricio Barahona.
"PyGenStability: Multiscale community detection with generalized Markov Stability."
arXiv preprint arXiv:2303.05385 (2023).
"""
