import bagpype
import networkx as nx
from scipy.io import savemat
from scipy.sparse import csr_matrix

def parse_pdb_to_prot(pdb_code: str) -> bagpype.molecules.Protein:

    # create `Protein` object in BagPype
    myprot = bagpype.molecules.Protein()

    # specify relative path to pdb file
    # if download=False, then pdb file must already be stored locally
    pdb_file = f'./PDBs/{pdb_code}.pdb'
    parser = bagpype.parsing.PDBParser(pdb_file, download=False)

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

"""
Notes
-----
.. [1]
Florian Song, Mauricio Barahona & Sophia N. Yaliraki.
BagPype: A Python package for the construction of atomistic, energy-weighted graphs from biomolecular structures.
(2022).
"""
