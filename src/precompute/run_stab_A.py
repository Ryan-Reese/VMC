import pygenstability
from pygenstability.plotting import plot_scan
from scipy.sparse import csr_matrix
from scipy.io import loadmat
from datetime import datetime
from typing import Any, Dict, Optional

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
    Documentation can be found in PyGenStability [1]_.
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

"""
Notes
-----
.. [1]
Arnaudon, Alexis, Dominik J. Schindler, Robert L. Peach, Adam Gosztolai, Maxwell Hodges, Michael T. Schaub, and Mauricio Barahona.
"PyGenStability: Multiscale community detection with generalized Markov Stability."
arXiv preprint arXiv:2303.05385 (2023).
"""
