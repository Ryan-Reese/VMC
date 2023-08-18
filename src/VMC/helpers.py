"""Miscellaneous helper functions used to assist methods in the `Protein` class.

"""

from tqdm import tqdm
from typing import Any, Dict, List, Tuple, Union

def verify_scale(results: Dict[str, Any], scale: int) -> Tuple[bool, int]:

    """Verifies that the Markov time scale to be visualised is present in the Markov results.

    Implemented within the `Protein.visualise_A/M()` methods.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scale : int
        The Markov time scale whose communities are to be visualised.

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

    Implemented within the `Protein.visualise_multi_A/M` methods.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scales : Union[List[int], range]
        The list or range of Markov time scales whose communities are to be visualised.

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

