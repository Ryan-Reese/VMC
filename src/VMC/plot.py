"""Module containing the functions used to plot the results of Markov Stability analysis using MatPlotLib.

The main function of the module `create_plot()` is implemented within the `Protein.plot()` method.
The functions whose names begin with underscores are modified from functions found in the `plotting` module of PyGenStability [1]_.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.axes import Axes
from typing import Any, Dict, List, Union

def _get_scales(results: Dict[str, Any],
                scale_axis: bool) -> Union[np.ndarray, List[float]]:

    """Returns the Markov timescales of the Markov Stability analysis.

    Used to determine the time values displayed on the x-axis of the graph.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scale_axis : bool
        If True, display the actual Markov timescales on the x-axis.
        If False, simply number the Markov timescales on the x-axis.

    Returns
    -------
    Union[np.ndarray, List[float]]
        The Markov timescales to be displayed on the x-axis.

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
                      scales: Union[np.ndarray, List]) -> None:

    """Plot the number of Markov communities at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot the number of communities.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the x-axis.

    """
    # plot the number of the communities
    ax.plot(scales,
            results["number_of_communities"],
            linestyle="solid",
            linewidth=1.0)
    # label y-axis
    ax.set_ylabel("# communities")
    # set x/y-axes limits
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["number_of_communities"]) * 1.1))

def _plot_stability(results: Dict[str, Any],
                    ax: Axes,
                    scales: Union[np.ndarray, List]) -> None:

    """Plot the value of Markov Stability at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot Markov Stability.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the x-axis.

    """
    # plot Markov stability
    ax.plot(scales,
            results["stability"],
            linestyle="solid",
            linewidth=1.0)
    # label y-axis
    ax.set_ylabel("Stability")
    # set x/y-axes limits
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["stability"]) * 1.1))

def _plot_NVI(results: Dict[str, Any],
              ax: Axes,
              scales: Union[np.ndarray, List]) -> None:

    """Plot the variation of information at each Markov timescale.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to plot the variation of information.

    scales : Union[np.ndarray, List]
        The Markov timescales displayed on the x-axis.

    """
    # plot variation of information
    ax.plot(scales,
            results["NVI"],
            linestyle="solid",
            linewidth=1.0)
    # label y-axis
    ax.set_ylabel(r"NVI")
    # create horizontal line at NVI=1
    ax.axhline(1, ls="--", lw=0.5, c="C2")
    # set x/y-axes limits
    ax.set(xlim=(scales[0],
                 scales[-1]),
           ylim=(0.0,
                 np.max(results["NVI"]) * 1.1))

def label_x_axis(results: Dict[str, Any],
                 ax: Axes,
                 scale_axis: bool) -> None:

    """Label the x-axis of the graph with the Markov timescales.

    Parameters
    ----------
    results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    ax : Axes
        The axes on which to label the Markov timescales.

    scale_axis : bool
        If True, display the actual Markov timescales on the x-axis.
        If False, simply number the Markov timescales on the x-axis.

    """
    if not scale_axis:
        ax.set_xlabel(r"time scale")
    else:
        if results["run_params"]["log_scale"]:
            ax.set_xlabel(r"$log_{10}(t)$")
        else:
            ax.set_xlabel(r"t")

def plot_optimal_scales(results: Dict[str, Any],
                        axes: List[Axes],
                        scales: Union[np.ndarray, List]) -> None:

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
        The Markov timescales displayed on the x-axis.

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
                scale_axis: bool) -> None:

    """Plots the results of the Markov Stability analysis.

    The main function of the module, implemented within the 'Protein.plot()` method.

    Parameters
    ----------
    Markov_results : Dict[str, Any]
        The Markov_results of the `Protein` object.
        Typically stored in `Protein.results.Markov_results`.

    scale_axis : bool
        If True, display the actual Markov timescales on the x-axis.
        If False, simply number the Markov timescales on the x-axis.

    """
    # get Markov time scales
    scales = _get_scales(Markov_results, scale_axis=scale_axis)

    # initialise the figure
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
    _plot_NVI(Markov_results, ax2, scales)

    # label x-axis on subplot 2
    label_x_axis(Markov_results, ax2, scale_axis=scale_axis)

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
