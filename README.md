[![PyPI version](https://badge.fury.io/py/VisualiseMarkovCommunities.svg)](https://badge.fury.io/py/VisualiseMarkovCommunities)
[![License](https://img.shields.io/github/license/Ryan-Reese/VMC)](./LICENSE)

# VMC - *Visualise Markov Communities*
This Python package allows users to easily colour communities of atoms within proteins using the molecular visualization system *PyMOL*[^1]. 
These multiscale communities are found by performing computational analysis using Markov Stability on the energy-weighted atomistic graph representation of the protein. 
For background information on this method of graph partitioning and examples of its use in the unsupervised identification of protein substructures please see [^2][^3][^4][^5] and [^6][^7][^8][^9][^10] respectively.

<img src="/assets/header_1.png" alt="header_1" width="400"> <img src="/assets/header_2.png" alt="header_2" width="400">

## Usage
This package was designed so that results from the packages [*Bagpype*](https://github.com/FlorianSong/BagPype)[^11] and [*PyGenStability*](https://github.com/barahona-research-group/PyGenStability/tree/master)[^12] can be directly fed in. The package also includes [scripts](./src/VMC/precompute_A.py) that simplify the process of obtaining these results. An illustration of the intended scheme is shown below:

<img src="/assets/scheme.png" alt="scheme" width="">


## Installation
The official PyPI version of *VMC* can be installed with

```python
pip install VisualiseMarkovCommunities
```

> [!IMPORTANT]
> Ensure that VMC is installed into *PyMOL*’s bundled version of Python

## Tutorial + Demo Run
### Initialisation
To initialise *VMC* into the current *PyMOL* session, simply use

```python
import VMC
```

A grey background indicates that the package has been successfully imported.

<img src="/assets/demo_init.png" alt="demo_init" width="400">

### Loading a protein
To load a protein's structure from its PDB file, simply create a `Protein` object by specifying its PDB code and use the built-in method `load_PDB()`

```python
from VMC import Protein
prot_name = Protein(pdb_code: str)
prot_name.load_PDB()
```

The name of the `Protein` object does not matter. In this tutorial, it is named `ADK` due to the PDB code belonging to a conformation of Adenylate Kinase in Aquifex Aeolicus[^13].
```python
from VMC import Protein
ADK = Protein("2rh5")
ADK.load_PDB()
```

<img src="/assets/demo_load_protein.png" alt="demo_load_protein" width="400">

> [!NOTE]
> If a local PDB file is used, ensure that the PDB file is saved in the relative directory `./PDBs`.
> Otherwise, *VMC* will raise an error when the specified PDB file cannot be found.

### Loading the results of Markov Stability
Before atoms communities can be visualised, results must be loaded into the `Protein` object. 
This is achieved using the method `load_results()`

```python
prot_name.load_results(matrix_type: str, constructor: str, datetime: str)
```

The currently-supported inputs for this method are

| matrix_type | constructor | datetime |
| :---: | :---: | :---: |
| `A` (energy-weighted adjacency matrix) | `linearized` `continuous_combinatorial` `continuous_normalized` | DDMMYY-hh_mm |
| Constructed using [*Bagpype*](https://github.com/FlorianSong/BagPype)[^11] | Chosen within [*PyGenStability*](https://github.com/barahona-research-group/PyGenStability/tree/master)[^12] | Used to distinguish between multiple results files whose other parameters are identical |

For instance,

```python
ADK.load_results(matrix_type = “A”, constructor = “linearized”, datetime = “040823-19_40”)
```

<img src="/assets/demo_load_results_print.png" alt="demo_load_results_print" width="">

> [!NOTE]
> Results obtained using the included [precomputation scripts](./src/VMC/precompute_A.py) are automatically saved in the relative directory `./pygenstability`.
> By default, *VMC* searches for the results file in this folder when attempting to load results, hence using the included scripts is recommended.

### Visualising Atom Communities
We finally arrive at the crux of the package!

To colour and visualise atoms by their community, simply use the method `visualise_A()`.
The only parameter that **has** to be specified is the Markov timescale that would like to be visualised. Optional parameters for `visualise_A()` are specified in the [`protein`](src/VMC/protein.py) module.

```python
ADK.visualise_A(scale = 90)
```

<img src="/assets/demo_visualise_A.png" alt="demo_visualise_A" width="400"> <img src="/assets/demo_visualise_A_prot.png" alt="demo_visualise_A_prot" width="400">

To visualise atom communities at multiple Markov timescales, the method `visualise_multi_A()` can be used

```python
ADK.visualise_multi_A(scales = range(20))
```

If `image = True` then snapshots of the atom communities at each Markov timescale are named and saved in .png format within the relative directory `./pymol_images`.

Snapshots saved in `./pymol_images` can then be compiled into an animation using the method `compile()`, which is saved in the relative directory `./pymol_videos`.
This allows us to visualise how atom communities change as the Markov timescale increases. For instance,

```python
ADK.visualise_multi_A(scales = [0,25,50,75], image = True)
ADK.compile()
```

<img src="/assets/demo_visualise_A_directory.png" alt="demo_visualise_A_directory" width="400"> <img src="/assets/demo_compile_directory.png" alt="demo_compile_directory" width="400"> 

Importantly, the results of Markov Stability analysis can be plotted and viewed without exiting *PyMOL* using the method `plot()`

```python
ADK.plot()
```

<img src="/assets/demo_plot.png" alt="demo_plot" width="400"> <img src="/assets/plot.png" alt="plot" width="400"> 

The most robust partitions as found by [*PyGenStability*](https://github.com/barahona-research-group/PyGenStability/tree/master)[^12] are automatically labelled with their scale number and number of communities.

## Examples


https://github.com/Ryan-Reese/VMC/assets/109569773/e6d1757f-9a0c-4988-95e3-8344365bc9b1


## Functions

Currently built-in methods for the `Protein` class are:
| Method | Function |
| :---: | :--- |
| `load_PDB()` | Loads the structure of the protein into PyMol using its PDB file |
| `load_results()` | Loads the Markov Stability results into the `Protein` object |
| `image()` | Saves a PNG image of the protein into the `./pymol_images/` directory |
| `compile()` | Compiles the images in the `./pymol_images/` directory into an MP4V video |
| `plot()` | Plots the Markov Stability results using MatPlotLib |
| `visualise_A()` | Colour atoms based on their communities at a specific Markov timescale |
| `visualise_multi_A()` | Colour atoms based on their communities at multiple Markov timescales |

## Contributors
 - Ryan Reese, Yaliraki Group, Department of Chemistry, Imperial College London
   - GitHub: [`Ryan-Reese`](https://github.com/Ryan-Reese)

## *References*
[^1]: Schrodinger, LLC. The PyMOL Molecular Graphics System, Version 2.0 (2017).

[^2]: Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. Stability of graph communities across time scales. Proceedings of the National Academy of Sciences 107, 12755-12760 (2010). https://doi.org:10.1073/pnas.0903215107
  
[^3]: Schaub, M. T., Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. Markov Dynamics as a Zooming Lens for Multiscale Community Detection: Non Clique-Like Communities and the Field-of-View Limit. PLOS ONE 7, e32210 (2012). https://doi.org:10.1371/journal.pone.0032210
  
[^4]: Delvenne, J.-C., Schaub, M. T., Yaliraki, S. N. & Barahona, M. The stability of a graph partition: A dynamics-based framework for community detection. Dynamics On and Of Complex Networks 2, 221-242 (2013). https://doi.org:10.48550/arXiv.1308.1605
  
[^5]: Lambiotte, R., Delvenne, J. C. & Barahona, M. Random Walks, Markov Processes and the Multiscale Modular Organization of Complex Networks. IEEE Transactions on Network Science and Engineering 1, 76-90 (2014). https://doi.org:10.1109/TNSE.2015.2391998

[^6]: Delmotte, A., Tate, E. W., Yaliraki, S. N. & Barahona, M. Protein multi-scale organization through graph partitioning and robustness analysis: application to the myosin-myosin light chain interaction. Phys Biol 8, 055010 (2011). https://doi.org:10.1088/1478-3975/8/5/055010

[^7]: Amor, B., Yaliraki, S. N., Woscholski, R. & Barahona, M. Uncovering allosteric pathways in caspase-1 using Markov transient analysis and multiscale community detection. Mol Biosyst 10, 2247-2258 (2014). https://doi.org:10.1039/c4mb00088a

[^8]: Bacik, K. A., Schaub, M. T., Beguerisse-Díaz, M., Billeh, Y. N. & Barahona, M. Flow-Based Network Analysis of the Caenorhabditis elegans Connectome. PLOS Computational Biology 12, e1005055 (2016). https://doi.org:10.1371/journal.pcbi.1005055

[^9]: Zhang, H., Salazar, J. D. & Yaliraki, S. N. Proteins across scales through graph partitioning: application to the major peanut allergen Ara h 1. Journal of Complex Networks 6, 679-692 (2017). https://doi.org:10.1093/comnet/cnx052

[^10]: Peach, R. L. et al. Unsupervised Graph-Based Learning Predicts Mutations That Alter Protein Dynamics. bioRxiv (2019). https://doi.org:10.1101/847426 

[^11]: Song, F., Barahona, M. & Yaliraki, S. N. BagPype: A Python package for the construction of atomistic, energy-weighted graphs from biomolecular structures.  (2022). https://doi.org:10.5281/zenodo.6326081

[^12]: Arnaudon, A. et al. barahona-research-group/PyGenStability.  (2023). https://doi.org:10.5281/zenodo.7898442

[^13]: Henzler-Wildman, K. A. et al. Intrinsic motions along an enzymatic reaction trajectory. Nature 450, 838-844 (2007). https://doi.org:10.1038/nature06410
