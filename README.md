[![PyPI version](https://badge.fury.io/py/VisualiseMarkovCommunities.svg)](https://badge.fury.io/py/VisualiseMarkovCommunities)
[![License](https://img.shields.io/github/license/Ryan-Reese/VMC)](./LICENSE)

# VMC - *Visualise Markov Communities*

This *Python* package allows users to easily colour communities of atoms within proteins using the molecular visualization system *PyMOL* [1]. 
These multiscale communities are found by performing computational analysis using Markov Stability on the energy-weighted atomistic graph representation of the protein. 
For background information on this method of graph partitioning and examples of its use in the unsupervised identification of protein substructures please see [2-5] and [6-10] respectively.

<img src="/assets/header_1.png" alt="header_1" width="400"> <img src="/assets/header_2.png" alt="header_2" width="400">

## Usage

This package was written so that results from the packages [*Bagpype*](https://github.com/FlorianSong/BagPype) and [*PyGenStability*](https://github.com/barahona-research-group/PyGenStability/tree/master) can be directly fed in.

## Functions

## Demo Run

## Contributors
 - Ryan Reese, Yaliraki Group, Department of Chemistry, Imperial College London
   - GitHub: [`Ryan-Reese`](https://github.com/Ryan-Reese)

## References
[1] Schrodinger, LLC. The PyMOL Molecular Graphics System, Version 2.0 (2017).

[2] Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. Stability of graph communities across time scales. Proceedings of the National Academy of Sciences 107, 12755-12760 (2010). https://doi.org:10.1073/pnas.0903215107

[3] Schaub, M. T., Delvenne, J.-C., Yaliraki, S. N. & Barahona, M. Markov Dynamics as a Zooming Lens for Multiscale Community Detection: Non Clique-Like Communities and the Field-of-View Limit. PLOS ONE 7, e32210 (2012). https://doi.org:10.1371/journal.pone.0032210

[4] Delvenne, J.-C., Schaub, M. T., Yaliraki, S. N. & Barahona, M. The stability of a graph partition: A dynamics-based framework for community detection. Dynamics On and Of Complex Networks 2, 221-242 (2013). https://doi.org:10.48550/arXiv.1308.1605

[5] Lambiotte, R., Delvenne, J. C. & Barahona, M. Random Walks, Markov Processes and the Multiscale Modular Organization of Complex Networks. IEEE Transactions on Network Science and Engineering 1, 76-90 (2014). https://doi.org:10.1109/TNSE.2015.2391998

[6] Delmotte, A., Tate, E. W., Yaliraki, S. N. & Barahona, M. Protein multi-scale organization through graph partitioning and robustness analysis: application to the myosin-myosin light chain interaction. Phys Biol 8, 055010 (2011). https://doi.org:10.1088/1478-3975/8/5/055010

[7] Amor, B., Yaliraki, S. N., Woscholski, R. & Barahona, M. Uncovering allosteric pathways in caspase-1 using Markov transient analysis and multiscale community detection. Mol Biosyst 10, 2247-2258 (2014). https://doi.org:10.1039/c4mb00088a

[8] Bacik, K. A., Schaub, M. T., Beguerisse-Díaz, M., Billeh, Y. N. & Barahona, M. Flow-Based Network Analysis of the Caenorhabditis elegans Connectome. PLOS Computational Biology 12, e1005055 (2016). https://doi.org:10.1371/journal.pcbi.1005055

[9] Zhang, H., Salazar, J. D. & Yaliraki, S. N. Proteins across scales through graph partitioning: application to the major peanut allergen Ara h 1. Journal of Complex Networks 6, 679-692 (2017). https://doi.org:10.1093/comnet/cnx052

[10] Peach, R. L. et al. Unsupervised Graph-Based Learning Predicts Mutations That Alter Protein Dynamics. bioRxiv (2019). https://doi.org:10.1101/847426 

[11] Song, F., Barahona, M. & Yaliraki, S. N. BagPype: A Python package for the construction of atomistic, energy-weighted graphs from biomolecular structures.  (2022). https://doi.org:10.5281/zenodo.6326081

[12] Arnaudon, A. et al. barahona-research-group/PyGenStability.  (2023). https://doi.org:10.5281/zenodo.7898442

