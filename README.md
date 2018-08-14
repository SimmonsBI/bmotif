
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/SimmonsBI/bmotif.svg?branch=master)](https://travis-ci.org/SimmonsBI/bmotif)
[![codecov](https://codecov.io/gh/SimmonsBI/bmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/SimmonsBI/bmotif)

## Overview

`bmotif` is software for motif analyses of bipartite networks. It can
count occurrences of motifs in bipartite networks, as well as the number
of times each node or link appears in each unique node or link position
within motifs (a node or link’s structural role). `bmotif` supports
weighted as well as unweighted networks: the mean weight of motifs can
be calculated, as well as the standard deviation of motifs mean weights;
weighted versions of node and link position counts are also supported
(see below). As well as R, core functionality is also available in
[MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and
[Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was
originally developed to analyse bipartite species interaction networks
in ecology but its methods are general and can be applied to any
bipartite graph.

## Installation

To install the released version from CRAN:

``` r
install.packages("bmotif")
```

To install the development version from GitHub:

``` r
install.packages("devtools") # install the devtools package
devtools::install_github("SimmonsBI/bmotif") # install bmotif
```

## Use

`bmotif` considers all 44 unique bipartite motifs up to six nodes.
Within these motifs there are 148 unique node positions and 106 unique
link positions.

`bmotif` has three functions, with can all be used with binary or
weighted networks:

1.  `mcount`
      - **Unweighted**: counts how many times each motif occurs in a
        bipartite network
      - **Weighted**: calculates the mean weight of motifs and the
        standard deviation of their weights
2.  `node_positions`: counts the number of times each node in a network
    occurs in each of the unique node positions within the motifs.
      - **Unweighted**: counts the number of times each node in a
        network occurs in each of the unique node positions within the
        motifs.  
      - **Weighted**: calculates a range of weighted metrics, such as
        the mean link strength of a each node in each position or the
        contribution of each node to each motif’s total weight
3.  `link_positions`:
      - **Unweighted**: counts the number of times each link in a
        network occurs in each of the unique link positions within the
        motifs.
      - **Weighted**: calculates the number of times each link in a
        network occurs in each of the unique link positions within the
        motifs, multiplied by each link’s strength

The motifs corresponding to each motif ID and the node positions
corresponding to each motif node position ID can be found in **Simmons,
B. I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di
Clemente, R. bmotif: a package for counting motifs in bipartite
networks. bioRxiv. doi: 10.1101/302356**. These were defined following
Baker et al (2015) Appendix 1 Figure A27.

## License

The code is released under the MIT license (see LICENSE file).

## Citation

If you use the package in your work, please cite: Simmons, B. I.,
Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di Clemente, R.
bmotif: a package for counting motifs in bipartite networks. bioRxiv.
doi: 10.1101/302356

## References

Baker, N.J., Kaartinen, R., Roslin, T. and Stouffer, D.B., 2015.
Species’ roles in food webs show fidelity across a highly variable oak
forest. Ecography, 38(2), pp.130-139.
