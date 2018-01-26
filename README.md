
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/SimmonsBI/bmotif.svg?branch=master)](https://travis-ci.org/SimmonsBI/bmotif) [![codecov](https://codecov.io/gh/SimmonsBI/bmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/SimmonsBI/bmotif)

Overview
--------

bmotif counts occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. bmotif was intended for use in ecology but its methods are general and can be applied to any bipartite graph.

Installation
------------

To install the released version from CRAN:

``` r
install.packages("bmotif")
```

To install the development version from GitHub:

``` r
install.packages("devtools") # install the devtools package
devtools::install_github("SimmonsBI/bmotif") # install bmotif
```

Use
---

bmotif considers all 44 unique bipartite motifs up to six nodes. Within these motifs there are 148 unique positions.

`bmotif` has two functions: `mcount` and `positions`. `mcount` counts occurrences of motifs in a bipartite network. `positions` counts the number of times each node in a network occurs in each of the positions within the motifs.
