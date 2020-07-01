This is an update to fix the noLD issue, as requested by Brian Ripley.

## Test environments
* Debian Linux, R-devel, GCC, no long double support (rhub)
* local OS X install, R 4.0.2
* ubuntu 16.04.6 (on travis-ci), release and devel
* win-builder, release and devel
* Debian Linux, R-devel, GCC ASAN/UBSAN (rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTES:

* checking CRAN incoming feasibility ... NOTE
The Date field is over a month old.

I have now updated the date field to be current
  
* checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    libs   3.7Mb
  
It seems that on Linux architectures (Ubuntu and Fedora rhub configurations, though not on the rhub Debian configuration), the check returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OS X. My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp. Indeed, some functions of the bmotif package have been written in C++ using Rcpp. They are needed to perform calculation of motif frequencies. Without the speed up gained from those C++ functions, this package would become impractical.

## Downstream dependencies
There are currently no known downstream dependencies for this package.
