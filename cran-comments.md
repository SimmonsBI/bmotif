## Test environments
* local OS X install, R 3.6.1
* ubuntu 16.04.6 (on travis-ci), release and devel
* win-builder, release (attempted to run devel, but did not receieve an email with the results despite multiple attempts so ran locally instead - see below)
* local Windows 10 install, devel
* R-hub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit; Ubuntu Linux 16.04 LTS, R-release, GCC; Fedora Linux, R-devel, clang, gfortran; Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE for R-rhub Windows Server 2008 R2 SP1, R-devel, 32/64 bit; and the local Windows 10 devel install:

* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_x64' 'tests_i386' 'tests_x64'
    'bmotif-Ex_i386.Rout' 'bmotif-Ex_x64.Rout' 'examples_i386'

  This seems to be a result of the CMD check rather than an actual problem.
  
There was 1 NOTE for R-hub Ubuntu Linux 16.04 LTS, R-release, GCC:
  
* checking installed package size ... NOTE
  installed size is  6.6Mb
  sub-directories of 1Mb or more:
    libs   4.3Mb
    
And 1 NOTE for R-hub Fedora Linux, R-devel, clang, gfortran:
    
* checking installed package size ... NOTE
  installed size is  5.7Mb
  sub-directories of 1Mb or more:
    libs   3.4Mb
  
It seems that on Linux architectures (though not on Debian Linux, R-devel, GCC ASAN/UBSAN), the check returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OS X. My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp. Indeed, some functions of the bmotif package have been written in C++ using Rcpp. They are needed to perform calculation of motif frequencies. Without the speed up gained from those C++ functions, this package would become impractical.

## Downstream dependencies
There are currently no known downstream dependencies for this package.
