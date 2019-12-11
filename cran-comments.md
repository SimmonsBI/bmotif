## Test environments
* local OS X install, R 3.6.1
* ubuntu 16.04.6 (on travis-ci), release and devel
* win-builder, release
* rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit; Ubuntu Linux 16.04 LTS, R-release, GCC; Fedora Linux, R-devel, clang, gfortran; Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE for R-rhub Windows Server 2008 R2 SP1, R-devel, 32/64 bit:

* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_x64' 'tests_i386' 'tests_x64'
    'bmotif-Ex_i386.Rout' 'bmotif-Ex_x64.Rout' 'examples_i386'

  This seems to be a result of the R-hub check rather than an actual problem.

## Downstream dependencies
There are currently no known downstream dependencies for this package.
