# Network Comparison
An R package implementing the Netdis and NetEMD alignment-free network comparison measures.   


### :warning: BETA: Package under construction (pre-release) :warning:
Until this package hits release 1.0 anything can change with no notice.

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![GitHub release](https://img.shields.io/github/release/alan-turing-institute/network-comparison.svg)](https://github.com/alan-turing-institute/network-comparison/releases/latest)
[![Travis](https://img.shields.io/travis/alan-turing-institute/network-comparison/master.svg)](https://travis-ci.org/alan-turing-institute/network-comparison/branches)
[![Appveyor](https://ci.appveyor.com/api/projects/status/jn1a36c22vjw1l4d/branch/master?svg=true)](https://ci.appveyor.com/project/alan-turing-institute/network-comparison/branch/master)

[![Codecov](https://img.shields.io/codecov/c/github/alan-turing-institute/network-comparison/master.svg)](https://codecov.io/gh/alan-turing-institute/network-comparison?branch=master)
[![license](https://img.shields.io/github/license/alan-turing-institute/network-comparison.svg)](https://github.com/alan-turing-institute/network-comparison/edit/master/LICENSE)
[![Github All Releases](https://img.shields.io/github/downloads/alan-turing-institute/network-comparison/total.svg)](https://github.com/alan-turing-institute/network-comparison/releases/latest)

## Usage
See "Quick start" vignettes in documentation for example usage.

### Parallel processing support
The `gdd_for_all_graphs` method will make use of multiple threads to calculate graphlet-based degree distributions for multiple graphs in parallel. However, this is only supported on unix-like systems (e.g. Linux and Mac OSX), as the underlying `mcapply` method requires system support for forks. The `gdd_for_all_graphs` method will run on Windows, but the number of cores will be restricted to 1, regardless of the value of `mc.cores` provided to `gdd_for_all_graphs` or set in the R environment.

## Installing package from source
When published to the CRAN package repository, the library and all documentation
will be installed in the standard manner using `install.packages("netdist")`,
which will also install all dependencies from CRAN. However, prior to publication
on CRAN, the package can be installed from Github using the `devtools` package.

### Prerequisites
- Install the `devtools` package using `install.packages("devtools")`

### Installing package and core dependencies needed to use it
Core dependencies are listed in the `Imports` section of the package
`DEPENDENCIES` file and will be automatically installed when installing /
updating the package using `devtools` as described below.

Install the latest package from GitHub using:
  - `devtools::install_github("alan-turing-institute/network-comparison")`

### Installing additional optional dependencies
Additional optional packages required to generate documentation, run tests or
develop the package are listed in the `Suggests` section of the package
`DEPENDENCIES`  file. These will not be automatically installed by
`devtools::install()` and, if required, each optional package will need to be
separately installed using `devtools::install("<package_name>")`.

Example optional libraries are:

- `testthat` for running tests via `devtools::test()` or generating package
documentation using `devtools::document()`
- `knitr` for building the package long-form documentation vignettes using
`devtools::build_vignettes()`

## Documentation
You can browse a list of long form examples illustrating how to use the package
using `browseVignettes(package = "netdist")`.

You can list the functions available in the package with `library(help = "netdist")`
and get more detailed help on individual functions using `?function_name` (e.g.
`?net_emd`). In RStudio, typing `?netdist::` should also provide a drop down list
of functions you can select to load the more detailed help.

## References
#### NetEMD
Anatol E. Wegner, Luis Ospina-Forero, Robert E. Gaunt, Charlotte M. Deane, Gesine Reinert; Identifying networks with common organizational principles. [arXiv:1704.00387](https://arxiv.org/abs/1704.00387) [stat.ML] (Open Access pre-print)

#### NetDis
Waqar Ali, Tiago Rito, Gesine Reinert, Fengzhu Sun, Charlotte M. Deane; Alignment-free protein interaction network comparison. Bioinformatics 2014; 30 (17): i430-i437. https://doi.org/10.1093/bioinformatics/btu447 (Open Access)
