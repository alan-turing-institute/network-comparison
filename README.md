# netdist
An R package implementing the NetEMD network comparison measure

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

## Usage
See "Quick start" vignette in documentation for example usage.

## Installing package from source
When published to the CRAN package repository, the library and all documentation
will be installed in the standard manner using `install.packages("netdist")`,
which will also install all dependencies from CRAN. However, prior to publication
on CRAN, some additional effort is required to install the package from source.

### Prerequisites
- Install the `devtools` package using `install.packages("devtools")`

### Installing package and core dependencies needed to use it
Core dependencies are listed in the `Imports` section of the package 
`DEPENDENDCIES` file and will be automatically installed when installing / 
updating the package using `devtools` as described below. Note that this does
not build the package documentation.

  - Fetch the package source from GitHub 
  ("https://github.com/alan-turing-institute/network-comparison.git")
  - Open an R terminal in the folder containing the package source
  - Install the package using `devtools::install()`
  
### Installing additional optional dependencies
Additional optional packages required to generate documentation, run tests or 
develop the package are listed in the `Suggests` section of the apckage 
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
