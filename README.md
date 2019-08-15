# lacunaritycovariance 
*Author: Kassel Liam Hingee*

This directory contains the source code of the R package *lacunaritycovariance*. This R package is for estimating gliding box lacunarity and other random closed set properties from binary coverage maps (images composed of binary-valued pixels). 

#### Table of Contents
[Installation](#Installation)
  + [From GitHub using remotes package](#from-github-using-remotes-package)
  + [From .tar.gz file](#from-.tar.gz-file)
  + [Install from source code not in .tar.gz form](#install-from-source-code-not-in-.tar.gz-form)

[Manual pages](#manual-pages)

[Git repository branches](#git-repository-branches)

## Installation
### From GitHub using remotes package
From inside an R interactive session run:

    library(remotes)
    install_github("kasselhingee/lacunaritycovariance", ref = "release")

### From .tar.gz file
Inside an R session run

    install.packages("<PATH>", repos = NULL, type = "source")

where `<PATH>` is the path to a .tar.gz file containing the contents of this repository.

### Install from source code not in .tar.gz form
First run the following to build a .tar.gz file.

    R CMD build --compact-vignettes=gs+qpdf . 

This should create a file `lacunaritycovariance-<VERSION>.tar.gz`. 
Where `<VERSION>` is a string on numbers separated by periods and hyphens, for example 1.0-0.
Then run 

    R CMD INSTALL lacunaritycovariance-<VERSION>.tar.gz

## Manual pages
The manual pages for each function (in man/) and the file NAMESPACE have been generated using roxygen2. To edit these, edit the specially formatted comments in the files in the directory R/, then regenerate the manual pages and NAMESPACE function by running `roxygen2::roxygenise()` from an R session with the working directory at the root of the package directory.


## Git repository branches
It is intended that the *release* branch is consistent with the package available on CRAN (except for brief periods whilst the new versions are submitted to CRAN and awaiting approval). Editing of the package should occur in either the *master* branch or other branches. 

