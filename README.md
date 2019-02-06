# stationaryracsinference 

This is an R package for estimating gliding box lacunarity and other RACS properties. This directory contains the source code of this package.

## Installation
Install from .tar.gz form of this directory: Inside an R session run
    install.packages("<PATH>", repos = NULL, type = "source")

where <PATH> is the path to the .tar.gz file.

*or*

To install from this source code in uncompressed form:
    R CMD build . 
This should create a file "stationaryracsinference-***.tar.gz"
Then run 
    R CMD INSTALL stationaryracsinference-***.tar.gz

*or* 

Install from GitHub using devtools package:
    library(devtools)
    install_github("kasselhingee/racsstats", ref = "release", auth_token = "****")

where auth_token is required until the GitHub repository is made public. If you have a GitHub account then, after getting permission from Kassel Hingee, you can obtain an auth_token by going to this website: https://github.com/settings/tokens.
 Copy the personal access token from your browser and use it as the auth_token argument of install_github.

## Further Notes
### Manual pages
The manual pages for each function (in man/) and the file NAMESPACE have been generated using roxygen2. To edit these, edit the specially formatted comments in the files in the directory R/, then regenerate the manual pages and NAMESPACE function by running roxygen2::roxygenise() from an R session with the working directory at the root of the package directory.


### Git repository branches
This source code is the release branch of a git repository. For those familiar with git, it is intended that all developments and modifications of the code are performed in the *working* branch and then merged into this *release* branch.
