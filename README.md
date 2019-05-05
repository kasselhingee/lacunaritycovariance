# lacunaritycovariance 
*Author: Kassel Liam Hingee*

This directory contains the source code of the R package *lacunaritycovariance*. This R package is for estimating gliding box lacunarity and other random closed set properties from binary coverage maps (images composed of binary-valued pixels). 

## Installation
### From GitHub using devtools package:
From inside an R interactive session run:

    library(devtools)
    install_github("kasselhingee/racsstats", ref = "release", auth_token = "****")

where auth_token is required whilst the GitHub repository is private. 
*I will soon make this repository public so the auth_token will not be needed.* 
 If you have a GitHub account then, after getting permission from me (Kassel Hingee), you can obtain an auth_token by going to this website: https://github.com/settings/tokens.
 Copy the personal access token from your browser and use it as the auth_token argument of install_github.

### From .tar.gz file
Inside an R session run

    install.packages("<PATH>", repos = NULL, type = "source")

where `<PATH>` is the path to a .tar.gz file containing the contents of this repository.

### Install from source code not in .tar.gz form
First run the following to build a .tar.gz file.

    R CMD build . 

This should create a file `lacunaritycovariance-<VERSION>.tar.gz`. 
Where `<VERSION>` is a string on numbers separated by periods and hyphens, for example 0.5-2.
Then run 

    R CMD INSTALL lacunaritycovariance-<VERSION>.tar.gz

## Further Notes
### Manual pages
The manual pages for each function (in man/) and the file NAMESPACE have been generated using roxygen2. To edit these, edit the specially formatted comments in the files in the directory R/, then regenerate the manual pages and NAMESPACE function by running `roxygen2::roxygenise()` from an R session with the working directory at the root of the package directory.


### Git repository branches
This source code is the release branch of a git repository. For those familiar with git: it is intended that all developments and modifications of the code are performed in the *working* branch and then merged into this *release* branch.

