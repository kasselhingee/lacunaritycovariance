
# INSTALLATION
This package uses roxygen2 and desc to generate R manual files, and the package DESCRIPTION file.
You can install the release version as is using devtools (see README).
To install this development version you will need roxygen2 and desc installed, and then run the following commands in R:

## TO BUILD DESCRIPTION FILE FROM DESCRIPTION_editable:
library(desc)
desc <- description$new("DESCRIPTION_editable")
desc2 <- desc$normalize()
desc2$write(file = "DESCRIPTION")

## TO BUILD MAN AND NAMESPACE FILES
run roxygen2::roxygenise() in R from the root directory of this package.

The package can now be installed using R CMD INSTALL or devtools' install().
