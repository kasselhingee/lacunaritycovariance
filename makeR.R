library(desc)
desc <- description$new("DESCRIPTION_editable")
desc2 <- desc$normalize()
desc2$write(file = "DESCRIPTION")

library(roxygen2)
roxygenise()

