library(devtools)
library(profvis)
library(testthat)

library(devtools)
install(upgrade_dependencies = FALSE)
install(quick = TRUE, upgrade_dependencies = FALSE)

test_check()
getwd()

#timing of testing
profout <- profvis(test("."), interval = 0.1)
#on 31 July the above took 143900ms = 2.4 minutes on my uwa computer.
#After reducing rbdr resolution the time is down to 70 seconds

system.time(test("."))

#running specific tests
profout <- profvis(test(".", filter = "imbinaryops"), interval = 0.1)
htmlwidgets::saveWidget(profout, "testprofile.html")
test(".", filter = "rbdr")

#test plotting with vdiffr:
library(vdiffr)
manage_cases(, filter = "addfvtomap")

check(build_args = list("--compact-vignettes=gs+qpdf"))
