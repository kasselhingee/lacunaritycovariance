library(devtools)
check(env_vars = list(NOT_CRAN = TRUE))
check_win_release()
check_win_devel()
build()
