# Contributing to lacunaritycovariance

Thank you for helping to improve lacunaritycovariance. Please feel free to email me through kassel.hingee@gmail.com. This is the first open source project I've coordinated so I anticipate many questions not covered by the documentation.

#### Table of Contents
[Files that are machine generated](#files-that-are-machine-generated)

[Testing package](#testing-package)

[Code style conventions](#code-style-conventions)

[To-Do File](#to-do-file)

## Files that are machine generated.
The package contains a number of files that are machine generated:
  + 'DESCRIPTION' (generated from 'DESCRIPTION_editable'). Please only edit the file 'DESCRIPTION_editable', and then rebuild 'DESCRIPTION'. See below. This system is in place so that small changes in the description of the package don't require heavy manual reformatting of the text flow.
  + 'NAMESPACE'  (by 'roxygen2')
  + All files in 'man/' (generated by roxygen2). Do no edit these files directly, instead edit the corresponding files in 'R/' and then apply 'roxygen2'. See below.

To create all the machine-generated files run from command line:

    R CMD BATCH makeR.R

The R script, makeR.R, combines the following steps for building the DESCRIPTION file and help files:
  + To generate 'DESCRIPTION' run the following in 'R'
(in the future this step should be combined with other package building steps):

        library(desc)
        desc <- description$new("DESCRIPTION_editable")
        desc2 <- desc$normalize()
        desc2$write(file = "DESCRIPTION")

  + To generate the 'NAMESPACE' file and all files in 'man/' (except 'lacunaritycovariance-package.Rd'), within an R session run:

        library(roxygen2)
        roxygenise()

## Testing package
The package is tested using the package 'testthat'. Tests can be found in the directory 'tests/testthat/'. Tests can be run with the following in R:

    library(devtools)
    test()

For quick checking or checking on CRAN many tests use smaller images, or are skipped. Such tests are most reliably invoked using 'R CMD check'. Setting the R environment variables 'NOT_CRAN' to 'false' is not sufficient as devtools appears to run the test helpers at full resolution regardless (devtools might be setting temporarily setting NOT_CRAN to 'true' for this part?).

## Code style conventions
I used the 'lintr' package  to check the styling of my code. lintr settings are below. Please leave spaces beside binary operators and after commas.

'xi' or 'xiim' should always be used to mean a binary map or binary image.

For function names there is not a strong convention yet. However they somewhat follow these principles:
  + no camel case
  + acronyms/abbreviations are lower case
  + implementations of algorithms/operations that have a capital in their name may have corresponding capitals in the function name
  + periods in function names are ok

The lintr settings used for this package are:

    lacunaritycovariance_linters <- list(
      assignment_linter = assignment_linter,
      single_quotes_linter = single_quotes_linter,
      absolute_paths_linter = absolute_paths_linter,
      no_tab_linter = no_tab_linter,
      line_length_linter = line_length_linter,
      commas_linter = commas_linter,
      infix_spaces_linter = infix_spaces_linter,
      spaces_left_parentheses_linter = spaces_left_parentheses_linter,
      spaces_inside_linter = spaces_inside_linter,
      open_curly_linter = open_curly_linter,
      closed_curly_linter = closed_curly_linter(allow_single_line = TRUE),
      camel_case_linter = camel_case_linter,
      snake_case_linter = snake_case_linter,
      multiple_dots_linter = multiple_dots_linter,
      object_length_linter = object_length_linter,
      object_usage_linter = object_usage_linter,
      trailing_whitespace_linter = trailing_whitespace_linter,
      trailing_blank_lines_linter = trailing_blank_lines_linter,
      commented_code_linter = commented_code_linter
    )

These linters can be used to examine a file from within R by:

    library(lintr)
    lint(<filename>, linters = lacunaritycovariance_linters)

## To-Do File
The file 'todo.md' is a list of various things I'd like to do to improve the package. Mostly improving readability and filling missing parts of the package's functionality. Other contributors such as you are not expected to follow this list, though I won't complain if you do!
