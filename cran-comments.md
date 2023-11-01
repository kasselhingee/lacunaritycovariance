Attempt to resolve issue with M1 Mac check not finding 'envi' package.
I've removed dependence on envi completely.

The CRAN M1 Mac check also can't find the terra package, but I can't see why that would be the case and my check on using https://mac.r-project.org/macbuilder/results/1698796277-d81b9fffe0030511/ found no issue with terra.

The note in windows checks relates to a URL that has always been in the DESCRIPTION. The URL works on usual browsers.


