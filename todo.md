# TO DO BEFORE PUBLISHING:

----

## Really required:

- build using R CMD build --compact-vignettes=gs+qpdf  ...

- Change function names to reflect names `empirical GBL' or gliding box estimator, and plugin moment estimators  ACT 2h

- glbtrad description to: empirical Gliding Box, proposed by A and C [1]. And else where ACT 1h

- consistent naming ACT 1h

- More Picka references ACT: 30min

- tradcovarest require obswin. ACT 1h

- change function file names!

- version number with non-leading zeros

- check vignettes and paper source code still run smoothly

Total ACT: 8h 

----

## Great for more polished publication:
- In pixel contag, what is the `double count' method as described in the FRAGSTATS manual on adjacency matrix, and is it different to what is in this package? 

- better indexing of topics. ACT 1h

- tests of summary properties of grain.lib. ACT 3h

- move covar.grainlib function to another help file  ACT 2h

- disc state contagion function that acts on binary map directly. ACT 2h

- Discussion on when to use the functions? Take 50% from thesis intro on pros and cons of summary functions. ACT 3h (to first smooth draft)

- Make it make sense to someone who isn't a spatial stats expert, to an ecologist, that I want to use it. ACT 3h (to first smooth draft)

##**feature request: a test that simulations rbllognorm correspond to theoretical values for the model. ACT 1h

##** add mean area contagion. And contagion for categorical maps. ACT 20h

- polish vignette? ACT 2h (to first smooth draft)


Total ACT: 37h

----------------------
Completed:
D- unit tests ACT 4h

D- make sure everything works on logically valued images (TRUE, FALSE, NA) created say using eval.im ACT 2h

D- syphon off estimation of cover in small window ACT 1.5h

D- syphon off helper functions never for release ACT 1.5h

D- scdcontag() --> contagdiscstate()? ACT 0.5h

D- setcov and imcov maintain units of input image - pass on to spatstat ACT 1h
