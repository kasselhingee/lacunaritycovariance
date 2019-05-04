# TO DO BEFORE PUBLISHING:

change package filename to something narrower: covargbl  ? lacunaritycovariance
- I'm doing this because it is now clearer that I won't be expanding on the capacity of the package in the next few years.
- The more descriptive name here the better.
- A zinger of a name would be nice. lacunaritycovariance?
- But contagion estimators are included?

Remake description file

In centred covariance (and other places) the quotation marks are not in the write direction.

in !im! format --> as an !im! object
sed -i "s/in \\\code{im} format/as an \\\code{im} object/g"
sed -i "s/in \\\code{owin} format/as an \\\code{im} object/g"

delete the plotting things in examples

----

## Really required:

- build using R CMD build --compact-vignettes=gs+qpdf  ...

- More Picka references ACT: 30min

- watch use of GBL when it should be gliding box lacunarity

- update lacunarity paper reference (e.g. it is published..)

- version number with non-leading zeros

- check vignettes and paper source code still run smoothly

Total ACT: 8h 

----

## Great for more polished publication:
- Clean up mix ups between estimates and estimator in language.

- Referencing in coverageprob.R is numeric, whilst most of other pages are using author-year

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

- caps incorrect thing in ccvc_byconv



Total ACT: 37h

----------------------
Completed:
D- unit tests ACT 4h

D- make sure everything works on logically valued images (TRUE, FALSE, NA) created say using eval.im ACT 2h

D- syphon off estimation of cover in small window ACT 1.5h

D- syphon off helper functions never for release ACT 1.5h

D- scdcontag() --> contagdiscstate()? ACT 0.5h

D- setcov and imcov maintain units of input image - pass on to spatstat ACT 1h

D- Change function names to reflect names `empirical GBL' or gliding box estimator, and plugin moment estimators  ACT 2h

D- tradcovarest --> plugincvc

D- 'trad' method name changed to 'plugin'

D- glbtrad description to: empirical Gliding Box, proposed by A and C [1]. And else where ACT 1h

D- tradcovarest require obswin. ACT 1h

D- 'mvlgb or GBLgb to gblemp'

D- second order props has a notyetimplemented argument

D- cross references (links) rather than just code{}

D- change function file names to match function names. (e.g. MVL and tradcovarest)

