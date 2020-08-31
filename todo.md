# TO DO:

- faster running of examples for gblg, gbl

- check handling of na.replace in integration_trad, and elsewhere

- Add estimator selection option to gblg function: function gblg doesn't allow selection of estimator, this is inconsistent with other gbl functions. Futhermore function 'gbl' help suggests it does.

- Make gbl() return a single fv object that contains normed and unnormed gbl estimates. 

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

