# GmAMisc 1.1.1:

* minor fixes to the help documentation.

# GmAMisc 1.1.0:

* minor corrections to the components in the list returned by the `NNa()` function; 
* change in the output of the `perm.t.test()` function, which now produces a frequency distribution histogram; 1-sided permuted p-values are now also reported; user-defined labels for the two samples being tested can be used.
* `distRandSign()`: in the map of from- and to-features, from-features are given a colour according to whether or not they are closer or more distant than expected to the nearest to-feature; users can now 'export' as a shapefile the input dataset featuring 2 new fields: one storing each feature's distance to the nearest to-feature; one containing a string indicating if the corresponding feature is closer or more distant than expected.
* New functions added: `distDiffTest()`, `pointsCovarDistr()`.

# GmAMisc 1.0.0:

first release to CRAN.
