# sppmix
Mixture models for spatial point patterns

The sppmix package implements classes and methods for modeling spatial point process data using Poisson point processes where the intensity surface is assumed to be a multiple of a finite additive mixture of normal components. 

We utilize the ppp and owin classes from the spatstat package in order to describe a point pattern. We further introduce a normmix class for handling 2d mixtures of bivariate normal components, and help us build the Poisson point process intensity surface.

Data augmentation MCMC (DAMCMC) and Birth-Death MCMC (BDMCMC) are the two main methods we have implemented for estimating the Poisson intensity surface, in a Bayesian framework.

The MCMC algorithms are implemented in C++ using Rcpp and RcppArmadillo, and were optimized after extensive testing, meaning that this approach is significantly faster than some other implementations of mixture models.

Plotting is accomplished using the rgl package in order to create 3d plots of the intensity surfaces. In addition, the fields and ggplot packages were used for 2d plots.

To learn more about sppmix, start with the basic tutorials at http://www.stat.missouri.edu/~amicheas/sppmix/sppmix_tutorial_links.html
