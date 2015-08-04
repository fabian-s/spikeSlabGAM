`spikeSlabGAM`
===================

[ ![](https://travis-ci.org/fabian-s/spikeSlabGAM.svg?branch=master) ](https://travis-ci.org/fabian-s/spikeSlabGAM)

An R-package for Bayesian variable selection, model choice, and regularized estimation for (spatial) generalized additive mixed regression models via stochastic search variable selection with spike-and-slab priors.

- Fits **additive models** for **Gaussian, Binary/Binomial and Poisson responses** 
- (Correlated) **random effects** 
- Automagically performs **variable selection** and **function selection** (i.e., do I need this effect at all, is it linear or is it non-linear?), also for interactions between multiple covariates. 
- Yields **marginal posterior inclusion probabilities** for each term as well as **posterior model probabilities** and (model-averaged) effect estimates.
- Convenient **formula-based model specification**
- Fully Bayesian via MCMC, multiple parallelized chains for diagnostics and faster mixing, sampler implemented in `C`. 

--------------------------------------------------------------------------------
### References:

Short & applied intro (also the vignette that comes with the package, with some minor modifications...):

> Fabian Scheipl. (2011) `spikeSlabGAM`: Bayesian Variable Selection, Model Choice and Regularization for Generalized Additive Mixed Models in R. *Journal of Statistical Software*, 43(14). [[pdf]](http://www.jstatsoft.org/v43/i14)

More theory, simulation studies and real-world case studies:

> Fabian Scheipl, Ludwig Fahrmeir, Thomas Kneib (2012). Spike-and-Slab Priors for Function Selection in Structured Additive Regression Models. *Journal of the American Statistical Association*, 107(500), 1518-1532. [[pdf on arXiv]](http://arxiv.org/abs/1105.5250)

