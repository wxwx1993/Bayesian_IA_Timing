# optimIA
The R Package for the manuscript "Optimizing Interim Analysis Timing for Bayesian Adaptive Commensurate Designs". (https://arxiv.org/abs/1905.07456). An innovative Bayesian adaptive commensurate design that borrows adaptively from 
historical information and also uses a payoff function to optimize the timing of the study’s interim analysis.

Code: 

BCDOT_function is functions for implementing the Bayesian algorithm for the commensurate design approach to optimize interim analysis timing.

commensurate.stan & commensurate_noborrow.stan is to conduct essential Bayesian analyses for the algorithm approch. The commensurate_noborrow.stan is aimming to calculate Effective Historical Sample Size.

Data:

adult is the fixed 50 sample size (hypothetical) historical adult data with balanced arms used for analysis.
adult20 is the fixed 20 sample size (hypothetical) historical adult data with balanced arms used for analysis.
adult100 is the fixed 100 sample size (hypothetical) historical adult data with balanced arms used for analysis.


