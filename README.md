# Bayesian_IA_Timing
Code for "Optimizing Interim Analysis Timing for Bayesian Adaptive Commensurate Designs". A innovative Bayesian adaptive commensurate designs that borrows adaptively from historical information and also uses a payoff function to optimize the timing of the study’s interim analysis.

Code:
BCDOT_function is functions for implementing the Bayesian algorithm for the commensurate design approach to optimize interim analysis timing.
commensurate.stan & commensurate_noborrow.stan is to conduct essential Bayesian analyses for the algorithm approch. The commensurate_noborrow.stan is aimming to calculate Effective Historical Sample Size.

Data:
adult.RData is the fix (hypothetical) historical adult data used for analysis.
