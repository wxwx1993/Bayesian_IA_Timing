# optimIA
## The R Package for the manuscript "Optimizing Interim Analysis Timing for Bayesian Adaptive Commensurate Designs". 
An innovative Bayesian adaptive commensurate design that is particularly useful for rare disease trials since it has the capability to borrow adaptively from historical information and also uses a payoff function to optimize the timing of the study’s interim analysis. (https://arxiv.org/abs/1905.07456).

## Installation
```r
library("devtools")
install_github("wxwx1993/Bayesian_IA_Timing")
library("optimIA")
```

## Example
```r
optimIA(historical=adult, 
        w=0.5, 
        IA_num=seq(0,40,2), 
        N=40, 
        p_l=0.25, 
        p_0=0.975, 
        p_u=c(0.96,0.998,0.002), 
        delta_0=0, 
        delta_1=20, 
        delta_design=25, 
        delta_min=15, 
        type1_level=0.05,  
        sim_sd=25, 
        niter=5000, 
        nrep=5000, 
        mc.cores=6)
```
        
## Code
BCDOT_function is functions for implementing the Bayesian algorithm for the commensurate design approach to optimize interim analysis timing.

commensurate.stan & commensurate_noborrow.stan is to conduct essential Bayesian analyses for the algorithm approch. The commensurate_noborrow.stan is aimming to calculate Effective Historical Sample Size.

## Data
adult is the fixed 50 sample size (hypothetical) historical adult data with balanced arms used for analysis.
adult20 is the fixed 20 sample size (hypothetical) historical adult data with balanced arms used for analysis.
adult100 is the fixed 100 sample size (hypothetical) historical adult data with balanced arms used for analysis.

## Reference
Wu, X., Xu, Y. and Carlin, B.P., 2020. Optimizing interim analysis timing for Bayesian adaptive commensurate designs. Statistics in Medicine, 39(4), pp.424-437.

