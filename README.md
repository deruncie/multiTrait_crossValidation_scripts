Scripts for performing simulations described in paper: [Pitfalls and Remedies for Cross Validation with Multi-trait Genomic Prediction Methods](http://biorxiv.org/cgi/content/short/595397v1). The three scripts:

```
Sim_GP.R
Sim_GP_knownG.R
Sim_GP_big.R
```
all run the same simulation-style. The first is the main set of simulations described in the paper. The middle assumes that `G` and `R` are known, and so don't need to be estimated. The final one uses a larger sample size. 