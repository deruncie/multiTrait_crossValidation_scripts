Scripts for performing simulations described in paper: [Pitfalls and Remedies for Cross Validation with Multi-trait Genomic Prediction Methods](). The three scripts:

```
Sim_GP.R
Sim_GP_knownG.R
Sim_GP_big.R
```
all run the same simulation-style. The first is the main set of simulations described in the paper. The middle assumes that `G` and `R` are known, and so don't need to be estimated. The final one uses a larger sample size. The first and last use a modified version of the `relmatLmer` function from the [lme4qtl](https://github.com/variani/lme4qtl) package. The modification is necessary to use correlated random effects. The modified function is loaded by `source(my_relmatLmer.R)`.