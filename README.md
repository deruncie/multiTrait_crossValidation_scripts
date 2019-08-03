Scripts for performing simulations described in paper: [Pitfalls and Remedies for Cross Validation with Multi-trait Genomic Prediction Methods](http://biorxiv.org/cgi/content/short/595397v1). 

The scripts `test_*.R` run different versions of the simulations:

1. known G,R:
   a. simulate data, hold out 10%, calculate cor(u,u\_hat), cor(y,u\_hat) analytically, for single, CV1, CV2
   b. simulate data, hold out 10%, simulate cor(u,u\_hat), cor(y,u\_hat) for single, CV1, CV2
2. baseline
   a. simulate data, hold out 10%, estimate cor(u,u\_hat), cor(y,u\_hat) for single, CV1, CV2
3. parametric
   a. simulate data, hold out 50%, estimate cor\_g(y,y\_hat)/sqrt(h2_y_hat)
4. semi-parametric
   a. simulate data, hold out 10%, calculate cor(u,u\_hat), cor(y,u\_hat) for single, CV1, CV2 with correction
5. non-parametric
   a. simulate data, hold out 10%, split 2 clones, estimate cor(u,u\_hat), cor(y,u\_hat) for single, CV1, CV2, each with 2nd clone. Adjust residual variance
   b. simulate data, hold out 20% (10% validation + nearest relatives), estimate cor(u,u\_hat), cor(y,u_hat) for single, CV1, CV2, each with relative
   
Data from Lopez-Cruz, M., J. Crossa, D. Bonnett, S. Dreisigacker, J. Poland,
et al., 2015 Increased Prediction Accuracy in Wheat Breeding Trials Using a Marker × Environment Interaction Genomic Se- lection Model. G3: Genes|Genomes|Genetics 5: 569–582. is included in `FleS3.rdata`. It was downloaded directly from the paper website.