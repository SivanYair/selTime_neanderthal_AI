# selTime_neanderthal_AI
Inferring the time that selection favored Neanderthal introgressed alleles in modern human populations: simulations, statistical inference method, imputation of ancient samples, examples

# Overview
This repository contains code associated with [Yair, Lee, and Coop (2020)](https://doi.org/10.1101/2020.10.04.325183). We provide a brief overview of the contents of each directory below. See those directories and
comments in scripts for more details.

1. **Imputation**: We provide scripts showing how we impute genotypes for ancient samples from genotype likelihoods at fewer sites.
  + [imputation](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/imputation) contains scripts using Beagle 4.1 to impute genotypes and remove sites with low imputation accuracy.
  + [region-info](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/region-info) contains files that list windows of analysis for inference and exploring allele frequencies in putative windows of adaptive introgression.
2. **Simulations and Method**: We provide scripts showing how we simulated data for
    method validation, in addition to how we run the method on those data.
    + [specification_files](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/specification_files) contains scripts that describe the demographic history and sample sizes used for simulations and inference.
    + [selection_simulations](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/selection_simulations) contains scripts to generate simulations under different scenarios of selection. This directory includes files that detail each scenario that we simulated and scripts that obtain linked neutral allele frequencies from tree sequences recorded in SLiM.
    + [method](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/method) contains scripts to run the method on simulations.
