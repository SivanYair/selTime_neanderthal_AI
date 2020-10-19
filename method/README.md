This folder contains all of the scripts necessary to run the method on a particular
region and partition site. The scripts should work with the example simulated data
provided in the [selection_simulations](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/selection_simulations) folder.

For a detailed introduction of the steps taken to run the method,
see this [dmc example](https://github.com/kristinmlee/dmc/blob/master/dmc_example.md) from Lee and Coop (2017).
Our approach here differs slightly, but uses the same general logic.

# Overview of files
 1. **[calcNeutralF.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/calcNeutralF.R)**: calculates neutral **F** matrix used in the method from genome-wide neutral allele frequency data. We provide the matrix for demography version A and sample size option 1 in [neutralF_versionA_ss1.RDS](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/neutralF_versionA_ss1.RDS) (version refers to demography, ss refers to sample size option).
 2. **[master_script.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/master_script.R)**: This script is called from the command line. It loads in the raw data, and sources all other scripts, which are listed below in the order that they are sourced.
 3. **[prep_data_perRun.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/prep_data_perRun.R)**: filters allele frequency data, defines admixture graph parameters, and changes neutral **F** matrix to represent partitioned (rather than whole) populations.
 4. **[genSelMatrices_fxns.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/genSelMatrices_fxns.R)**: defines functions to calculate predicted **F<sup>(N)</sup>** or **F<sup>(S)</sup>** matrix for each of the parameters they are a function of
 5. **[genSelMatrices_exec.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/genSelMatrices_exec.R)**: uses above defined functions to predict **F<sup>(N)</sup>** under
    different genetic distances from the selected site (_r_), and to predict **F<sup>(S)</sup>** under
    different times between admixture and selection (_t<sub>b</sub>_), selection coefficients (_s_),
    and genetic distances from the selected site (_r_).
 6. **[calcCompLike_fxns.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/calcCompLike_fxns.R)**: defines functions to calculate the composite likelihood
    of a given combination of parameters of each function **F<sup>(N)</sup>** and **F<sup>(S)</sup>**.
 7. **[calcCompLike_exec.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/calcCompLike_exec.R)**: uses above defined functions to calculate composite likelihoods and saves these results to a file.
 8. **[plotCompLikeResults.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method/plotCompLikeResults.R)**: plots results and prints estimates of _t<sub>b</sub>_ and _s_.

# Details to run the method
You can run the method on the data from the first run of the Scenario1b simulation as follows:

```
Rscript --vanilla master_script.R ../selection_simulations/args_files/Scenario1b.txt 1 1
```

The first argument is the path to the file detailing the scenario that was simulated, the second argument is the sample size option, and the third argument is the run of the simulated scenario.
