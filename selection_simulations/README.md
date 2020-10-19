This directory contains code to run simulations that record tree sequences, and to
add mutations to and obtain allele frequencies from those tree sequences.

The file [run_slim_and_get_freqs.txt](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/selection_simulations/run_slim_and_get_freqs.txt) shows how to run all of the following steps:

### 1. Make and run simulation scripts
  + [make_selection_slim_script.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/selection_simulations/make_selection_slim_script.R) contains code to write a slim script, given the
     name of the scenario that we aim to simulate. The details of the scenario are
     described in the args_file directory. The number corresponds to a specific combination of selected populations and the time between admixture and the selection onset. The letter corresponds to a specific combination of the selection coefficient and targeted final frequency of the selected allele.  Each scenario has a text file that contains
     lines in the following order describing:
    + the letter referring to the demography version to use (see [../specification_files](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/specification_files))
    + the sample size option to use (see [../specification_files](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/specification_files))
    + the selection coefficient (_s_)
    + the time between admixture and selection (_t<sub>b</sub>_)
    + populations (space-separated) in which the Neanderthal allele is positively selected
  + The slim scripts are written to the [slim_scripts/version{letter}](https://github.com/SivanYair/selTime_neanderthal_AI/tree/main/selection_simulations/slim_scripts/versionA) directory. They need to be run on the command line with defined variable 'trees_file' (path to file that tree sequences are recorded to)


### 2. Get allele frequencies from tree sequences
  + [selection_getFreqs_fromSlimTrees.py](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/selection_simulations/selection_getFreqs_fromSlimTrees.py) contains code to add neutral mutations and
    calculate population allele frequencies according to the sample sizes provided for each population.
    It outputs data frames to feather files.
  + [convert_feather_to_RDS.R](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/selection_simulations/convert_feather_to_RDS.R) contains code to convert the feather files produced in the last step to RDS files

Note that output is provided for Scenario1b's first run. This output is used as input in the [../method](https://github.com/SivanYair/selTime_neanderthal_AI/blob/main/method) directory.
